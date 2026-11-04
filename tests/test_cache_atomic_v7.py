"""r406 PROTOTYPE tests (isolated copy only): CEA cache atomic replacement + write errors, v5 (v4 tests kept;
round-57 regressions at the end).  Originally: v4 (cache-only diff on
da4e95d).  v3 behaviour that is kept, plus regressions for the round-56 review:
xattrs (user.*) kept or refused non-destructively; directory open / file fsync / directory fsync / close classified
separately (only ENOTSUP / EOPNOTSUPP from fsync(dir_fd) is the warning; EINVAL is exit 3); fd-based attributes;
temp/target name checks against a symlink swap; truthful cleanup report."""
import errno
import hashlib
import importlib.util
import json
import os
import shutil
import stat
import subprocess
import sys
import threading
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BODY = sorted(ROOT.glob("GRIBS_v*.py"))[0]
G3 = ROOT / "cases" / "golden" / "g3_low_condensed.json"
spec = importlib.util.spec_from_file_location("gribs_cache_v7", BODY)
g = importlib.util.module_from_spec(spec); sys.modules[spec.name] = g; spec.loader.exec_module(g)
ACL = "system.posix_acl_access"
SUPP = [x for x in os.getgroups() if x != os.getgid()]
needs_acl = pytest.mark.skipif(shutil.which("setfacl") is None, reason="setfacl not available")
REAL_OPEN, REAL_FSYNC, REAL_UNLINK, REAL_SETXATTR, REAL_RENAME, REAL_CLOSE = (
    os.open, os.fsync, os.unlink, os.setxattr, os.rename, os.close)


def cfg(tmp, name, out, cache=None, rebuild=True):
    d = json.loads(G3.read_text()); d["output"]["directory"] = str(out)
    d["thermochemistry"]["cea_python"].update(cache_directory=str(cache) if cache else "", rebuild_cache=rebuild)
    p = tmp / f"{name}.json"; p.write_text(json.dumps(d)); return p


def cli(c):
    e = dict(os.environ); e.pop("GRIBS_PROTO_M3", None)
    return subprocess.run([sys.executable, str(BODY), "--config", str(c)], capture_output=True, text=True, env=e, timeout=900)


def rec(stderr):
    for line in stderr.splitlines():
        if line.startswith("{"):
            try:
                r = json.loads(line)
            except ValueError:
                continue
            if "code" in r:
                return r
    return None


def xall(p):
    return {k: os.getxattr(str(p), k) for k in os.listxattr(str(p))}


def attrs(p):
    st = os.stat(p); return (stat.S_IMODE(st.st_mode), st.st_gid, xall(p))


def existing(tmp_path, mode=0o640, content="old"):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"; f.write_text(content); os.chmod(f, mode); return f


def write(f, text="new"):
    return g._cache_atomic_write_text(f, text)


def failed(e):
    return g._cache_write_failed(e)


@pytest.fixture
def checkpoint(monkeypatch):
    """install callbacks per stage: checkpoint[stage] = fn(tmp_path)"""
    calls = {}

    def cp(stage, dfd, tmp):
        fn = calls.get(stage)
        if fn:
            fn(Path(os.readlink(f"/proc/self/fd/{dfd}")) / tmp)
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)
    return calls


# ===================================================================== CLI: format, errors (kept from v3)
@pytest.mark.integration
def test_default_run_format_unchanged(tmp_path):
    out = tmp_path / "out"
    r = cli(cfg(tmp_path, "a", out))
    assert r.returncode == 0, r.stderr[-1500:]
    s = json.loads((out / "summary.json").read_text())
    assert "run_termination" not in s["results"] and "run_identity" not in s
    assert "W_THERMO_CACHE_DURABILITY_NOT_GUARANTEED" not in [w["code"] for w in s["results"]["warnings"]]
    assert sorted(p.name for p in out.iterdir()) == sorted(
        ["summary.json", "summary.txt", "time_history.csv", "ballistics.png", "thermo_cache"])
    assert [p.name for p in (out / "thermo_cache").iterdir() if p.name.endswith(".tmp")] == []


def _ref(tmp_path):
    c = tmp_path / "refcache"; assert cli(cfg(tmp_path, "ref", tmp_path / "refout", c)).returncode == 0
    return sorted(c.glob("thermo_*.json"))[0]


@pytest.mark.integration
def test_read_only_file_exit3_before_replace(tmp_path):
    f = _ref(tmp_path); os.chmod(f, 0o444); h = hashlib.sha256(f.read_bytes()).hexdigest()
    r = cli(cfg(tmp_path, "w", tmp_path / "o", f.parent)); x = rec(r.stderr)
    assert r.returncode == 3 and x["code"] == "E_THERMO_CACHE_WRITE_FAILED"
    assert x["context"]["phase"] == "before_replace" and x["context"]["step"] == "inspect_existing"
    assert x["context"]["replaced"] is False and x["context"]["temp_cleanup"]["status"] == "not_created"
    assert hashlib.sha256(f.read_bytes()).hexdigest() == h and sorted(p.name for p in f.parent.iterdir()) == [f.name]


@pytest.mark.integration
def test_read_only_directory_exit3_no_fallback(tmp_path):
    f = _ref(tmp_path); ino = f.stat().st_ino; os.chmod(f.parent, 0o555)
    try:
        r = cli(cfg(tmp_path, "w", tmp_path / "o", f.parent)); x = rec(r.stderr)
        assert r.returncode == 3 and x["context"]["phase"] == "before_replace"
        assert x["context"]["step"] == "create_temp" and x["context"]["errno"] == errno.EACCES
        assert f.stat().st_ino == ino
    finally:
        os.chmod(f.parent, 0o755)


@pytest.mark.integration
def test_cache_directory_cannot_be_created(tmp_path):
    ro = tmp_path / "ro"; ro.mkdir(); os.chmod(ro, 0o555)
    try:
        r = cli(cfg(tmp_path, "w", tmp_path / "o", ro / "sub")); x = rec(r.stderr)
        assert r.returncode == 3 and x["context"]["phase"] == "create_directory"
    finally:
        os.chmod(ro, 0o755)


# ===================================================================== review 2: operation stages separated
def _is_dir_fd(fd):
    return stat.S_ISDIR(os.fstat(fd).st_mode)


@pytest.mark.integration
@pytest.mark.parametrize("eno", [errno.ENOTSUP, getattr(errno, "EOPNOTSUPP", errno.ENOTSUP)])
def test_directory_fsync_unsupported_warns_and_continues(tmp_path, monkeypatch, capsys, eno):
    def fs(fd):
        if _is_dir_fd(fd):
            raise OSError(eno, os.strerror(eno))
        return REAL_FSYNC(fd)
    monkeypatch.setattr(os, "fsync", fs)
    out = tmp_path / "o"; cache = tmp_path / "c"
    rc = g.main(["--config", str(cfg(tmp_path, "w", out, cache))])
    err = capsys.readouterr().err
    assert rc == 0 and "durability after a crash or power loss is not guaranteed" in err
    w = [w for w in json.loads((out / "summary.json").read_text())["results"]["warnings"]
         if w["code"] == "W_THERMO_CACHE_DURABILITY_NOT_GUARANTEED"]
    assert len(w) == 1 and w[0]["context"]["errno"] == eno and w[0]["context"]["replaced"] is True


@pytest.mark.integration
@pytest.mark.parametrize("eno", [errno.EINVAL, errno.EIO, errno.EBADF])
def test_directory_fsync_einval_and_real_failures_exit3(tmp_path, monkeypatch, capsys, eno):
    def fs(fd):
        if _is_dir_fd(fd):
            raise OSError(eno, os.strerror(eno))
        return REAL_FSYNC(fd)
    monkeypatch.setattr(os, "fsync", fs)
    rc = g.main(["--config", str(cfg(tmp_path, "w", tmp_path / "o", tmp_path / "c"))])
    x = rec(capsys.readouterr().err)
    assert rc == 3 and x["code"] == "E_THERMO_CACHE_WRITE_FAILED"
    assert x["context"]["phase"] == "after_replace" and x["context"]["step"] == "fsync_directory"
    assert x["context"]["replaced"] is True and "previous cache is gone" in x["message"]
    assert not (tmp_path / "o" / "summary.json").exists()


@pytest.mark.fast
@pytest.mark.parametrize("eno", [errno.EINVAL, errno.ENOTSUP, errno.EACCES])
def test_directory_open_failure_is_not_the_warning(tmp_path, monkeypatch, eno):
    f = existing(tmp_path); d = str(f.parent)

    def op(path, flags, mode=0o777, *a, dir_fd=None, **k):
        if flags & os.O_DIRECTORY and str(path) == d:
            raise OSError(eno, "simulated directory open failure")
        return REAL_OPEN(path, flags, mode, *a, dir_fd=dir_fd, **k)
    monkeypatch.setattr(os, "open", op)
    with pytest.raises(g._CacheWriteError) as e:
        write(f)
    assert (e.value.phase, e.value.step, e.value.replaced) == ("before_replace", "open_directory", False)
    assert e.value.cleanup == {"status": "not_created", "created": False} and f.read_text() == "old"


@pytest.mark.fast
@pytest.mark.parametrize("eno", [errno.ENOTSUP, errno.EINVAL, errno.EIO])
def test_file_fsync_failure_is_not_the_warning(tmp_path, monkeypatch, eno):
    f = existing(tmp_path)

    def fs(fd):
        if not _is_dir_fd(fd):
            raise OSError(eno, os.strerror(eno))
        return REAL_FSYNC(fd)
    monkeypatch.setattr(os, "fsync", fs)
    with pytest.raises(g._CacheWriteError) as e:
        write(f)
    assert (e.value.phase, e.value.step, e.value.replaced) == ("before_replace", "fsync_file", False)
    assert e.value.cleanup["status"] == "removed" and f.read_text() == "old"
    assert sorted(p.name for p in f.parent.iterdir()) == ["x.json"]


@pytest.mark.fast
def test_directory_close_failure_after_replace(tmp_path, monkeypatch):
    f = existing(tmp_path)

    def cl(fd):
        dir_fd = _is_dir_fd(fd)
        REAL_CLOSE(fd)
        if dir_fd:
            raise OSError(errno.EIO, "simulated close failure")
    monkeypatch.setattr(os, "close", cl)
    with pytest.raises(g._CacheWriteError) as e:
        write(f)
    assert (e.value.phase, e.value.step, e.value.replaced) == ("after_replace", "close_directory", True)
    assert f.read_text() == "new"


# ===================================================================== review 1: extended attributes
@pytest.mark.fast
def test_user_xattrs_kept(tmp_path):
    f = existing(tmp_path)
    os.setxattr(f, "user.gribs.note", b"keep me"); os.setxattr(f, "user.bin", b"\x00\x01\xff")
    before = xall(f)
    assert write(f) is None
    assert f.read_text() == "new" and xall(f) == before


@pytest.mark.fast
def test_xattr_that_cannot_be_reproduced_refuses_non_destructively(tmp_path, monkeypatch):
    """security.* cannot be created without privilege here; simulated by refusing to set one user.* name."""
    f = existing(tmp_path); os.setxattr(f, "user.protected", b"v"); before = (f.read_bytes(), attrs(f), f.stat().st_ino)

    def sx(target, name, value, *a, **k):
        if name == "user.protected":
            raise PermissionError(errno.EPERM, "simulated: attribute cannot be set")
        return REAL_SETXATTR(target, name, value, *a, **k)
    monkeypatch.setattr(os, "setxattr", sx)
    with pytest.raises(g._CacheWriteError) as e:
        write(f)
    assert (e.value.phase, e.value.step) == ("before_replace", "set_attributes")
    assert e.value.cleanup["status"] == "removed"
    assert (f.read_bytes(), attrs(f), f.stat().st_ino) == before
    assert sorted(p.name for p in f.parent.iterdir()) == ["x.json"]
    assert "existing cache, if any, was not modified" in str(failed(e.value))


@pytest.mark.fast
def test_silently_dropped_xattr_is_detected(tmp_path, monkeypatch):
    """setxattr 'succeeds' but the attribute is not there -> verify_attributes refuses."""
    f = existing(tmp_path); os.setxattr(f, "user.ghost", b"v")

    def sx(target, name, value, *a, **k):
        if name == "user.ghost":
            return None
        return REAL_SETXATTR(target, name, value, *a, **k)
    monkeypatch.setattr(os, "setxattr", sx)
    with pytest.raises(g._CacheWriteError) as e:
        write(f)
    assert e.value.step == "verify_attributes" and f.read_text() == "old" and xall(f) == {"user.ghost": b"v"}


# ===================================================================== review 3: fd-based attributes, name checks
@pytest.mark.fast
def test_symlink_swap_before_rename_is_refused(tmp_path, checkpoint):
    f = existing(tmp_path); victim = tmp_path / "victim.json"; victim.write_text("NOT A CACHE"); os.chmod(victim, 0o600)

    def swap(tmp):
        REAL_UNLINK(tmp); os.symlink(victim, tmp)
    checkpoint["before_replace"] = swap
    with pytest.raises(g._CacheWriteError) as e:
        write(f)
    assert (e.value.phase, e.value.step) == ("before_replace", "verify_temp_name")
    assert e.value.cleanup["status"] == "left_untouched_not_this_runs_file"
    assert not f.is_symlink() and f.read_text() == "old"
    assert victim.read_text() == "NOT A CACHE" and stat.S_IMODE(victim.stat().st_mode) == 0o600


@pytest.mark.fast
def test_symlink_swap_before_attributes_does_not_touch_the_other_file(tmp_path, checkpoint):
    """fd-based fchmod/fchown/xattr act on this run's file, not on what the name points to."""
    f = existing(tmp_path, 0o640); os.setxattr(f, "user.a", b"1")
    victim = tmp_path / "victim.json"; victim.write_text("V"); os.chmod(victim, 0o600)

    def swap(tmp):
        os.rename(tmp, tmp.with_name(tmp.name + ".moved")); os.symlink(victim, tmp)
    checkpoint["temp_created"] = swap
    with pytest.raises(g._CacheWriteError) as e:
        write(f)
    assert e.value.step == "verify_temp_name"
    assert stat.S_IMODE(victim.stat().st_mode) == 0o600 and xall(victim) == {} and victim.read_text() == "V"
    assert f.read_text() == "old"
    moved = [p for p in f.parent.iterdir() if p.name.endswith(".moved")]
    assert len(moved) == 1 and stat.S_IMODE(moved[0].stat().st_mode) == 0o640 and xall(moved[0]) == {"user.a": b"1"}


@pytest.mark.fast
def test_swap_after_rename_is_reported_not_success(tmp_path, monkeypatch):
    f = existing(tmp_path); victim = tmp_path / "victim.json"; victim.write_text("NOT A CACHE")

    def rn(src, dst, *a, **k):
        REAL_RENAME(src, dst, *a, **k)
        REAL_UNLINK(dst, dir_fd=k.get("dst_dir_fd")); os.symlink(victim, dst, dir_fd=k.get("dst_dir_fd"))
    monkeypatch.setattr(os, "rename", rn)
    with pytest.raises(g._CacheWriteError) as e:
        write(f)
    assert (e.value.phase, e.value.step, e.value.replaced) == ("after_replace", "verify_target", True)


# ===================================================================== review 4: truthful cleanup report
@pytest.mark.fast
def test_cleanup_failure_is_reported_and_primary_error_kept(tmp_path, checkpoint, monkeypatch):
    f = existing(tmp_path)

    def boom(tmp):
        raise OSError(errno.EIO, "primary failure")
    checkpoint["before_replace"] = boom

    def ul(path, *a, **k):
        if str(path).endswith(".tmp"):
            raise PermissionError(errno.EACCES, "cannot remove")
        return REAL_UNLINK(path, *a, **k)
    monkeypatch.setattr(os, "unlink", ul)
    with pytest.raises(g._CacheWriteError) as e:
        write(f)
    te = failed(e.value)
    assert e.value.original.errno == errno.EIO and te.context["error"].endswith("primary failure")
    assert te.context["temp_cleanup"]["status"] == "failed" and "cannot remove" in te.context["temp_cleanup"]["error"]
    assert "could NOT be removed" in str(te) and "temporary file was removed" not in str(te)
    left = [p for p in f.parent.iterdir() if p.name.endswith(".tmp")]
    assert len(left) == 1 and te.context["temp_cleanup"]["path"] == str(left[0])


@pytest.mark.fast
@pytest.mark.parametrize("stage", ["temp_created", "attributes_set", "before_replace"])
@pytest.mark.parametrize("exc", [OSError(errno.ENOSPC, "simulated"), KeyboardInterrupt()])
def test_only_own_temp_removed(tmp_path, checkpoint, stage, exc):
    f = existing(tmp_path, 0o644)
    leftover = f.parent / ".x.json.999.deadbeef0000.tmp"; leftover.write_text("killed run")

    def boom(tmp):
        raise exc
    checkpoint[stage] = boom
    with pytest.raises(g._CacheWriteError if isinstance(exc, OSError) else KeyboardInterrupt) as e:
        write(f)
    if isinstance(exc, OSError):
        assert e.value.cleanup["status"] == "removed"
    else:
        assert any("cache temp cleanup" in n for n in getattr(e.value, "__notes__", []))
    assert f.read_text() == "old" and leftover.read_text() == "killed run"
    assert sorted(p.name for p in f.parent.iterdir()) == sorted([leftover.name, "x.json"])


# ===================================================================== permissions during the write (kept from v3)
def _window(f, monkeypatch, checkpoint):
    want = attrs(f); seen = []

    def fs(fd):
        if not _is_dir_fd(fd):
            st = os.fstat(fd); p = os.readlink(f"/proc/self/fd/{fd}")
            seen.append(("fsync", st.st_size, (stat.S_IMODE(st.st_mode), st.st_gid, xall(p))))
        return REAL_FSYNC(fd)
    monkeypatch.setattr(os, "fsync", fs)
    checkpoint["attributes_set"] = lambda tmp: seen.append(("attributes_set", tmp.stat().st_size, attrs(tmp)))
    checkpoint["before_replace"] = lambda tmp: seen.append(("before_replace", tmp.stat().st_size, attrs(tmp)))
    write(f, '{"new": true}\n' * 100)
    return want, seen


@pytest.mark.fast
@pytest.mark.parametrize("um", [0o022, 0o002, 0o000])
def test_no_looser_window_mode(tmp_path, monkeypatch, checkpoint, um):
    f = existing(tmp_path, 0o640); os.setxattr(f, "user.k", b"v")
    old = os.umask(um)
    try:
        want, seen = _window(f, monkeypatch, checkpoint)
    finally:
        os.umask(old)
    assert [s[0] for s in seen] == ["attributes_set", "fsync", "before_replace"] and seen[0][1] == 0
    assert all(s[2] == want for s in seen) and attrs(f) == want


@pytest.mark.fast
@needs_acl
def test_no_looser_window_extended_acl(tmp_path, monkeypatch, checkpoint):
    f = existing(tmp_path, 0o600); subprocess.run(["setfacl", "-m", f"u:{os.getuid()}:rw,g::r,o::-", str(f)], check=True)
    want, seen = _window(f, monkeypatch, checkpoint)
    assert seen[0][1] == 0 and all(s[2] == want for s in seen) and attrs(f) == want


@pytest.mark.fast
@pytest.mark.skipif(not SUPP, reason="no supplementary group")
def test_no_looser_window_group(tmp_path, monkeypatch, checkpoint):
    f = existing(tmp_path, 0o640); os.chown(f, -1, SUPP[0])
    want, seen = _window(f, monkeypatch, checkpoint)
    assert seen[0][1] == 0 and all(s[2] == want for s in seen) and attrs(f) == want


@pytest.mark.fast
@needs_acl
def test_new_file_like_open_under_default_acl(tmp_path):
    d = tmp_path / "c"; d.mkdir()
    subprocess.run(["setfacl", "-d", "-m", f"g:{os.getgid()}:rw,o::-", str(d)], check=True)
    with open(d / "ref.json", "w") as fh:
        fh.write("x")
    write(d / "x.json", "x")
    assert attrs(d / "x.json") == attrs(d / "ref.json")


@pytest.mark.fast
@needs_acl
def test_existing_without_acl_in_default_acl_directory(tmp_path):
    f = existing(tmp_path, 0o640)
    subprocess.run(["setfacl", "-d", "-m", f"g:{os.getgid()}:rw,o::-", str(f.parent)], check=True)
    want = attrs(f); write(f)
    assert attrs(f) == want


@pytest.mark.fast
def test_new_file_umask_like_open(tmp_path):
    d = tmp_path / "c"; d.mkdir(); old = os.umask(0o027)
    try:
        write(d / "x.json", "x")
        with open(d / "ref.json", "w") as fh:
            fh.write("x")
    finally:
        os.umask(old)
    assert attrs(d / "x.json") == attrs(d / "ref.json")


@pytest.mark.fast
def test_other_owner_refused_simulated(tmp_path, monkeypatch):
    """Simulated with a replaced os.getuid (no second uid available); real foreign ownership is untested."""
    f = existing(tmp_path, 0o664)
    monkeypatch.setattr(os, "getuid", lambda: os.stat(f).st_uid + 4242)
    with pytest.raises(g._CacheWriteError) as e:
        write(f)
    assert e.value.step == "inspect_existing" and e.value.original.errno == errno.EPERM
    assert f.read_text() == "old" and sorted(p.name for p in f.parent.iterdir()) == ["x.json"]


@pytest.mark.fast
def test_hard_link_alias_not_updated_recorded(tmp_path):
    f = existing(tmp_path, 0o644); os.link(f, tmp_path / "alias.json")
    write(f)
    assert f.read_text() == "new" and (tmp_path / "alias.json").read_text() == "old"


@pytest.mark.fast
def test_symlinked_cache_file_written_through(tmp_path):
    t = tmp_path / "t"; t.mkdir(); (t / "real.json").write_text("old"); d = tmp_path / "d"; d.mkdir()
    os.symlink(t / "real.json", d / "x.json")
    write(d / "x.json")
    assert (d / "x.json").is_symlink() and (t / "real.json").read_text() == "new"
    assert [p.name for p in t.iterdir()] == ["real.json"]


@pytest.mark.fast
def test_no_mutable_hook_globals():
    assert not hasattr(g, "_CACHE_WRITE_HOOK") and not hasattr(g, "_CACHE_WRITE_STAGE_HOOK")
    assert g._cache_write_checkpoint("before_replace", -1, "x") is None


# ===================================================================== round 57 (v5)
CONTENT = '{"same": "bytes"}\n' * 200


def _fds_to(path):
    r = []
    for fd in os.listdir("/proc/self/fd"):
        try:
            if os.readlink(f"/proc/self/fd/{fd}") == str(path):
                r.append(int(fd))
        except OSError:
            pass
    return r


def _run_threads(names, target):
    res = {}

    def run(n):
        try:
            res[n] = ("ok", target(n))
        except g._CacheWriteError as e:
            res[n] = ("error", (e.phase, e.step, e.replaced))
        except BaseException as e:          # noqa: BLE001 - recorded for the assertion
            res[n] = ("exc", repr(e))
    ts = [threading.Thread(target=run, args=(n,), name=n) for n in names]
    for t in ts:
        t.start()
    for t in ts:
        t.join(60)
    assert all(not t.is_alive() for t in ts)
    return res


@pytest.mark.fast
def test_two_writers_rename_between_rename_and_verify_same_bytes(tmp_path, monkeypatch):
    """Contract: cooperative writers of identical bytes all succeed; the name ends with one complete file."""
    f = existing(tmp_path, 0o644); os.setxattr(f, "user.k", b"v")
    a_renamed, b_done = threading.Event(), threading.Event()

    def cp(stage, dfd, tmp):
        if threading.current_thread().name == "A" and stage == "after_replace":
            a_renamed.set(); assert b_done.wait(30)
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)

    def target(n):
        if n == "B":
            assert a_renamed.wait(30)
        try:
            return write(f, CONTENT)
        finally:
            if n == "B":
                b_done.set()
    res = _run_threads(["A", "B"], target)
    assert res == {"A": ("ok", None), "B": ("ok", None)}
    assert f.read_text() == CONTENT and xall(f) == {"user.k": b"v"} and stat.S_IMODE(f.stat().st_mode) == 0o644
    assert sorted(p.name for p in f.parent.iterdir()) == ["x.json"]


@pytest.mark.fast
@pytest.mark.parametrize("hold", ["after_replace", "before_replace"])
def test_four_writers_synchronised(tmp_path, monkeypatch, hold):
    """4 writers stop together at one stage (all renamed / all ready to rename), then continue: all succeed."""
    f = existing(tmp_path, 0o640)
    barrier = threading.Barrier(4, timeout=30)

    def cp(stage, dfd, tmp):
        if stage == hold:
            barrier.wait()
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)
    res = _run_threads([f"W{i}" for i in range(4)], lambda n: write(f, CONTENT))
    assert all(v == ("ok", None) for v in res.values()), res
    assert f.read_text() == CONTENT and stat.S_IMODE(f.stat().st_mode) == 0o640
    assert sorted(p.name for p in f.parent.iterdir()) == ["x.json"] and _fds_to(f.parent) == []


@pytest.mark.fast
def test_four_writers_one_held_after_rename_others_finish(tmp_path, monkeypatch):
    f = existing(tmp_path, 0o644)
    w0_renamed, others_done = threading.Event(), threading.Event()
    done = []

    def cp(stage, dfd, tmp):
        if threading.current_thread().name == "W0" and stage == "after_replace":
            w0_renamed.set(); assert others_done.wait(30)
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)

    def target(n):
        if n != "W0":
            assert w0_renamed.wait(30)
        try:
            return write(f, CONTENT)
        finally:
            if n != "W0":
                done.append(n)
                if len(done) == 3:
                    others_done.set()
    res = _run_threads([f"W{i}" for i in range(4)], target)
    assert all(v == ("ok", None) for v in res.values()), res
    assert f.read_text() == CONTENT


@pytest.mark.fast
def test_concurrent_writer_with_different_bytes_is_not_success(tmp_path, monkeypatch):
    f = existing(tmp_path, 0o644)
    a_renamed, b_done = threading.Event(), threading.Event()

    def cp(stage, dfd, tmp):
        if threading.current_thread().name == "A" and stage == "after_replace":
            a_renamed.set(); assert b_done.wait(30)
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)

    def target(n):
        if n == "B":
            assert a_renamed.wait(30)
        try:
            return write(f, CONTENT if n == "A" else "DIFFERENT\n")
        finally:
            if n == "B":
                b_done.set()
    res = _run_threads(["A", "B"], target)
    assert res["B"] == ("ok", None) and res["A"] == ("error", ("after_replace", "verify_target", True))
    assert f.read_text() == "DIFFERENT\n"                       # one complete file remains


@pytest.mark.fast
def test_identical_regular_file_placed_after_rename_is_accepted(tmp_path, monkeypatch):
    f = existing(tmp_path)

    def rn(src, dst, *a, **k):
        REAL_RENAME(src, dst, *a, **k)
        d = k.get("dst_dir_fd"); REAL_UNLINK(dst, dir_fd=d)
        fd = os.open(dst, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644, dir_fd=d)
        os.write(fd, CONTENT.encode()); REAL_CLOSE(fd)
    monkeypatch.setattr(os, "rename", rn)
    assert write(f, CONTENT) is None and f.read_text() == CONTENT


@pytest.mark.fast
def test_first_fstat_eio_temp_identified_later_and_removed(tmp_path, monkeypatch):
    f = existing(tmp_path); real = os.fstat; hit = {"n": 0}

    def fs(fd):
        if os.readlink(f"/proc/self/fd/{fd}").endswith(".tmp") and hit["n"] == 0:
            hit["n"] += 1; raise OSError(errno.EIO, "injected")
        return real(fd)
    monkeypatch.setattr(os, "fstat", fs)
    with pytest.raises(g._CacheWriteError) as e:
        write(f)
    assert e.value.step == "create_temp" and e.value.cleanup["status"] == "removed" and e.value.cleanup["created"]
    assert sorted(p.name for p in f.parent.iterdir()) == ["x.json"] and f.read_text() == "old"


@pytest.mark.fast
def test_unidentifiable_temp_is_reported_and_left(tmp_path, monkeypatch):
    f = existing(tmp_path); real = os.fstat

    def fs(fd):
        if os.readlink(f"/proc/self/fd/{fd}").endswith(".tmp"):
            raise OSError(errno.EIO, "injected, persistent")
        return real(fd)
    monkeypatch.setattr(os, "fstat", fs)
    with pytest.raises(g._CacheWriteError) as e:
        write(f)
    c = e.value.cleanup
    assert e.value.step == "create_temp" and c["status"] == "left_unidentified" and c["created"] is True
    left = [p for p in f.parent.iterdir() if p.name.endswith(".tmp")]
    assert len(left) == 1 and c["path"] == str(left[0]) and f.read_text() == "old"
    msg = str(failed(e.value))
    assert "was NOT removed" in msg and "temporary file was removed" not in msg and str(left[0]) in msg
    assert _fds_to(f.parent) == [] and not any(os.readlink(f"/proc/self/fd/{x}").endswith(".tmp")
                                              for x in os.listdir("/proc/self/fd") if os.path.exists(f"/proc/self/fd/{x}"))


@pytest.mark.fast
def test_keyboard_interrupt_after_rename_closes_directory_fd(tmp_path, monkeypatch):
    f = existing(tmp_path)

    def cp(stage, dfd, tmp):
        if stage == "after_replace":
            raise KeyboardInterrupt
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)
    with pytest.raises(KeyboardInterrupt) as e:
        write(f)
    assert any("already holds the new file" in n for n in getattr(e.value, "__notes__", []))
    assert _fds_to(f.parent) == [] and f.read_text() == "new"


@pytest.mark.fast
def test_close_failure_does_not_hide_primary_keyboard_interrupt(tmp_path, monkeypatch):
    f = existing(tmp_path)

    def cp(stage, dfd, tmp):
        if stage == "after_replace":
            raise KeyboardInterrupt

    def cl(fd):
        is_dir = _is_dir_fd(fd)
        REAL_CLOSE(fd)
        if is_dir:
            raise OSError(errno.EIO, "close failed")
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)
    monkeypatch.setattr(os, "close", cl)
    with pytest.raises(KeyboardInterrupt):
        write(f)
    assert _fds_to(f.parent) == []


@pytest.mark.fast
@pytest.mark.parametrize("stage", ["temp_created", "attributes_set", "before_replace", "after_replace"])
@pytest.mark.parametrize("exc", [OSError(errno.EIO, "injected"), KeyboardInterrupt()])
def test_no_descriptor_leak_on_any_interruption(tmp_path, monkeypatch, stage, exc):
    f = existing(tmp_path)

    def cp(s, dfd, tmp):
        if s == stage:
            raise exc
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)
    before = set(os.listdir("/proc/self/fd"))
    with pytest.raises((g._CacheWriteError, KeyboardInterrupt)):
        write(f)
    after = set(os.listdir("/proc/self/fd"))
    assert _fds_to(f.parent) == [] and len(after - before) <= 0, after - before


# ===================================================================== round 58 (v6)
import signal
import time

CHILD = Path(__file__).with_name("_cache_v7_child.py")


@pytest.fixture
def no_hang():
    def handler(*a):
        raise AssertionError("blocked for more than 5 s")
    old = signal.signal(signal.SIGALRM, handler); signal.alarm(5)
    yield
    signal.alarm(0); signal.signal(signal.SIGALRM, old)


# --- 1: replacement by an identical cooperative writer between inspection steps
@pytest.mark.fast
def test_existing_replaced_identically_during_inspection(tmp_path, monkeypatch):
    f = existing(tmp_path, 0o640); os.setxattr(f, "user.k", b"v"); done = {"b": None}

    def op(p, flags, mode=0o777, *a, dir_fd=None, **k):
        if (done["b"] is None and str(p) == f.name and dir_fd is not None
                and not flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT)):
            monkeypatch.setattr(os, "open", REAL_OPEN)
            try:
                done["b"] = write(f, CONTENT)                       # B completes in between
            finally:
                monkeypatch.setattr(os, "open", op)
        return REAL_OPEN(p, flags, mode, *a, dir_fd=dir_fd, **k)
    monkeypatch.setattr(os, "open", op)
    assert write(f, CONTENT) is None                                # A
    monkeypatch.setattr(os, "open", REAL_OPEN)
    assert done["b"] is None and f.read_text() == CONTENT
    assert stat.S_IMODE(f.stat().st_mode) == 0o640 and xall(f) == {"user.k": b"v"}
    assert sorted(p.name for p in f.parent.iterdir()) == ["x.json"]


# --- 2: new cache published without clobbering; first publisher's attributes stay
def _two_children(tmp_path, first, second, texts=None, flags=()):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"; sync = tmp_path / "s"; sync.mkdir()
    texts = texts or {first: CONTENT, second: CONTENT}
    ps = {}
    for um in (first, second):
        tf = tmp_path / f"text_{um}"; tf.write_text(texts[um])
        ps[um] = subprocess.Popen([sys.executable, str(CHILD), str(BODY), str(f), um, str(sync), um, str(tf), *flags],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    for um in (first, second):
        t0 = time.monotonic()
        while not (sync / f"{um}.ready").exists():
            assert ps[um].poll() is None, ps[um].stderr.read()
            assert time.monotonic() - t0 < 60
            time.sleep(0.01)
    res = {}
    for um in (first, second):
        (sync / f"{um}.go").touch(); o, e = ps[um].communicate(timeout=60); res[um] = json.loads(o.strip().splitlines()[-1])
    return f, res


@pytest.mark.integration
@pytest.mark.parametrize("first,second,mode", [("077", "022", 0o600), ("022", "077", 0o644)])
def test_new_cache_first_publisher_attributes_stay(tmp_path, first, second, mode):
    f, res = _two_children(tmp_path, first, second)
    assert res == {first: {"result": "ok"}, second: {"result": "ok"}}
    assert f.read_text() == CONTENT and stat.S_IMODE(f.stat().st_mode) == mode
    assert sorted(p.name for p in f.parent.iterdir()) == ["x.json"]


@pytest.mark.integration
def test_new_cache_different_bytes_second_refused_first_kept(tmp_path):
    f, res = _two_children(tmp_path, "022", "077", {"022": "FIRST\n", "077": CONTENT})
    assert res["022"] == {"result": "ok"}
    assert res["077"]["result"] == "error" and res["077"]["phase"] == "before_replace"
    assert res["077"]["step"] == "verify_concurrent_existing" and res["077"]["cleanup"]["status"] == "removed"
    assert f.read_text() == "FIRST\n" and sorted(p.name for p in f.parent.iterdir()) == ["x.json"]


@pytest.mark.fast
def test_existing_strict_cache_not_broadened_by_lax_umask(tmp_path):
    f = existing(tmp_path, 0o600); old = os.umask(0o022)
    try:
        write(f, CONTENT)
    finally:
        os.umask(old)
    assert stat.S_IMODE(f.stat().st_mode) == 0o600


@pytest.mark.fast
def test_new_cache_threads_released_together(tmp_path, monkeypatch):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"; barrier = threading.Barrier(4, timeout=30)

    def cp(stage, dfd, tmp):
        if stage == "before_replace":
            barrier.wait()
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)
    res = _run_threads([f"W{i}" for i in range(4)], lambda n: write(f, CONTENT))
    assert all(v == ("ok", None) for v in res.values()), res
    assert f.read_text() == CONTENT and sorted(p.name for p in d.iterdir()) == ["x.json"] and _fds_to(d) == []


@pytest.mark.fast
def test_link_unsupported_publishes_with_rename_noreplace(tmp_path, monkeypatch):
    """v7: no silent clobbering rename; renameat2(RENAME_NOREPLACE) is used instead."""
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"

    def nolink(*a, **k):
        raise OSError(errno.EPERM, "simulated: hard links not supported")
    monkeypatch.setattr(os, "link", nolink)
    monkeypatch.setattr(os, "rename", lambda *a, **k: pytest.fail("clobbering rename used"))
    assert write(f, CONTENT) is None and f.read_text() == CONTENT
    assert sorted(p.name for p in d.iterdir()) == ["x.json"] and f.stat().st_nlink == 1


@pytest.mark.fast
def test_identical_existing_found_but_own_temp_not_removable(tmp_path, checkpoint, monkeypatch):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"

    def other_publishes(tmp):
        f.write_text(CONTENT)
    checkpoint["before_replace"] = other_publishes

    def ul(path, *a, **k):
        if str(path).endswith(".tmp"):
            raise PermissionError(errno.EACCES, "cannot remove")
        return REAL_UNLINK(path, *a, **k)
    monkeypatch.setattr(os, "unlink", ul)
    with pytest.raises(g._CacheWriteError) as e:
        write(f, CONTENT)
    assert (e.value.phase, e.value.step, e.value.replaced) == ("adopt_existing", "remove_temp", False)
    assert e.value.cleanup["status"] == "failed" and f.read_text() == CONTENT


# --- 3: a cache this user could not read back is not published
@pytest.mark.fast
def test_new_cache_unreadable_by_umask_not_published(tmp_path):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"; old = os.umask(0o466)
    try:
        with pytest.raises(g._CacheWriteError) as e:
            write(f, CONTENT)
    finally:
        os.umask(old)
    assert (e.value.phase, e.value.step, e.value.replaced) == ("before_replace", "check_readable", False)
    assert e.value.original.errno == errno.EACCES and e.value.cleanup["status"] == "removed"
    assert list(d.iterdir()) == []


@pytest.mark.fast
def test_existing_write_only_cache_refused_compat_change(tmp_path):
    """Documented compatibility change: HEAD rewrites a write-only cache in place; v6 refuses (cannot read it)."""
    f = existing(tmp_path, 0o200)
    try:
        with pytest.raises(g._CacheWriteError) as e:
            write(f)
        assert e.value.step == "inspect_existing" and e.value.original.errno == errno.EACCES
    finally:
        os.chmod(f, 0o600)
    assert f.read_text() == "old"


# --- 4: special files never block
@pytest.mark.fast
def test_fifo_after_rename_does_not_block(tmp_path, monkeypatch, no_hang):
    f = existing(tmp_path)

    def rn(src, dst, *a, **k):
        REAL_RENAME(src, dst, *a, **k); dd = k.get("dst_dir_fd")
        REAL_UNLINK(dst, dir_fd=dd); os.mkfifo(dst, 0o644, dir_fd=dd)
    monkeypatch.setattr(os, "rename", rn)
    t0 = time.monotonic()
    with pytest.raises(g._CacheWriteError) as e:
        write(f, CONTENT)
    assert time.monotonic() - t0 < 2
    assert (e.value.phase, e.value.step, e.value.replaced) == ("after_replace", "verify_target", True)
    assert stat.S_ISFIFO(os.lstat(f).st_mode) and _fds_to(f.parent) == []


@pytest.mark.fast
def test_fifo_after_link_publish_does_not_block(tmp_path, monkeypatch, no_hang):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"; real_link = os.link

    def ln(src, dst, *a, **k):
        real_link(src, dst, *a, **k); dd = k.get("dst_dir_fd")
        REAL_UNLINK(dst, dir_fd=dd); os.mkfifo(dst, 0o644, dir_fd=dd)
    monkeypatch.setattr(os, "link", ln)
    with pytest.raises(g._CacheWriteError) as e:
        write(f, CONTENT)
    assert (e.value.phase, e.value.step) == ("after_replace", "verify_target")
    assert sorted(p.name for p in d.iterdir()) == ["x.json"]           # our second name was removed


@pytest.mark.fast
def test_fifo_at_name_before_publish_new_does_not_block(tmp_path, checkpoint, no_hang):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"
    checkpoint["before_replace"] = lambda tmp: os.mkfifo(f, 0o644)
    with pytest.raises(g._CacheWriteError) as e:
        write(f, CONTENT)
    assert (e.value.phase, e.value.step) == ("before_replace", "verify_concurrent_existing")
    assert e.value.cleanup["status"] == "removed" and stat.S_ISFIFO(os.lstat(f).st_mode)


@pytest.mark.fast
def test_fifo_as_existing_cache_does_not_block(tmp_path, no_hang):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"; os.mkfifo(f, 0o644)
    with pytest.raises(g._CacheWriteError) as e:
        write(f, CONTENT)
    assert (e.value.phase, e.value.step) == ("before_replace", "inspect_existing")
    assert e.value.cleanup == {"status": "not_created", "created": False}



# ===================================================================== round 59 (v7)
def _dir_fsync_counter(monkeypatch):
    counts = {}

    def fs(fd):
        if _is_dir_fd(fd):
            n = threading.current_thread().name
            counts[n] = counts.get(n, 0) + 1
        return REAL_FSYNC(fd)
    monkeypatch.setattr(os, "fsync", fs)
    return counts


def _adoption_race(tmp_path, monkeypatch, b_extra=None):
    """A publishes a NEW cache by link and is held before its directory fsync; B (which inspected before A's
    publication) then adopts it.  b_extra(stage) may inject failures into B."""
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"
    a_pub, b_done, b_inspected = threading.Event(), threading.Event(), threading.Event()

    def cp(stage, dfd, tmp):
        n = threading.current_thread().name
        if n == "A" and stage == "temp_created":
            assert b_inspected.wait(30)
        if n == "A" and stage == "after_replace":
            a_pub.set(); assert b_done.wait(30)
        if n == "B" and stage == "before_replace":
            b_inspected.set(); assert a_pub.wait(30)
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)

    def target(n):
        try:
            return write(f, CONTENT)
        finally:
            if n == "B":
                b_done.set()

    def run(n):
        try:
            res[n] = ("ok", target(n))
        except g._CacheWriteError as e:
            res[n] = ("error", e)
        except BaseException as e:          # noqa: BLE001 - recorded for the assertion
            res[n] = ("exc", e)
    res = {}
    ts = [threading.Thread(target=run, args=(n,), name=n) for n in ("A", "B")]
    for t in ts:
        t.start()
    for t in ts:
        t.join(60)
    return f, res


@pytest.mark.fast
def test_adopting_writer_fsyncs_directory_itself(tmp_path, monkeypatch):
    counts = _dir_fsync_counter(monkeypatch)
    f, res = _adoption_race(tmp_path, monkeypatch)
    assert res["A"] == ("ok", None) and res["B"] == ("ok", None)
    assert counts.get("B") == 1 and counts.get("A") == 1
    assert f.read_text() == CONTENT and sorted(p.name for p in f.parent.iterdir()) == ["x.json"]
    assert _fds_to(f.parent) == []


@pytest.mark.fast
@pytest.mark.parametrize("eno", [errno.ENOTSUP, getattr(errno, "EOPNOTSUPP", errno.ENOTSUP)])
def test_adopting_writer_fsync_unsupported_warns_not_replaced(tmp_path, monkeypatch, eno):
    def fs(fd):
        if _is_dir_fd(fd) and threading.current_thread().name == "B":
            raise OSError(eno, os.strerror(eno))
        return REAL_FSYNC(fd)
    monkeypatch.setattr(os, "fsync", fs)
    f, res = _adoption_race(tmp_path, monkeypatch)
    assert res["A"] == ("ok", None)
    kind, warn = res["B"]
    assert kind == "ok" and warn["replaced"] is False and warn["adopted_existing"] is True
    assert warn["durability"] == "not_guaranteed" and warn["errno"] == eno


@pytest.mark.fast
@pytest.mark.parametrize("eno", [errno.EIO, errno.EINVAL])
def test_adopting_writer_fsync_failure_exit3_not_replaced(tmp_path, monkeypatch, eno):
    def fs(fd):
        if _is_dir_fd(fd) and threading.current_thread().name == "B":
            raise OSError(eno, os.strerror(eno))
        return REAL_FSYNC(fd)
    monkeypatch.setattr(os, "fsync", fs)
    f, res = _adoption_race(tmp_path, monkeypatch)
    kind, e = res["B"]
    assert kind == "error" and (e.phase, e.step, e.replaced) == ("adopt_existing", "fsync_directory", False)
    assert e.cleanup["status"] == "removed" and "did not change the cache name" in str(failed(e))
    assert f.read_text() == CONTENT and _fds_to(f.parent) == []


@pytest.mark.fast
def test_adopting_writer_close_failure_not_suppressed(tmp_path, monkeypatch):
    def cl(fd):
        is_dir = _is_dir_fd(fd)
        REAL_CLOSE(fd)
        if is_dir and threading.current_thread().name == "B":
            raise OSError(errno.EIO, "injected close failure")
    monkeypatch.setattr(os, "close", cl)
    f, res = _adoption_race(tmp_path, monkeypatch)
    kind, e = res["B"]
    assert kind == "error" and (e.phase, e.step, e.replaced) == ("adopt_existing", "close_directory", False)
    assert res["A"] == ("ok", None) and f.read_text() == CONTENT


@pytest.mark.fast
def test_adopting_writer_keyboard_interrupt_closes_fd(tmp_path, monkeypatch):
    real = REAL_FSYNC

    def fs(fd):
        if _is_dir_fd(fd) and threading.current_thread().name == "B":
            raise KeyboardInterrupt
        return real(fd)
    monkeypatch.setattr(os, "fsync", fs)
    f, res = _adoption_race(tmp_path, monkeypatch)
    kind, e = res["B"]
    assert kind == "exc" and isinstance(e, KeyboardInterrupt)
    assert any("was not changed by this run" in n for n in e.__notes__)
    assert res["A"] == ("ok", None) and f.read_text() == CONTENT and _fds_to(f.parent) == []


# --- link unsupported: renameat2(RENAME_NOREPLACE) or explicit failure; never clobber
@pytest.mark.fast
def test_link_unsupported_other_content_appears_not_clobbered(tmp_path, checkpoint, monkeypatch):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"

    def nolink(*a, **k):
        raise OSError(errno.EPERM, "simulated")
    monkeypatch.setattr(os, "link", nolink)
    checkpoint["before_replace"] = lambda tmp: f.write_text("OTHER WRITER\n")
    with pytest.raises(g._CacheWriteError) as e:
        write(f, CONTENT)
    assert (e.value.phase, e.value.step) == ("before_replace", "verify_concurrent_existing")
    assert f.read_text() == "OTHER WRITER\n" and e.value.cleanup["status"] == "removed"


@pytest.mark.fast
@pytest.mark.parametrize("eno", [errno.ENOSYS, errno.EINVAL, errno.ENOTSUP])
def test_no_link_and_no_noreplace_fails_explicitly(tmp_path, monkeypatch, eno):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"

    def nolink(*a, **k):
        raise OSError(errno.EPERM, "simulated")

    def nonr(*a, **k):
        raise OSError(eno, "simulated renameat2 failure")
    monkeypatch.setattr(os, "link", nolink)
    monkeypatch.setattr(g, "_cache_rename_noreplace", nonr)
    monkeypatch.setattr(os, "rename", lambda *a, **k: pytest.fail("clobbering rename used"))
    with pytest.raises(g._CacheWriteError) as e:
        write(f, CONTENT)
    assert (e.value.phase, e.value.step, e.value.replaced) == ("before_replace", "publish_new_noreplace", False)
    assert "no non-clobbering way" in str(e.value.original) and "nothing was published" in str(failed(e.value))
    assert list(d.iterdir()) == []


@pytest.mark.integration
@pytest.mark.parametrize("first,second,mode", [("077", "022", 0o600), ("022", "077", 0o644)])
def test_no_link_children_first_publisher_attributes_stay(tmp_path, first, second, mode):
    f, res = _two_children(tmp_path, first, second, flags=("nolink",))
    assert res == {first: {"result": "ok"}, second: {"result": "ok"}}
    assert f.read_text() == CONTENT and stat.S_IMODE(f.stat().st_mode) == mode
    assert sorted(p.name for p in f.parent.iterdir()) == ["x.json"]


@pytest.mark.integration
def test_no_link_children_different_content_second_refused(tmp_path):
    f, res = _two_children(tmp_path, "022", "077", {"022": "FIRST\n", "077": CONTENT}, flags=("nolink",))
    assert res["022"] == {"result": "ok"}
    assert (res["077"]["result"], res["077"]["step"]) == ("error", "verify_concurrent_existing")
    assert f.read_text() == "FIRST\n" and sorted(p.name for p in f.parent.iterdir()) == ["x.json"]


@pytest.mark.integration
def test_no_link_no_noreplace_children_nothing_published(tmp_path):
    f, res = _two_children(tmp_path, "022", "077", flags=("nolink", "noreplace_unavailable"))
    assert all(v["result"] == "error" and v["step"] == "publish_new_noreplace" for v in res.values()), res
    assert not f.exists() and list(f.parent.iterdir()) == []


# --- link publication: the temporary second name is cleaned up and reported truthfully
@pytest.mark.fast
@pytest.mark.parametrize("exc", [OSError(errno.EIO, "injected"), KeyboardInterrupt()])
def test_error_right_after_link_removes_temp_name(tmp_path, monkeypatch, exc):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"

    def cp(stage, dfd, tmp):
        if stage == "after_replace":
            raise exc
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)
    if isinstance(exc, OSError):
        with pytest.raises(g._CacheWriteError) as e:
            write(f, CONTENT)
        assert (e.value.phase, e.value.step, e.value.replaced) == ("after_replace", "after_rename", True)
        assert e.value.cleanup["status"] == "removed"
        assert "temporary file was removed" in str(failed(e.value))
    else:
        with pytest.raises(KeyboardInterrupt) as e:
            write(f, CONTENT)
        assert any("'status': 'removed'" in n for n in e.value.__notes__)
    assert f.read_text() == CONTENT and f.stat().st_nlink == 1
    assert sorted(p.name for p in d.iterdir()) == ["x.json"] and _fds_to(d) == []


@pytest.mark.fast
def test_error_after_link_temp_removal_failure_reported(tmp_path, monkeypatch):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"

    def cp(stage, dfd, tmp):
        if stage == "after_replace":
            raise OSError(errno.EIO, "injected")

    def ul(path, *a, **k):
        if str(path).endswith(".tmp"):
            raise PermissionError(errno.EACCES, "cannot remove")
        return REAL_UNLINK(path, *a, **k)
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)
    monkeypatch.setattr(os, "unlink", ul)
    with pytest.raises(g._CacheWriteError) as e:
        write(f, CONTENT)
    c = e.value.cleanup
    assert c["status"] == "failed" and "cannot remove" in c["error"]
    left = [p for p in d.iterdir() if p.name.endswith(".tmp")]
    assert len(left) == 1 and c["path"] == str(left[0]) and os.stat(left[0]).st_ino == f.stat().st_ino
    assert f.read_text() == CONTENT and "could NOT be removed" in str(failed(e.value))


@pytest.mark.fast
def test_temp_name_removal_failure_in_normal_flow_is_error(tmp_path, monkeypatch):
    d = tmp_path / "c"; d.mkdir(); f = d / "x.json"

    def ul(path, *a, **k):
        if str(path).endswith(".tmp"):
            raise PermissionError(errno.EACCES, "cannot remove")
        return REAL_UNLINK(path, *a, **k)
    monkeypatch.setattr(os, "unlink", ul)
    with pytest.raises(g._CacheWriteError) as e:
        write(f, CONTENT)
    assert (e.value.phase, e.value.step) == ("after_replace", "remove_temp_name")
    assert e.value.cleanup["status"] == "failed" and f.read_text() == CONTENT


@pytest.mark.fast
def test_rename_publication_reports_not_applicable(tmp_path, monkeypatch):
    f = existing(tmp_path)

    def cp(stage, dfd, tmp):
        if stage == "after_replace":
            raise OSError(errno.EIO, "injected")
    monkeypatch.setattr(g, "_cache_write_checkpoint", cp)
    with pytest.raises(g._CacheWriteError) as e:
        write(f, CONTENT)
    assert e.value.cleanup["status"] == "not_applicable" and f.read_text() == CONTENT
