"""r434 PROTOTYPE tests (isolated copies only): cache v7 - what the README may and may not promise.
Two synchronised writer PROCESSES (own umask), held at chosen stages:
 * no-clobber path (both observed "no cache"): the first publisher's file stays (inode kept), the later one adopts;
 * a writer that observes the cache as EXISTING rebuilds it: same bytes, but the inode changes (attributes copied);
 * writers with DIFFERENT bytes are outside the premise: a later rebuild replaces the earlier cache (no refusal).
"""
import json
import os
import stat
import subprocess
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BODY = sorted(ROOT.glob("GRIBS_v*.py"))[0]
CHILD = Path(__file__).with_name("_cache_v7_sync_child.py")
CONTENT = '{"same": "bytes"}\n' * 200
pytestmark = pytest.mark.integration


class Writers:
    def __init__(self, tmp):
        self.tmp = tmp; self.d = tmp / "c"; self.d.mkdir(); self.f = self.d / "x.json"
        self.sync = tmp / "s"; self.sync.mkdir(); self.p = {}

    def start(self, name, umask, text, hold="-", *flags):
        tf = self.tmp / f"text_{name}"; tf.write_text(text)
        self.p[name] = subprocess.Popen([sys.executable, str(CHILD), str(BODY), str(self.f), umask, str(self.sync), name,
                                         str(tf), hold, *flags], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    def wait_held(self, name):
        t0 = time.monotonic()
        while not (self.sync / f"{name}.held").exists():
            assert self.p[name].poll() is None, self.p[name].stderr.read()[-1500:]
            assert time.monotonic() - t0 < 60
            time.sleep(0.01)

    def release(self, name):
        (self.sync / f"{name}.go").touch()

    def result(self, name):
        o, e = self.p[name].communicate(timeout=120)
        assert o.strip(), e[-1500:]
        return json.loads(o.strip().splitlines()[-1])

    def inode(self):
        return os.stat(self.f).st_ino


def test_sequential_rebuild_same_bytes_changes_inode_keeps_attributes(tmp_path):
    w = Writers(tmp_path)
    w.start("A", "077", CONTENT); assert w.result("A")["result"] == "ok"
    ino_a = w.inode(); assert stat.S_IMODE(os.stat(w.f).st_mode) == 0o600
    w.start("B", "022", CONTENT); assert w.result("B")["result"] == "ok"     # B observes an EXISTING cache
    assert w.inode() != ino_a                                             # rebuilt: new inode
    assert w.f.read_text() == CONTENT and stat.S_IMODE(os.stat(w.f).st_mode) == 0o600
    assert sorted(p.name for p in w.d.iterdir()) == ["x.json"]


def test_synchronised_publish_then_rebuild_changes_inode(tmp_path):
    """A publishes a NEW cache (link) and is held before verify/fsync; B starts afterwards, observes the cache as
    existing, rebuilds it (rename).  Both succeed with identical bytes; the name ends on B's inode."""
    w = Writers(tmp_path)
    w.start("A", "077", CONTENT, "after_replace"); w.wait_held("A")
    ino_a = w.inode()
    w.start("B", "022", CONTENT); rb = w.result("B")
    ino_b = w.inode()
    w.release("A"); ra = w.result("A")
    assert ra == {"result": "ok", "value": None} and rb == {"result": "ok", "value": None}
    assert ino_b != ino_a and w.inode() == ino_b
    assert w.f.read_text() == CONTENT and stat.S_IMODE(os.stat(w.f).st_mode) == 0o600   # copied from A's file
    assert sorted(p.name for p in w.d.iterdir()) == ["x.json"]


def test_synchronised_both_observed_absent_first_file_stays(tmp_path):
    """Contrast: both observed 'no cache' (held before publishing) -> no-clobber path: A's inode stays, B adopts."""
    w = Writers(tmp_path)
    w.start("A", "077", CONTENT, "before_replace"); w.start("B", "022", CONTENT, "before_replace")
    w.wait_held("A"); w.wait_held("B")
    w.release("A"); assert w.result("A")["result"] == "ok"
    ino_a = w.inode()
    w.release("B"); assert w.result("B")["result"] == "ok"
    assert w.inode() == ino_a and stat.S_IMODE(os.stat(w.f).st_mode) == 0o600
    assert sorted(p.name for p in w.d.iterdir()) == ["x.json"]


@pytest.mark.parametrize("flags", [(), ("nolink",)])
def test_both_absent_first_file_stays_also_without_hard_links(tmp_path, flags):
    w = Writers(tmp_path)
    w.start("A", "022", CONTENT, "before_replace", *flags); w.start("B", "077", CONTENT, "before_replace", *flags)
    w.wait_held("A"); w.wait_held("B")
    w.release("A"); assert w.result("A")["result"] == "ok"
    ino_a = w.inode()
    w.release("B"); assert w.result("B")["result"] == "ok"
    assert w.inode() == ino_a and stat.S_IMODE(os.stat(w.f).st_mode) == 0o644


def test_different_bytes_outside_premise_rebuild_replaces(tmp_path):
    """Outside the premise: a writer that observes the cache as existing replaces it even with different bytes;
    there is no refusal on this path."""
    w = Writers(tmp_path)
    w.start("A", "022", "FIRST\n"); assert w.result("A")["result"] == "ok"
    w.start("B", "022", "SECOND\n"); assert w.result("B")["result"] == "ok"
    assert w.f.read_text() == "SECOND\n"


def test_different_bytes_outside_premise_held_publisher_fails_after_replacement(tmp_path):
    w = Writers(tmp_path)
    w.start("A", "022", "FIRST\n", "after_replace"); w.wait_held("A")
    w.start("B", "022", "SECOND\n"); assert w.result("B")["result"] == "ok"
    w.release("A"); ra = w.result("A")
    assert (ra["result"], ra["phase"], ra["step"], ra["replaced"]) == ("error", "after_replace", "verify_target", True)
    assert w.f.read_text() == "SECOND\n"
