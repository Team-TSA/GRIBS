"""Test helper (not collected): one cache writer in its own process (own umask).
argv: <body.py> <cache file> <umask octal> <sync dir> <name> <text file> [nolink] [noreplace_unavailable]
  nolink: os.link raises EPERM (file system without hard links, simulated)
  noreplace_unavailable: renameat2(RENAME_NOREPLACE) reported unavailable (ENOSYS, simulated)
Waits at the 'before_replace' checkpoint until <sync>/<name>.go exists; prints a JSON result."""
import importlib.util, json, os, sys, time
from pathlib import Path
body, cache, um, sync, name, textf = sys.argv[1:7]; flags = set(sys.argv[7:])
spec = importlib.util.spec_from_file_location("gribs_child", body); m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m; spec.loader.exec_module(m)
os.umask(int(um, 8)); S = Path(sync)
def cp(stage, dfd, tmp):
    if stage == "before_replace":
        (S / f"{name}.ready").touch()
        while not (S / f"{name}.go").exists():
            time.sleep(0.01)
m._cache_write_checkpoint = cp
import errno
if "nolink" in flags:
    def _nolink(*a, **k):
        raise OSError(errno.EPERM, "simulated: hard links not supported")
    os.link = _nolink
if "noreplace_unavailable" in flags:
    def _nonr(*a, **k):
        raise OSError(errno.ENOSYS, "simulated: renameat2 not available")
    m._cache_rename_noreplace = _nonr
try:
    m._cache_atomic_write_text(Path(cache), Path(textf).read_text()); print(json.dumps({"result": "ok"}))
except m._CacheWriteError as e:
    print(json.dumps({"result": "error", "phase": e.phase, "step": e.step, "cleanup": e.cleanup}))
