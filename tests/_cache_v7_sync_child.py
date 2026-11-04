"""r434 PROTOTYPE test helper (not collected): one cache writer in its own process, held at a chosen stage.
argv: <body.py> <cache file> <umask octal> <sync dir> <name> <text file> <hold stage|-> [nolink]
At checkpoint <hold stage> (temp_created / attributes_set / before_replace / after_replace): touch
<sync>/<name>.held, then wait for <sync>/<name>.go.  Prints a JSON result line."""
import errno, importlib.util, json, os, sys, time
from pathlib import Path
body, cache, um, sync, name, textf, hold = sys.argv[1:8]; flags = set(sys.argv[8:])
spec = importlib.util.spec_from_file_location("gribs_sync_child", body); m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m; spec.loader.exec_module(m)
os.umask(int(um, 8)); S = Path(sync)
def cp(stage, dfd, tmp):
    if stage == hold:
        (S / f"{name}.held").touch()
        while not (S / f"{name}.go").exists():
            time.sleep(0.01)
m._cache_write_checkpoint = cp
if "nolink" in flags:
    def _nolink(*a, **k):
        raise OSError(errno.EPERM, "simulated: hard links not supported")
    os.link = _nolink
try:
    v = m._cache_atomic_write_text(Path(cache), Path(textf).read_text())
    print(json.dumps({"result": "ok", "value": v}))
except m._CacheWriteError as e:
    print(json.dumps({"result": "error", "phase": e.phase, "step": e.step, "replaced": e.replaced,
                      "cleanup": e.cleanup}))
