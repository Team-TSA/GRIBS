"""README.md の相対リンクの参照先が tree にあることを確かめる。"""

import re
from pathlib import Path

_LINK = re.compile(r"!?\[[^\]]*\]\(([^)]*)\)")
_SCHEME = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*:")


def test_readme_relative_links_exist():
    root = Path(__file__).resolve().parent.parent
    readme = root / "README.md"
    text = readme.read_text(encoding="utf-8")

    checked = []
    missing = []
    for raw in _LINK.findall(text):
        target = raw.strip()
        if _SCHEME.match(target) or target.startswith("#"):
            continue
        target = target.split(None, 1)[0] if target else ""
        if target.startswith("<") and target.endswith(">"):
            target = target[1:-1]
        target = target.split("#", 1)[0].split("?", 1)[0]
        if not target:
            continue
        checked.append(target)
        if not (readme.parent / target).exists():
            missing.append(target)

    assert len(checked) >= 5, f"relative links checked: {len(checked)}"
    assert not missing, f"missing link targets: {missing}"
