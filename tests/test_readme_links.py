"""README.md の相対リンクの参照先が tree にあることを確かめる。"""

import re
from pathlib import Path

_LINK_OPEN = re.compile(r"!?\[[^\]]*\]\(")
_SCHEME = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*:")
_ESCAPE = re.compile(r"\\([!-/:-@\[-`{-~])")


def _unclosed(text):
    """`](<` で始まり、同じ行に閉じの `>` が無いものの行番号を返す。"""
    lines = []
    for m in _LINK_OPEN.finditer(text):
        i = m.end()
        while i < len(text) and text[i] in " \t":
            i += 1
        if not text.startswith("<", i):
            continue
        end = i + 1
        while end < len(text) and text[end] not in ">\n":
            end += 2 if text[end] == "\\" else 1
        if not (end < len(text) and text[end] == ">"):
            lines.append(text.count("\n", 0, i) + 1)
    return lines


def _targets(text):
    """`](` の後の target を取り出す（タイトルは含めない）。"""
    for m in _LINK_OPEN.finditer(text):
        i = m.end()
        while i < len(text) and text[i] in " \t":
            i += 1
        if text.startswith("<", i):
            end = i + 1
            while end < len(text) and text[end] not in ">\n":
                end += 2 if text[end] == "\\" else 1
            if end < len(text) and text[end] == ">":
                yield text[i + 1 : end]
            continue
        depth = 0
        end = i
        while end < len(text) and text[end] != "\n":
            ch = text[end]
            if ch == "\\":
                end += 2
                continue
            if ch == "(":
                depth += 1
            elif ch == ")":
                if depth == 0:
                    break
                depth -= 1
            elif ch in " \t" and depth == 0:
                break
            end += 1
        yield text[i : min(end, len(text))]


def test_readme_relative_links_exist():
    root = Path(__file__).resolve().parent.parent
    readme = root / "README.md"
    text = readme.read_text(encoding="utf-8")

    checked = []
    missing = []
    for raw in _targets(text):
        target = raw.strip()
        if _SCHEME.match(target) or target.startswith("#"):
            continue
        if target.startswith("<") and target.endswith(">"):
            target = target[1:-1]
        target = target.split("#", 1)[0].split("?", 1)[0]
        target = _ESCAPE.sub(r"\1", target)
        if not target:
            continue
        checked.append(target)
        if not (readme.parent / target).exists():
            missing.append(target)

    unclosed = [f"README.md:{n}: unclosed <" for n in _unclosed(text)]

    assert not unclosed, f"unclosed angle brackets: {unclosed}"
    assert len(checked) >= 5, f"relative links checked: {len(checked)}"
    assert not missing, f"missing link targets: {missing}"
