"""README.md の相対リンクの参照先が tree にあることを確かめる。"""

import re
from pathlib import Path

from test_markdown_table_columns import DOCS, FENCE_OPEN, FENCE_RUN

_LINK_OPEN = re.compile(r"!?\[[^\]]*\]\(")
_SCHEME = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*:")
_ESCAPE = re.compile(r"\\([!-/:-@\[-`{-~])")


def _scan_angle(text, start):
    """`<` の位置 start から、閉じの `>` または改行（無ければ末尾）の位置を返す。"""
    end = start + 1
    while end < len(text) and text[end] not in ">\n":
        end += 2 if text[end] == "\\" and text[end + 1 : end + 2] != "\n" else 1
    return end


def _unclosed(text):
    """`](<` で始まり、同じ行に閉じの `>` が無いものの行番号を返す。"""
    lines = []
    for m in _LINK_OPEN.finditer(text):
        i = m.end()
        while i < len(text) and text[i] in " \t":
            i += 1
        if not text.startswith("<", i):
            continue
        end = _scan_angle(text, i)
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
            end = _scan_angle(text, i)
            if end < len(text) and text[end] == ">":
                yield text[i + 1 : end]
            continue
        depth = 0
        end = i
        while end < len(text) and text[end] != "\n":
            ch = text[end]
            if ch == "\\" and text[end + 1 : end + 2] != "\n":
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


def test_backslash_before_newline_does_not_join_lines():
    text = "[a](x\\\n[b](y.md)\n"
    targets = list(_targets(text))
    assert all("\n" not in t for t in targets), targets
    assert "y.md" in targets, targets

    angle = "[a](<x\\\n[b](<y.md>)\n"
    assert _unclosed(angle) == [1]
    assert list(_targets(angle)) == ["y.md"]

    escaped = list(_targets("[a](x\\(1\\).md)"))
    assert [_ESCAPE.sub(r"\1", t) for t in escaped] == ["x(1).md"]

    assert list(_targets("[a](x\\")) == ["x\\"]


def test_link_destination_does_not_span_lines():
    for text in ("[a](x\ny.md)", "[a](x\\\ny.md)", "[a](<x\ny.md>)", "[a](<x\\\ny.md>)"):
        targets = list(_targets(text))
        assert all("\n" not in t for t in targets), (text, targets)

    assert list(_targets("[a](x.md)\n[a](<x y.md>)\n")) == ["x.md", "x y.md"]


def test_parentheses_in_destination_on_synthetic_input():
    def unescaped(text):
        return [_ESCAPE.sub(r"\1", t) for t in _targets(text)]

    assert list(_targets("[a](x(1).md)")) == ["x(1).md"]
    assert unescaped("[a](x\\).md)") == ["x).md"]
    assert unescaped("[a](x\\(.md)") == ["x(.md"]


_ANCHOR_LINK = re.compile(r"\]\(#([^)\s]*)\)")
_HEADING = re.compile(r"^ {0,3}#{1,6}[ \t]+(.*?)[ \t]*$")
_CODE_SPAN = re.compile(r"(`+)(.+?)\1")
_INLINE_LINK = re.compile(r"\[([^\]]*)\]\([^)]*\)")


def _outside_fences(text):
    """コードフェンスの外の (行番号, 行) を返す。"""
    fence = None
    for n, line in enumerate(text.split("\n"), 1):
        stripped = line.strip()
        if fence is None:
            m = FENCE_OPEN.match(stripped) if len(line) - len(line.lstrip(" ")) < 4 else None
            if m:
                fence = m.group(1)
                continue
            yield n, line
        else:
            m = FENCE_RUN.match(stripped)
            if m and m.group(1)[0] == fence[0] and len(m.group(1)) >= len(fence) and stripped == m.group(1):
                fence = None


def _slug(title):
    title = _INLINE_LINK.sub(r"\1", title)
    title = re.sub(r"\s+#+$", "", title).lower()
    title = re.sub(r"[^\w\- ]", "", title)
    return title.replace(" ", "-")


def _anchors(text):
    """見出しから GitHub の作り方でアンカーを作る。重なったら -1, -2 を付ける。"""
    seen = {}
    anchors = []
    for _, line in _outside_fences(text):
        m = _HEADING.match(line)
        if not m:
            continue
        base = _slug(m.group(1))
        count = seen.get(base, 0)
        seen[base] = count + 1
        anchors.append(base if count == 0 else f"{base}-{count}")
    return anchors


def _anchor_links(text):
    """コードフェンスとインラインコードの外の ](#...) を (行番号, 行き先) で返す。"""
    for n, line in _outside_fences(text):
        line = _CODE_SPAN.sub("", line)
        for m in _ANCHOR_LINK.finditer(line):
            yield n, m.group(1)


def test_in_document_anchor_links_match_headings():
    root = Path(__file__).resolve().parent.parent
    total = 0
    bad = []
    for rel in DOCS:
        text = (root / rel).read_text(encoding="utf-8")
        anchors = set(_anchors(text))
        for n, dest in _anchor_links(text):
            total += 1
            if dest not in anchors:
                bad.append(f"{rel}:{n}: #{dest}")
    assert total >= 1, f"anchor links checked: {total}"
    assert not bad, f"anchor links without a heading: {bad}"


def test_heading_anchor_rule_on_synthetic_input():
    assert _anchors("# Known Limitations\n") == ["known-limitations"]
    assert _anchors("## Chamber property cache\n") == ["chamber-property-cache"]
    assert _anchors("## Cache, writes (and) runs!\n") == ["cache-writes-and-runs"]
    assert _anchors("## 検証 status\n") == ["検証-status"]
    assert _anchors("## Dup\n\n## Dup\n\n## Dup\n") == ["dup", "dup-1", "dup-2"]
    assert _anchors("```\n# Not heading\n```\n## Real\n") == ["real"]
    assert list(_anchor_links("[a](#x) `[b](#y)` [c](#z)\n```\n[d](#w)\n```\n")) == [(1, "x"), (1, "z")]
