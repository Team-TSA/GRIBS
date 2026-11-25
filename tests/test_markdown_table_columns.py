"""User-facing Markdown documents: table header, delimiter and body rows have the same column count."""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

DOCS = [
    "README.md",
    "CONTRIBUTING.md",
    "Specification/README.md",
    "SECURITY.md",
    "CODE_OF_CONDUCT.md",
    ".github/pull_request_template.md",
    "CHANGELOG.md",
    "CHANGELOG_v0.5.0-alpha_A5.1-5_to_A7.md",
    "RELEASE_NOTES_v0.3.1-alpha.md",
    "RELEASE_CHECKLIST.md",
]

FENCE_OPEN = re.compile(r"^(`{3,}|~{3,})")
FENCE_RUN = re.compile(r"^(`+|~+)")
DELIM_CELL = re.compile(r"^\s*:?-+:?\s*$")


def _split_cells(line):
    """Split a table row on unescaped | outside inline code; outer pipes are dropped."""
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|") and not s.endswith("\\|"):
        s = s[:-1]
    cells = []
    cur = []
    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if c == "\\" and i + 1 < n:
            cur.append(s[i:i + 2])
            i += 2
            continue
        if c == "`":
            j = i
            while j < n and s[j] == "`":
                j += 1
            ticks = s[i:j]
            close = s.find(ticks, j)
            while close != -1 and close + len(ticks) < n and s[close + len(ticks)] == "`":
                # closing run must be exactly the same length
                k = close
                while k < n and s[k] == "`":
                    k += 1
                close = s.find(ticks, k)
            if close != -1:
                cur.append(s[i:close + len(ticks)])
                i = close + len(ticks)
                continue
            cur.append(ticks)
            i = j
            continue
        if c == "|":
            cells.append("".join(cur))
            cur = []
            i += 1
            continue
        cur.append(c)
        i += 1
    cells.append("".join(cur))
    return cells


def _is_delim(line):
    if "|" not in line:
        return False
    cells = _split_cells(line)
    return all(DELIM_CELL.match(c) for c in cells)


def _tables(text):
    """Return (header line no, header cols, [(line no, cols), ...]) per table; line numbers are 1-based.

    The first body entry is the delimiter row; the rest are the body rows.
    """
    lines = text.splitlines()
    out = []
    fence = None
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.lstrip()
        if fence is None:
            opener = FENCE_OPEN.match(stripped)
            if opener:
                fence = (opener.group(1)[0], len(opener.group(1)))
                i += 1
                continue
        else:
            closer = FENCE_RUN.match(stripped)
            if (
                closer
                and closer.group(1)[0] == fence[0]
                and len(closer.group(1)) >= fence[1]
                and stripped[closer.end():].strip() == ""
            ):
                fence = None
            i += 1
            continue
        if "|" in line and i + 1 < len(lines) and _is_delim(lines[i + 1]):
            rows = [(i + 2, len(_split_cells(lines[i + 1])))]
            j = i + 2
            while j < len(lines) and lines[j].strip() != "" and "|" in lines[j]:
                rows.append((j + 1, len(_split_cells(lines[j]))))
                j += 1
            out.append((i + 1, len(_split_cells(line)), rows))
            i = j
            continue
        i += 1
    return out


def _mismatches(text):
    bad = []
    for head_no, head_cols, rows in _tables(text):
        for no, cols in rows:
            if cols != head_cols:
                bad.append((no, cols, head_cols))
    return bad


@pytest.mark.parametrize("rel", DOCS)
def test_table_columns_match_header(rel):
    text = (ROOT / rel).read_text(encoding="utf-8")
    bad = [
        f"{rel}:{no}: {cols} columns, header has {head}"
        for no, cols, head in _mismatches(text)
    ]
    assert bad == [], "table row column count differs from header:\n" + "\n".join(bad)


def test_table_rule_on_synthetic_input():
    # (i) matching table: no mismatches
    ok = "| a | b |\n|---|:-:|\n| 1 | 2 |\n| 3 | 4 |\n"
    assert _mismatches(ok) == []
    assert _tables(ok) == [(1, 2, [(2, 2), (3, 2), (4, 2)])]
    # (ii) a body row with an extra column is reported by line number
    bad = "text\n\n| a | b |\n|---|---|\n| 1 | 2 |\n| 3 | 4 | 5 |\n"
    assert _mismatches(bad) == [(6, 3, 2)]
    # (iii) escaped pipe and pipe inside inline code do not add columns
    esc = "| a | b |\n|---|---|\n| x \\| y | z |\n"
    assert _mismatches(esc) == []
    code = "| a | b |\n|---|---|\n| `x | y` | z |\n"
    assert _mismatches(code) == []
    # (iv) tables inside code fences are not counted
    fenced = "```\n| a | b |\n|---|---|\n| 1 | 2 | 3 |\n```\n"
    assert _tables(fenced) == []
    tilde = "~~~~\n| a |\n|---|\n| 1 | 2 |\n~~~\n| a |\n|---|\n| 1 | 2 |\n~~~~\n"
    assert _tables(tilde) == []


def test_table_split_trailing_escaped_pipe_and_backtick_run_length():
    # (a) a trailing \| is cell text, not an outer pipe
    line_a = "| a | b \\|"
    cells_a = _split_cells(line_a)
    assert cells_a == [" a ", " b \\|"], f"input={line_a!r} result={cells_a!r}"
    # (b) inline code closes only on a backtick run of the same length
    line_b = "| `a ``` | b` | c |"
    cells_b = _split_cells(line_b)
    assert cells_b == [" `a ``` | b` ", " c "], f"input={line_b!r} result={cells_b!r}"


# MODEL_LIMITATIONS.md is intentionally not checked (human decision).
DOCS_EXCLUDED = {"MODEL_LIMITATIONS.md"}


def test_docs_cover_root_and_specification_markdown():
    found = {p.name for p in ROOT.glob("*.md")}
    found |= {f"Specification/{p.name}" for p in (ROOT / "Specification").glob("*.md")}
    missing = sorted(found - DOCS_EXCLUDED - set(DOCS))
    assert not missing, f"markdown files not in DOCS: {missing}"
