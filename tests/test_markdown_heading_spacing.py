"""User-facing Markdown documents: a heading must be preceded by a blank line."""

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

HEADING = re.compile(r"^#{1,6} ")
FENCE_OPEN = re.compile(r"^(`{3,}|~{3,})")
FENCE_RUN = re.compile(r"^(`+|~+)")


@pytest.mark.parametrize("rel", DOCS)
def test_heading_preceded_by_blank_line(rel):
    lines = (ROOT / rel).read_text(encoding="utf-8").splitlines()
    fence = None  # (char, length) of the open fence, or None outside a fence
    bad = []
    for i, line in enumerate(lines):
        stripped = line.lstrip()
        if fence is None:
            opener = FENCE_OPEN.match(stripped)
            if opener:
                fence = (opener.group(1)[0], len(opener.group(1)))
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
            continue
        if not HEADING.match(line):
            continue
        if i > 0 and lines[i - 1].strip() != "":
            bad.append(f"{rel}:{i + 1}: {line}")
    assert bad == [], "heading without a blank line before it:\n" + "\n".join(bad)
