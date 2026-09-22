"""Agent adapters must route to one canonical testing procedure."""
from pathlib import Path

from tests.vacuity import assert_discriminates

ROOT = Path(__file__).resolve().parents[1]
CANONICAL = "docs/directives/testing.md"
ADAPTERS = (ROOT / ".agents/skills/vdg-test/SKILL.md",
            ROOT / ".claude/commands/vdg-test.md")

def _thin_testing_adapter(text):
    body = text.split("---", 2)[-1]
    return CANONICAL in body and len(body.splitlines()) <= 8

def test_testing_adapters_route_to_one_canonical_directive():
    assert_discriminates(
        _thin_testing_adapter, [path.read_text() for path in ADAPTERS],
        ["Read docs/testing.md.", CANONICAL + "\n" + "copied rule\n" * 12],
        "testing adapter must use the canonical path without copying its checklist")

def _thin_codex_adapter(text):
    return "CLAUDE.md" in text and "canonical shared" in text and len(text.splitlines()) < 60

def test_agents_routes_to_the_canonical_shared_guide():
    assert_discriminates(
        _thin_codex_adapter, [(ROOT / "AGENTS.md").read_text()],
        ["Use local rules.", "CLAUDE.md is canonical shared.\n" + "copied rule\n" * 70],
        "AGENTS.md must remain a thin adapter to CLAUDE.md")

def _within_budget(item):
    size, limit = item
    return size <= limit

def test_auto_loaded_instructions_stay_compact():
    budgets = {"CLAUDE.md": 7000, "AGENTS.md": 1200}
    budgets.update({str(path.relative_to(ROOT)): 3000
                    for path in (ROOT / "docs/directives").glob("*.md")})
    assert_discriminates(_within_budget,
                         [(len((ROOT / path).read_bytes()), limit)
                          for path, limit in budgets.items()],
                         [(7001, 7000), (3001, 3000)],
                         "auto-loaded instruction byte budget")
