"""extract_affects_commits.py: ASCII-US/RS-delimited git-log extractor that
survives multi-line bodies with quotes and backslashes."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parents[1] / "extract_affects_commits.py"

# Import the helper directly for the in-process unit test.
sys.path.insert(0, str(SCRIPT.parent))
import extract_affects_commits as eac  # noqa: E402


@pytest.fixture
def tiny_repo(tmp_path):
    """A throw-away git repo with three commits including a tricky body."""
    repo = tmp_path / "repo"
    repo.mkdir()

    def git(*args, input_=None):
        return subprocess.run(
            ["git", *args], cwd=str(repo), capture_output=True, text=True,
            encoding="utf-8", input=input_, check=True,
        )

    git("init", "-q", "-b", "main")
    git("config", "user.email", "test@example.com")
    git("config", "user.name", "Test")
    git("config", "commit.gpgsign", "false")

    (repo / "a.txt").write_text("a")
    git("add", "a.txt")
    git("commit", "-q", "-m", "first commit")

    (repo / "b.txt").write_text("b")
    git("add", "b.txt")
    # Multi-line body with embedded quotes, a backslash, a JSON-like brace,
    # and an Affects: trailer. This is the case the fragile one-liner breaks on.
    tricky_body = (
        'fix(dmc) P1-1: Fisher "mid-p" convention\n\n'
        "Long body line with a backslash \\ and a brace { in it.\n"
        'Also a quote: "hello".\n\n'
        "Affects: lr@dmc_coverage, fisher@dmc_replicate\n"
    )
    git("commit", "-q", "-m", tricky_body)

    (repo / "c.txt").write_text("c")
    git("add", "c.txt")
    git("commit", "-q", "-m", "fix(dmr) P0-5: Stouffer math\n\nAffects: dmr_tile@dmr_coverage\n")

    return repo


def test_extract_round_trips_tricky_body(tiny_repo):
    commits = eac.extract("HEAD", repo_dir=tiny_repo)
    # Newest first.
    assert len(commits) == 3
    assert commits[0]["subject"].startswith("fix(dmr) P0-5")
    assert "Affects: dmr_tile@dmr_coverage" in commits[0]["body"]

    tricky = commits[1]
    assert tricky["subject"] == 'fix(dmc) P1-1: Fisher "mid-p" convention'
    assert "backslash \\" in tricky["body"]
    assert '"hello"' in tricky["body"]
    assert "Affects: lr@dmc_coverage, fisher@dmc_replicate" in tricky["body"]

    # Round-trips through json without escaping issues.
    payload = json.dumps(commits, ensure_ascii=False)
    decoded = json.loads(payload)
    assert decoded[1]["body"] == tricky["body"]


def test_cli_writes_file(tiny_repo, tmp_path):
    out = tmp_path / "commits.json"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "HEAD", "--out", str(out), "--repo", str(tiny_repo)],
        capture_output=True, text=True, encoding="utf-8", check=False,
    )
    assert result.returncode == 0, result.stderr
    data = json.loads(out.read_text(encoding="utf-8"))
    assert len(data) == 3
    assert any("P1-1" in c["subject"] for c in data)
