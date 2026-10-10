#!/usr/bin/env python3
"""Regression tests for the host-prep / pure-R boundary in module 06."""

import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
PURE_R = [
    HERE / "41_susie_coloc.R",
    HERE / "44_colocboost.R",
    HERE / "46_isotwas_train.R",
    HERE / "coloc_common.R",
]


def code_without_comments(path):
    lines = []
    for line in path.read_text().splitlines():
        # These module-06 workers do not use '#' inside string literals in a way
        # relevant to this static command-execution boundary check.
        lines.append(line.split("#", 1)[0])
    return "\n".join(lines)


def test_r_workers_do_not_execute_commands():
    banned = [
        r"\bsystem\s*\(",
        r"\bsystem2\s*\(",
        r"\bshell\s*\(",
        r"\bpipe\s*\(",
        r"\bfread\s*\([^\n]*\bcmd\s*=",
    ]
    for path in PURE_R:
        code = code_without_comments(path)
        for pattern in banned:
            assert re.search(pattern, code) is None, f"{path.name}: banned {pattern}"


def test_batch_workers_prepare_before_r():
    cases = [
        (HERE / "42a_run_coloc_shard.sh",
         "41_prepare_susie_coloc_inputs.py", "41_susie_coloc.R"),
        (HERE / "45a_run_colocboost_shard.sh",
         "44_prepare_colocboost_inputs.py", "44_colocboost.R"),
        (HERE / "47a_run_isotwas_shard.sh",
         "46_prepare_isotwas_inputs.py", "46_isotwas_train.R"),
    ]
    for worker, preparer, r_worker in cases:
        code = code_without_comments(worker)
        assert preparer in code
        assert r_worker in code
        assert code.index(preparer) < code.index(r_worker)
