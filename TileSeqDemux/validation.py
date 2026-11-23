"""Everything tileseq_demux checks before it does any work.

    ToolChecker           are bowtie2 and bowtie2-build installed and runnable?
    FastqFolderChecker    does each shared sample have exactly one R1 and one R2 file, and no clashing files?
    DemuxInputsValidator  runs those two checks and the SheetReader, and reports everything together
    InputValidationError  what DemuxInputsValidator raises; its message lists every problem by section

The command line calls DemuxInputsValidator first and does nothing else until it passes. The three checks are
independent: each one is a small object that can be replaced by a fake in a test, and the validator only gathers
what they report. The FASTQ folder is checked only when the sheet is valid, because it needs the sheet to know which
samples to look for.

Users never call this module. The command at the bottom is a developer check, for trying the validator on a real
sheet and folder before the command line exists (ticket 1.5, where `tileseq_demux --plan` replaces it):

    python -m TileSeqDemux.validation SHEET.csv FASTQ_DIR [--min-cover 0.9]
"""
import argparse
import itertools
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

from TileSeqDemux.rules import Rule, RuleSet
from TileSeqDemux.sheet import SheetError, SheetReader


class InputValidationError(Exception):
    """tileseq_demux cannot start. `sections` is a list of (heading, [problem, ...])."""

    HEADER = "tileseq_demux cannot start. Fix these and run it again:"

    def __init__(self, sections):
      self.sections = [(heading, list(problems)) for heading, problems in sections]

    def __str__(self):
        blocks = [self._format_section(heading, problems) for heading, problems in self.sections]
        return self.HEADER + "\n\n" + "\n\n".join(blocks)
    
    @staticmethod
    def _format_section(heading, problems):
        bullets = "\n".join("  - " + problem for problem in problems)
        return heading + "\n" + bullets


def _run_command(cmd, timeout=60):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=timeout)


def _last_line(text):
    """The last non-empty line of a tool's output. Warnings usually come first and the error last."""
    return next((line.strip() for line in reversed((text or "").splitlines()) if line.strip()), "no output")

@dataclass(frozen=True)
class ToolProbe:
    """What happened when one tool was looked for on PATH and asked for its version."""
    tool: str
    path: Optional[str]            # None when the tool is not on PATH
    status: Optional[int] = None   # its exit status, or None when it could not be run
    output: str = ""               # what it printed; its version is in here
    error: str = ""                # why it could not be run, or ""

    def details(self, context):
        return {"tool": self.tool, "path": self.path, "status": self.status, "error": self.error,
                "last_line": _last_line(self.output)}


TOOL_RULES = RuleSet(
    details=("tool", "path", "status", "error", "last_line"),
    rules=[
        Rule("not_found", "{tool} was not found on PATH.",
             lambda probe, context: probe.path is None),
        Rule("cannot_run", "{tool} ({path}) could not be run: {error}",
             lambda probe, context: probe.error != ""),
        Rule("bad_exit", "{tool} ({path}) exited with status {status} when asked for its version: {last_line}",
             lambda probe, context: probe.status != 0),
    ])

class ToolChecker:
    """Checks that the tools tileseq_demux calls are on PATH and run.

    Two jobs on the cluster failed this way during ticket 1.2: one ran under a Python without bowtie2 on its path.
    A job inherits the environment of the shell that submitted it, which is what the hint says.
    `run` and `which` can be replaced, which is how the tests avoid needing real tools."""

    HINT = ("Activate the environment that has these tools before running or submitting; "
            "an sbatch job inherits the environment of the shell that submitted it.")

    def __init__(self, tools=("bowtie2", "bowtie2-build"), run=None, which=None):
        self._tools = tuple(tools)
        self._run = run or _run_command
        self._which = which or shutil.which

    def check(self):
        """Return a list of problems, empty when every tool is usable."""
        problems = []
        for tool in self._tools:
            problems += TOOL_RULES.check(self._probe(tool), None).problems
        return problems

    def _probe(self, tool):
        """Look for `tool` on PATH and, if it is there, ask it for its version."""
        path = self._which(tool)
        if not path:
            return ToolProbe(tool, None)
        try:
            result = self._run([tool, "--version"])
        except (OSError, subprocess.TimeoutExpired) as err:
            return ToolProbe(tool, path, error=str(err))
        return ToolProbe(tool, path, result.returncode, result.stdout)


class FastqFolderChecker:
    """Checks the raw FASTQ folder against the parameter sheet. It only lists file names and never reads them.

    Each shared sample needs exactly one R1 and one R2 file. No file may already exist for a unique sample, because
    the demultiplexed output would collide with it. The folder is listed once and the names are grouped by sample ID,
    the first "_" token, which is how tileseq_mut reads them."""

    SUFFIX = ".fastq.gz"

    def check(self, sheet, folder):
        """Return a list of problems, empty when the folder is fine."""
        if not os.path.isdir(folder):
            return ["FASTQ folder %s does not exist" % folder]
        files = self._names_by_sample_id(folder)
        problems = []
        for shared, mate in itertools.product(sheet.shared_ids, (1, 2)):
            problems += self._mate_problems(files, shared, mate, folder)
        for unique in sheet.unique_ids:
            problems += self._clash_problems(files, unique, folder)
        return problems

    @staticmethod
    def _names_by_sample_id(folder):
        by_id = defaultdict(list)
        for name in sorted(os.listdir(folder)):
            if "_" in name:
                by_id[name.split("_")[0]].append(name)
        return by_id

    def _mate_problems(self, files, shared, mate, folder):
        marker = "_R%d_" % mate
        hits = [name for name in files.get(shared, []) if name.endswith(self.SUFFIX) and marker in name[len(shared) + 1:]]
        if len(hits) == 1:
            return []
        listed = (": " + ", ".join(hits)) if hits else ""
        return ["shared sample %s: expected one R%d file in %s, found %d%s" % (shared, mate, folder, len(hits), listed)]

    @staticmethod
    def _clash_problems(files, unique, folder):
        clash = files.get(unique, [])
        if not clash:
            return []
        return ["unique sample %s already has files in %s (%s); the demultiplexed output would collide with them"
                % (unique, folder, ", ".join(clash[:3]))]


class DemuxInputsValidator:
    """Runs every check and reports everything found in one message.

    Collaborators are passed in, so a test can replace any of them. `validate()` returns the ParameterSheet when all
    checks pass, so the caller never has to read the sheet a second time."""

    def __init__(self, sheet_path, fastq_dir, tool_checker=None, sheet_reader=None, folder_checker=None):
        self._sheet_path = sheet_path
        self._fastq_dir = fastq_dir
        self._tools = tool_checker or ToolChecker()
        self._reader = sheet_reader or SheetReader()
        self._folder = folder_checker or FastqFolderChecker()

    def validate(self):
        """Return the ParameterSheet, or raise InputValidationError listing every problem."""
        sheet, sheet_problems = self._read_sheet()
        checks = [("Tools", self._tool_problems()),
                  ("Parameter sheet (%s)" % self._sheet_path, sheet_problems),
                  ("Input files (%s)" % self._fastq_dir, self._folder_problems(sheet))]
        failed = [(heading, problems) for heading, problems in checks if problems]
        if failed:
            raise InputValidationError(failed)
        return sheet

    def _tool_problems(self):
        problems = self._tools.check()
        return problems + [ToolChecker.HINT] if problems else []

    def _read_sheet(self):
        """(sheet, []) when the sheet is valid, (None, problems) when it is not."""
        try:
            return self._reader.read(self._sheet_path), []
        except SheetError as err:
            return None, err.problems
        except OSError as err:
            return None, ["cannot read %s: %s" % (self._sheet_path, err.strerror or err)]

    def _folder_problems(self, sheet):
        if sheet is None:
            return ["not checked until the parameter sheet is valid"]
        return self._folder.check(sheet, self._fastq_dir)


def main(argv=None, tool_checker=None):
    parser = argparse.ArgumentParser(description="Developer check: validate a sheet and FASTQ folder and print the plan.")
    parser.add_argument("sheet")
    parser.add_argument("fastq_dir")
    parser.add_argument("--min-cover", type=float, help="override the sheet's Minimum coverage in the printed plan")
    args = parser.parse_args(argv)
    try:
      sheet = DemuxInputsValidator(args.sheet, args.fastq_dir, tool_checker=tool_checker).validate()
    except InputValidationError as err:
      print(err, file=sys.stderr)
      return 1
    for warning in sheet.warnings:
      print("warning:", warning)
    print("\n".join(sheet.describe(args.min_cover)))
    return 0


if __name__ == "__main__":
    sys.exit(main())