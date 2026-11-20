"""Check that the tools tileseq_demux needs are usable, and say plainly what is wrong when they are not.

A job inherits the environment of the shell that submitted it.
"""
import shutil
import subprocess

HINT = ("Activate the environment that has these tools before running or submitting; "
        "an sbatch job inherits the environment of the shell that submitted it.")
TSM_MARKERS = ("Unknown command", "there is no package called", "Execution halted", "command not found")


class ToolError(Exception):
    """One or more required tools are missing or broken. `problems` lists them all."""

    def __init__(self, problems):
        self.problems = list(problems)
        super().__init__("\n".join("  - " + p for p in self.problems) + "\n" + HINT)


def _run(cmd, timeout=60):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=timeout)


def _first_line(text):
    for line in (text or "").splitlines():
        if line.strip() and "libtinfo" not in line:          # the libtinfo warning is harmless noise on this cluster
            return line.strip()
    return "no output"


def _check_tsm(run, which):
    if not which("tsm"):
        return ["tsm was not found on PATH (tileseq_mut converts the sheet with `tsm csv2json`)."]
    try:
        res = run(["tsm", "csv2json", "-h"])
    except (OSError, subprocess.TimeoutExpired) as err:
        return ["tsm could not be run: %s" % err]
    marker = next((m for m in TSM_MARKERS if m in res.stdout), None)
    if res.returncode != 0 or marker:
        detail = next((ln.strip() for ln in res.stdout.splitlines() if marker and marker in ln), _first_line(res.stdout))
        return ["`tsm csv2json -h` does not work (%s). tsm needs R with the tileseqMave package visible on this node." % detail]
    return []


def check_tools(tools=("bowtie2", "bowtie2-build"), tsm=False, run=_run, which=shutil.which):
    """Return a list of problems (empty when everything is usable). `run` and `which` can be replaced in tests."""
    problems = []
    for tool in tools:
        path = which(tool)
        if not path:
            problems.append("%s was not found on PATH." % tool)
            continue
        try:
            res = run([tool, "--version"])
        except (OSError, subprocess.TimeoutExpired) as err:
            problems.append("%s (%s) could not be run: %s" % (tool, path, err))
            continue
        if res.returncode != 0:
            problems.append("%s (%s) exited with status %d when asked for its version: %s"
                            % (tool, path, res.returncode, _first_line(res.stdout)))
    if tsm:
        problems.extend(_check_tsm(run, which))
    return problems


def require_tools(tools=("bowtie2", "bowtie2-build"), tsm=False, run=_run, which=shutil.which):
    """Raise ToolError listing every problem, or return None."""
    problems = check_tools(tools, tsm, run, which)
    if problems:
        raise ToolError(problems)