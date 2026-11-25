"""The output folder of tileseq_demux: what a run will create, and how the folder goes from empty to complete.

    OutputPlan      every file and link a run will create, worked out before anything is written
    OutputFolder    the folder on disk: its state, the INCOMPLETE_RUN.txt marker, per-job commits and the links
    RunMode         what the user asked for: a fresh run, --resume or --force
    FolderState     what the folder already holds: nothing, an unfinished run, a complete run, or something else

The layout of a finished run:

    demux/
        97_S97_R1_001.fastq.gz ...     tile files, named the way tileseq_mut expects (naming.py)
        98_S98_R1_001.fastq.gz ...     links to the raw files of unshared samples and to Undetermined_*
        demux_manifest.tsv             written last
        demux_logs/                    demux_run.json, output_plan.json, jobs/<ID>/job.json, logs
        INCOMPLETE_RUN.txt             there while the run has not finished, removed as the very last step

While a job runs its tile files carry a .tmp suffix, which tileseq_mut ignores because it reads only *.fastq.gz.
commit_job renames a job's files once the job has passed its checks, and the job's own record (demux_logs/jobs/<ID>/
job.json) is written last: that file is the commit point. A file with a final name is therefore complete.

Nothing here deletes by pattern. Every deletion names an exact file from a plan, because .tmp is a common suffix and
other files may share the folder. Standard library only.
"""
import datetime
import itertools
import json
import os
import shutil
from collections import Counter
from dataclasses import asdict, dataclass
from enum import Enum
from typing import NamedTuple, Optional

from TileSeqDemux import naming

TMP_SUFFIX = ".tmp"
MARKER_NAME = "INCOMPLETE_RUN.txt"
MANIFEST_NAME = "demux_manifest.tsv"
LOGS_NAME = "demux_logs"
RUN_RECORD_NAME = "demux_run.json"
PLAN_NAME = "output_plan.json"
JOBS_NAME = "jobs"
JOB_RECORD_NAME = "job.json"

TILE = "tile"       # a file tileseq_demux writes
LINK = "link"       # a symbolic link to a raw file


class PlanError(Exception):
    """The output plan cannot be made from this sheet and raw folder."""


class OutputFolderError(Exception):
    """The output folder cannot be used as asked. The message says what to do."""


class AlreadyComplete(OutputFolderError):
    """The folder already holds a complete run, so --resume has nothing to do."""


class RunMode(Enum):
    NEW = "new"
    RESUME = "resume"
    FORCE = "force"


class FolderState(Enum):
    NEW = "new"                 # missing or empty
    INCOMPLETE = "incomplete"   # an unfinished run of ours
    COMPLETE = "complete"       # a finished run of ours
    FOREIGN = "foreign"         # something that is not ours, or not a folder

class FolderFacts(NamedTuple):
    """The three yes/no facts that settle what a folder holds."""
    missing_or_empty: bool    # there is nothing there yet
    made_by_us: bool          # it has a run record written by tileseq_demux
    finished: bool            # that run finished: it says complete, the manifest is there and the marker is gone


# Other combinations cannot happen, since a folder with a run record is not empty and a run cannot finish without a
# record. If one ever did, the safe reading is foreign: touch nothing.
_STATE_OF = {
    FolderFacts(missing_or_empty=True, made_by_us=False, finished=False): FolderState.NEW,
    FolderFacts(missing_or_empty=False, made_by_us=False, finished=False): FolderState.FOREIGN,
    FolderFacts(missing_or_empty=False, made_by_us=True, finished=False): FolderState.INCOMPLETE,
    FolderFacts(missing_or_empty=False, made_by_us=True, finished=True): FolderState.COMPLETE,
}

def _now():
    return datetime.datetime.now().isoformat(timespec="seconds")


@dataclass(frozen=True)
class PlannedFile:
    """One file or link of the plan. `name` is the final name inside the output folder."""
    kind: str
    name: str
    job: Optional[str] = None       # tile files: the shared sample ID of the job that writes them
    target: Optional[str] = None    # links: the absolute path of the raw file

    @property
    def tmp_name(self):
        return self.name + TMP_SUFFIX


class OutputPlan:
    """Every file and link a run will create. It is made from the sheet and a listing of the raw folder alone, so it
    can be printed (--dry-run) and tested without touching a disk."""

    def __init__(self, entries, raw_pairs):
        self.entries = tuple(entries)
        self.raw_pairs = {job: tuple(pair) for job, pair in raw_pairs.items()}   # job -> (raw R1 name, raw R2 name)
        self.jobs = tuple(self.raw_pairs)
        self._check_names_are_unique()

    @classmethod
    def build(cls, sheet, raw_dir, raw_names):
        """`raw_names` are the file names found in `raw_dir`. Raises PlanError if a shared sample has no single pair."""
        raw_dir = os.path.abspath(raw_dir)
        fastq_names = sorted(name for name in raw_names if name.endswith(naming.SUFFIX))
        raw_pairs = {group.shared: cls._raw_pair(group.shared, fastq_names) for group in sheet.groups}
        pooled = set(itertools.chain.from_iterable(raw_pairs.values()))
        tiles = itertools.chain.from_iterable(cls._tile_entries(group) for group in sheet.groups)
        links = [PlannedFile(LINK, name, target=os.path.join(raw_dir, name)) for name in fastq_names if name not in pooled]
        return cls(list(tiles) + links, raw_pairs)

    @staticmethod
    def _raw_pair(shared, fastq_names):
        pair = []
        for mate in (1, 2):
            hits = [name for name in fastq_names if naming.is_mate_file(name, shared, mate)]
            if len(hits) != 1:
                raise PlanError(f"shared sample {shared}: expected one R{mate} file in the raw folder, found {len(hits)}")
            pair.append(hits[0])
        return tuple(pair)

    @staticmethod
    def _tile_entries(group):
        for sample, mate in itertools.product(group.members, (1, 2)):
            yield PlannedFile(TILE, naming.tile_file_name(sample, mate), job=group.shared)

    def _check_names_are_unique(self):
        counts = Counter(entry.name for entry in self.entries)
        repeated = sorted(name for name, count in counts.items() if count > 1)
        if repeated:
            raise PlanError(f"these output names would be written twice: {', '.join(repeated)}")

    @property
    def tile_entries(self):
        return tuple(entry for entry in self.entries if entry.kind == TILE)

    @property
    def links(self):
        return tuple(entry for entry in self.entries if entry.kind == LINK)

    def tile_files(self, job):
        return tuple(entry for entry in self.tile_entries if entry.job == job)

    def describe(self):
        """Lines for --dry-run: what each job reads and writes, and how many links there will be."""
        lines = []
        for job in self.jobs:
            raw_r1, raw_r2 = self.raw_pairs[job]
            lines.append(f"job {job}: reads {raw_r1} and {raw_r2}")
            lines.append("  writes " + ", ".join(entry.name for entry in self.tile_files(job)))
        lines.append(f"{len(self.links)} links to raw files (unshared samples and other FASTQ files)")
        return lines

    def to_dict(self):
        return {"entries": [asdict(entry) for entry in self.entries],
                "raw_pairs": {job: list(pair) for job, pair in self.raw_pairs.items()}}

    @classmethod
    def from_dict(cls, data):
        return cls([PlannedFile(**item) for item in data["entries"]], data["raw_pairs"])


class OutputFolder:
    """The folder on disk. It knows what state the folder is in, what each run mode does in that state, and how a run
    moves it from empty to complete: prepare, then commit_job for each job, then finish. If a job fails, mark_failed."""

    def __init__(self, path, plan, rerun_hint="rerun the same command with --resume"):
        self.path = os.path.abspath(path)
        self.plan = plan
        self.rerun_hint = rerun_hint
        self.logs = os.path.join(self.path, LOGS_NAME)
        self.run_record = os.path.join(self.logs, RUN_RECORD_NAME)
        self.plan_record = os.path.join(self.logs, PLAN_NAME)
        self.marker = os.path.join(self.path, MARKER_NAME)
        self.manifest = os.path.join(self.path, MANIFEST_NAME)

    # -- what the folder holds ----------------------------------------------------------------------------------
    @property
    def state(self):
        made_by_us = os.path.isfile(self.run_record)
        facts = FolderFacts(missing_or_empty=self._is_missing_or_empty(), made_by_us=made_by_us,
                            finished=made_by_us and self._is_finished())
        return _STATE_OF.get(facts, FolderState.FOREIGN)

    def _is_missing_or_empty(self):
        return not os.path.exists(self.path) or (os.path.isdir(self.path) and not os.listdir(self.path))

    def _is_finished(self):
        """A run is finished when it says so, the manifest is there and the marker is gone."""
        return (self._read_json(self.run_record).get("status") == "complete"
                and os.path.isfile(self.manifest) and not os.path.exists(self.marker))

    def path_of(self, name):
        return os.path.join(self.path, name)

    def job_record_path(self, job):
        return os.path.join(self.logs, JOBS_NAME, job, JOB_RECORD_NAME)

    def completed_jobs(self):
        """The jobs of the plan that have a job record, which is written last."""
        return {job for job in self.plan.jobs if os.path.isfile(self.job_record_path(job))}

    # -- starting a run -----------------------------------------------------------------------------------------
    def prepare(self, mode):
        """Make the folder ready for a run of the given RunMode. Raises OutputFolderError when the folder cannot be
        used as asked, and AlreadyComplete when --resume finds nothing left to do."""
        self._refuse_if_raw_folder()
        self._actions()[(self.state, mode)]()

    def _actions(self):
        new, incomplete = FolderState.NEW, FolderState.INCOMPLETE
        complete, foreign = FolderState.COMPLETE, FolderState.FOREIGN
        return {
            (new, RunMode.NEW): self._start, (new, RunMode.RESUME): self._start, (new, RunMode.FORCE): self._start,
            (incomplete, RunMode.NEW): self._refuse_incomplete, (incomplete, RunMode.RESUME): self._resume,
            (incomplete, RunMode.FORCE): self._clear_and_start,
            (complete, RunMode.NEW): self._refuse_complete, (complete, RunMode.RESUME): self._already_complete,
            (complete, RunMode.FORCE): self._clear_and_start,
            (foreign, RunMode.NEW): self._refuse_foreign, (foreign, RunMode.RESUME): self._refuse_foreign,
            (foreign, RunMode.FORCE): self._refuse_foreign,
        }

    def _start(self):
        os.makedirs(self.logs, exist_ok=True)
        self._write_plan_record()
        self.update_run_record(status="running", started=_now())
        self._write_marker()

    def _resume(self):
        unfinished = [job for job in self.plan.jobs if job not in self.completed_jobs()]
        for job in unfinished:
            self._remove_tile_files(self.plan.tile_files(job), tmp=True, final=True)
        os.makedirs(self.logs, exist_ok=True)
        self._write_plan_record()
        resumes = self._read_json(self.run_record).get("resumes", 0) + 1
        self.update_run_record(status="running", resumed=_now(), resumes=resumes, failed_job=None, reason=None)
        self._write_marker()

    def _clear_and_start(self):
        self._clear_previous()
        self._start()

    def _refuse_incomplete(self):
        done, total = len(self.completed_jobs()), len(self.plan.jobs)
        raise OutputFolderError(f"{self.path} holds an unfinished run ({done} of {total} jobs complete; see {MARKER_NAME}). "
                                f"Rerun with --resume to finish the rest, or with --force to start over.")

    def _refuse_complete(self):
        raise OutputFolderError(f"{self.path} already holds a complete run. "
                                f"Use --force to start over, or choose another folder.")

    def _already_complete(self):
        raise AlreadyComplete(f"{self.path} already holds a complete run: nothing to do.")

    def _refuse_foreign(self):
        if os.path.exists(self.path) and not os.path.isdir(self.path):
            raise OutputFolderError(f"{self.path} exists and is not a folder.")
        raise OutputFolderError(f"{self.path} is not empty and was not made by tileseq_demux "
                                f"(it has no {LOGS_NAME}/{RUN_RECORD_NAME}), so nothing in it will be touched, "
                                f"not even with --force. Choose another folder or empty this one yourself.")

    def _refuse_if_raw_folder(self):
        here = os.path.realpath(self.path)
        inside = [entry.name for entry in self.plan.links if os.path.realpath(os.path.dirname(entry.target)) == here]
        if inside:
            raise OutputFolderError(f"{self.path} is the raw FASTQ folder. The output must be a different folder.")

    # -- during a run -------------------------------------------------------------------------------------------
    def update_run_record(self, **fields):
        """Add or change fields of demux_logs/demux_run.json. It is written whole and renamed into place."""
        record = self._read_json(self.run_record)
        record.update(fields)
        self._write_text(self.run_record, json.dumps(record, indent=2, sort_keys=True) + "\n")

    def commit_job(self, job):
        """Give a job's tile files their final names. Call this after the job has passed its checks, then write the
        job's record."""
        entries = self.plan.tile_files(job)
        missing = [entry.tmp_name for entry in entries if not os.path.isfile(self.path_of(entry.tmp_name))]
        if missing:
            raise OutputFolderError(f"job {job}: cannot commit, these files are missing: {', '.join(missing)}")
        for entry in entries:
            os.replace(self.path_of(entry.tmp_name), self.path_of(entry.name))

    def mark_failed(self, job, reason, keep_tmp=False):
        """Record that a job failed. The .tmp files go (unless keep_tmp), and so do any final-named files of the failing
        job, which can exist if it failed part way through its renames. Jobs that completed are left alone."""
        self._remove_tile_files(self.plan.tile_files(job), tmp=False, final=True)
        if not keep_tmp:
            self._remove_tile_files(self.plan.tile_entries, tmp=True, final=False)
        self.update_run_record(status="failed", failed_job=job, reason=reason, finished=_now(),
                               jobs_complete=len(self.completed_jobs()), jobs_total=len(self.plan.jobs))
        self._write_marker(failed_job=job, reason=reason)

    def write_report(self, name, text):
        """Write a report at the top of the folder (the manifest). It is written whole and renamed into place."""
        self._write_text(self.path_of(name), text)

    # -- ending a run -------------------------------------------------------------------------------------------
    def finish(self):
        """Link the raw files, mark the run complete and remove the marker, which is the very last step. The manifest
        must be written before this, and every job must have its record."""
        completed = self.completed_jobs()
        unfinished = [job for job in self.plan.jobs if job not in completed]
        if unfinished:
            raise OutputFolderError(f"cannot finish: no record for job(s) {', '.join(unfinished)}")
        self.link_raw()
        self.update_run_record(status="complete", finished=_now(),
                               jobs_complete=len(completed), jobs_total=len(self.plan.jobs))
        os.remove(self.marker)

    def link_raw(self):
        """Create the absolute symbolic links to the raw files. A link that already points to the right file is kept."""
        for entry in self.plan.links:
            path = self.path_of(entry.name)
            if os.path.islink(path) and os.readlink(path) == entry.target:
                continue
            if os.path.lexists(path):
                raise OutputFolderError(f"cannot link {entry.name}: a different file or link is already there")
            os.symlink(entry.target, path)

    # -- clearing, always by exact name -------------------------------------------------------------------------
    def _clear_previous(self):
        """Delete what a previous run made: the files of its recorded plan and of this plan, the manifest, the marker and
        demux_logs/. Anything else in the folder stays. Nothing is deleted if a link's name is held by a regular file."""
        both = itertools.chain(self._recorded_entries(), self.plan.entries)
        entries = {(entry.kind, entry.name): entry for entry in both}.values()
        # a name can be a tile file in one plan and a link in the other, and then it is deleted as a tile file
        tile_names = {entry.name for entry in entries if entry.kind == TILE}
        regular = [entry.name for entry in entries if entry.kind == LINK and entry.name not in tile_names
                   and os.path.exists(self.path_of(entry.name)) and not os.path.islink(self.path_of(entry.name))]
        if regular:
            raise OutputFolderError(f"{', '.join(regular)} are regular files, not links, "
                                    f"so they were left alone and nothing was deleted.")
        for entry in entries:
            self._remove_entry(entry)
        for path in (self.manifest, self.marker):
            self._remove_file(path)
        shutil.rmtree(self.logs, ignore_errors=True)

    def _remove_entry(self, entry):
        if entry.kind == LINK:
            if os.path.islink(self.path_of(entry.name)):
                os.remove(self.path_of(entry.name))
            return
        self._remove_tile_files([entry], tmp=True, final=True)

    def _remove_tile_files(self, entries, tmp, final):
        names = [entry.tmp_name for entry in entries if tmp] + [entry.name for entry in entries if final]
        for name in names:
            self._remove_file(self.path_of(name))

    @staticmethod
    def _remove_file(path):
        """Delete a regular file if it is there. Links and folders are never deleted here."""
        if os.path.isfile(path) and not os.path.islink(path):
            os.remove(path)

    # -- files of the folder itself -----------------------------------------------------------------------------
    def _write_plan_record(self):
        self._write_text(self.plan_record, json.dumps(self.plan.to_dict(), indent=2) + "\n")

    def _recorded_entries(self):
        """The entries of the plan an earlier run recorded, or none if there is no readable record."""
        try:
            return OutputPlan.from_dict(self._read_json(self.plan_record)).entries
        except (KeyError, TypeError, PlanError):
            return ()

    def _write_marker(self, failed_job=None, reason=None):
        lines = ["This run did not finish.", ""]
        if reason is not None:
            lines += [f"Job {failed_job} failed: {reason}",
                      f"{len(self.completed_jobs())} of {len(self.plan.jobs)} jobs finished and are complete.", ""]
        lines += [f"See {LOGS_NAME}/demux.log for details.",
                  f"Files ending in {TMP_SUFFIX} are incomplete and not usable. Files without it are complete.",
                  f"To finish: {self.rerun_hint}.",
                  "To start over, run it again with --force."]
        self._write_text(self.marker, "\n".join(lines) + "\n")

    @staticmethod
    def _write_text(path, text):
        temporary = path + TMP_SUFFIX
        with open(temporary, "w") as handle:
            handle.write(text)
        os.replace(temporary, path)

    @staticmethod
    def _read_json(path):
        try:
            with open(path) as handle:
                data = json.load(handle)
        except (OSError, ValueError):
            return {}
        return data if isinstance(data, dict) else {}