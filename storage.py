"""SQLite checkpoints for scoring runs and an optional shared response cache."""

from contextlib import closing, contextmanager
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import tempfile

from filelock import FileLock


RUN_APPLICATION_ID = 0x41473252  # AG2R
CACHE_APPLICATION_ID = 0x41473243  # AG2C
SCHEMA_VERSION = 1
SQLITE_TIMEOUT = 30.0


def run_lock(directory: Path) -> FileLock:
    return FileLock(str(directory / "run.sqlite-lock"), timeout=0)


@dataclass(frozen=True)
class Outcome:
    status: str
    payload: str | None = None
    row_count: int = 0
    error: str = ""
    source: str = ""
    retrieved_at_utc: str = ""


@dataclass(frozen=True)
class CachedScore:
    payload: str
    row_count: int
    condition_tsv: str
    retrieved_at_utc: str = ""


@contextmanager
def _new_database(path: Path, application_id: int):
    """Publish a fully initialized database while the caller holds its creation lock."""
    if os.path.lexists(path):
        raise FileExistsError(f"Database already exists: {path}")
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}-", dir=path.parent)
    os.close(descriptor)
    temporary = Path(name)
    try:
        with closing(sqlite3.connect(temporary, timeout=SQLITE_TIMEOUT)) as db:
            db.execute("PRAGMA synchronous = FULL")
            with db:
                db.execute("BEGIN IMMEDIATE")
                db.execute(f"PRAGMA application_id = {application_id}")
                db.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
                yield db
        if os.path.lexists(path):
            raise FileExistsError(f"Database already exists: {path}")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _open_database(path: Path, application_id: int, *, readonly: bool = False) -> sqlite3.Connection:
    mode = "ro" if readonly else "rw"
    db = sqlite3.connect(path.absolute().as_uri() + f"?mode={mode}", uri=True, timeout=SQLITE_TIMEOUT)
    try:
        found_id = db.execute("PRAGMA application_id").fetchone()[0]
        found_version = db.execute("PRAGMA user_version").fetchone()[0]
        if (found_id, found_version) != (application_id, SCHEMA_VERSION):
            raise ValueError(f"Unknown AG2 database or unsupported schema: {path}")
        if not readonly:
            db.execute("PRAGMA synchronous = FULL")
    except BaseException:
        db.close()
        raise
    return db


class RunStore:
    """Store one run's immutable plan, request counter, and latest outcomes."""

    def __init__(self, directory: Path):
        path = directory / "run.sqlite"
        self._db = _open_database(path, RUN_APPLICATION_ID)
        try:
            record = self._db.execute("SELECT plan FROM execution WHERE singleton = 1").fetchone()
            if record is None:
                raise ValueError(f"Run database has no plan: {path}")
            plan = json.loads(record[0])
            if not isinstance(plan, dict):
                raise ValueError(f"Invalid run plan: {path}")
            self._plan = plan
        except BaseException:
            self._db.close()
            raise

    @classmethod
    def create(cls, directory: Path, plan: dict) -> "RunStore":
        encoded = json.dumps(plan, ensure_ascii=False, allow_nan=False)
        with _new_database(directory / "run.sqlite", RUN_APPLICATION_ID) as db:
            db.execute("""CREATE TABLE execution (
                singleton INTEGER PRIMARY KEY CHECK (singleton = 1),
                plan TEXT NOT NULL,
                attempted_calls INTEGER NOT NULL CHECK (attempted_calls >= 0)
            )""")
            db.execute("INSERT INTO execution VALUES (1, ?, 0)", (encoded,))
            db.execute("""CREATE TABLE outcomes (
                input_id TEXT NOT NULL,
                condition TEXT NOT NULL,
                scorer TEXT NOT NULL,
                status TEXT NOT NULL,
                payload TEXT,
                row_count INTEGER NOT NULL CHECK (row_count >= 0),
                error TEXT NOT NULL,
                source TEXT NOT NULL,
                retrieved_at_utc TEXT NOT NULL,
                PRIMARY KEY (input_id, condition, scorer)
            ) WITHOUT ROWID""")
            db.execute("""CREATE TABLE conditions (
                condition TEXT PRIMARY KEY,
                condition_tsv TEXT NOT NULL
            ) WITHOUT ROWID""")
        return cls(directory)

    @property
    def plan(self) -> dict:
        return self._plan

    @property
    def attempted_calls(self) -> int:
        return self._db.execute("SELECT attempted_calls FROM execution WHERE singleton = 1").fetchone()[0]

    def start_attempt(self) -> None:
        with self._db:
            self._db.execute("UPDATE execution SET attempted_calls = attempted_calls + 1 WHERE singleton = 1")

    @property
    def completed_outcomes(self) -> int:
        return self._db.execute("SELECT count(*) FROM outcomes WHERE status != 'api_error'").fetchone()[0]

    def get_group(self, input_id: str, condition: str) -> dict[str, Outcome]:
        rows = self._db.execute("""SELECT scorer, status, payload, row_count, error, source, retrieved_at_utc
            FROM outcomes WHERE input_id = ? AND condition = ?""", (input_id, condition))
        return {scorer: Outcome(*fields) for scorer, *fields in rows}

    def record(self, input_id: str, condition: str, outcomes: dict[str, Outcome]) -> None:
        with self._db:
            self._db.executemany("""INSERT INTO outcomes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(input_id, condition, scorer) DO UPDATE SET
                    status = excluded.status, payload = excluded.payload,
                    row_count = excluded.row_count, error = excluded.error, source = excluded.source,
                    retrieved_at_utc = excluded.retrieved_at_utc""", (
                (input_id, condition, scorer, result.status, result.payload, result.row_count,
                 result.error, result.source, result.retrieved_at_utc)
                for scorer, result in outcomes.items()
            ))

    def set_conditions(self, mapping: dict[str, str]) -> None:
        with self._db:
            self._db.executemany("""INSERT INTO conditions VALUES (?, ?)
                ON CONFLICT(condition) DO UPDATE SET condition_tsv = excluded.condition_tsv""", mapping.items())

    def get_conditions(self) -> dict[str, str]:
        return dict(self._db.execute("SELECT condition, condition_tsv FROM conditions"))

    def close(self) -> None:
        self._db.close()


class SharedCache:
    """Reuse immutable scorer responses between runs that select the same request."""

    def __init__(self, directory: Path, *, readonly: bool = False):
        self._directory = directory.resolve()
        self._readonly = readonly
        path = self._directory / "scores.sqlite"
        if readonly:
            self._db = _open_database(path, CACHE_APPLICATION_ID, readonly=True) if os.path.lexists(path) else None
            return

        self._directory.mkdir(parents=True, exist_ok=True)
        with FileLock(str(self._directory / "scores.sqlite-init.lock"), timeout=SQLITE_TIMEOUT):
            if not os.path.lexists(path):
                with _new_database(path, CACHE_APPLICATION_ID) as db:
                    db.execute("""CREATE TABLE scores (
                        group_key TEXT NOT NULL,
                        scorer_key TEXT NOT NULL,
                        payload TEXT NOT NULL,
                        row_count INTEGER NOT NULL CHECK (row_count > 0),
                        condition_tsv TEXT NOT NULL,
                        retrieved_at_utc TEXT NOT NULL,
                        PRIMARY KEY (group_key, scorer_key)
                    ) WITHOUT ROWID""")
            (self._directory / "locks").mkdir(exist_ok=True)
            self._db = _open_database(path, CACHE_APPLICATION_ID)

    def contains(self, group_key: str, scorer_key: str) -> bool:
        return self._db is not None and self._db.execute(
            "SELECT 1 FROM scores WHERE group_key = ? AND scorer_key = ?", (group_key, scorer_key)
        ).fetchone() is not None

    def lookup(self, group_key: str, scorer_key: str) -> CachedScore | None:
        if self._db is None:
            return None
        row = self._db.execute("""SELECT payload, row_count, condition_tsv, retrieved_at_utc FROM scores
            WHERE group_key = ? AND scorer_key = ?""", (group_key, scorer_key)).fetchone()
        return None if row is None else CachedScore(*row)

    def save(self, group_key: str, entries: dict[str, CachedScore]) -> None:
        if self._readonly or self._db is None:
            raise ValueError("Cannot write to a read-only cache.")
        with self._db:
            self._db.executemany("""INSERT INTO scores VALUES (?, ?, ?, ?, ?, ?)
                ON CONFLICT(group_key, scorer_key) DO NOTHING""", (
                (group_key, scorer_key, entry.payload, entry.row_count, entry.condition_tsv,
                 entry.retrieved_at_utc)
                for scorer_key, entry in entries.items()
            ))

    def lock(self, group_key: str) -> FileLock:
        if self._readonly:
            raise ValueError("Cannot lock a read-only cache.")
        name = hashlib.sha256(group_key.encode("utf-8")).hexdigest()
        return FileLock(str(self._directory / "locks" / f"{name}.lock"))

    def close(self) -> None:
        if self._db is not None:
            self._db.close()
