"""COD Offline: bounded metadata API synchronization and local CIF archive import.

Only crystallographic *metadata* are indexed. No powder intensities are invented,
and none of the XRD scientific calculations are modified.

COD terms: CC0, credit the original crystallographic-data authors.
Official API: https://wiki.crystallography.net/RESTful_API/
"""
from __future__ import annotations

import csv
import io
import json
import math
import os
import re
import sqlite3
import time
import urllib.parse
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable

API_URL = "https://www.crystallography.net/cod/result"
CIF_BASE = "https://www.crystallography.net/cod/"
FULL_ARCHIVE_URL = "https://www.crystallography.net/archives/cod-cifs-mysql.zip"
USER_AGENT = "XRD-Analyzer-COD-Offline/20.3.5 (scientific research; rate-conscious)"
MAX_RESPONSE_BYTES = 48 * 1024 * 1024
MAX_CIF_BYTES = 6 * 1024 * 1024
MAX_RECORDS_PER_QUERY = 30000

# The exact column names used by the older XRD Analyzer CODDatabase connector.
TEXT_FIELDS = (
    "formula", "calcformula", "chemname", "mineral", "commonname", "sg", "sgHall",
    "doi", "journal", "method", "authors", "title", "volume", "issue",
    "firstpage", "lastpage", "compoundsource", "radiation", "radType",
    "radSymbol", "status", "flags", "date",
)
REAL_FIELDS = (
    "a", "b", "c", "alpha", "beta", "gamma", "vol", "celltemp",
    "diffrtemp", "cellpressure", "diffrpressure", "wavelength", "Rall",
    "Robs", "wRall", "wRobs", "gofall", "gofobs",
)
INTEGER_FIELDS = ("file", "sgNumber", "year", "nel")
DB_COLS = INTEGER_FIELDS + TEXT_FIELDS + REAL_FIELDS
REQUIRED_COD_FIELDS = {
    "file", "a", "b", "c", "alpha", "beta", "gamma", "sg", "sgNumber",
    "formula", "chemname", "mineral", "commonname", "doi", "journal",
    "year", "method", "nel",
}
_SYMBOLS = frozenset((
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
    "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
    "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr "
    "Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split())


class CODSyncError(RuntimeError):
    """Network, validation or database update issue."""


def default_cod_path() -> Path:
    home = os.getenv("LOCALAPPDATA")
    return Path(home) / "XRDAnalyzer" / "COD_Offline.db3" if home else Path.home() / ".xrd_analyzer" / "COD_Offline.db3"


def parse_elements(text: str) -> list[str]:
    """Accept e.g. 'Fe Cr Ni' and reject user-input URLs/SQL/formula expressions."""
    parts = re.split(r"[,\s;]+", str(text).strip())
    elements = [p[:1].upper() + p[1:].lower() for p in parts if p]
    if not elements or len(elements) > 8:
        raise ValueError("Podaj 1–8 symboli pierwiastków, np. Fe Cr Ni.")
    if len(set(elements)) != len(elements) or any(p not in _SYMBOLS for p in elements):
        raise ValueError("Niepoprawne lub powtórzone symbole pierwiastków.")
    return elements


def parse_ids(raw: str) -> list[int]:
    parts = re.split(r"[,\s;]+", str(raw).strip())
    parts = [p for p in parts if p]
    if not 1 <= len(parts) <= 100:
        raise ValueError("Podaj od 1 do 100 numerów COD.")
    if not all(re.fullmatch(r"\d{7}", p) for p in parts):
        raise ValueError("Każdy numer COD musi mieć dokładnie 7 cyfr.")
    out = [int(p) for p in parts]
    return list(dict.fromkeys(out))


def make_api_url(*, elements: Iterable[str] | None = None, ids: Iterable[int] | None = None,
                 strict: bool = False) -> str:
    params = [("format", "csv")]
    if (elements is None) == (ids is None):
        raise ValueError("Wybierz dokładnie jeden filtr: pierwiastki albo ID COD.")
    if elements is not None:
        symbols = parse_elements(" ".join(elements))
        params.extend((f"el{i}", symbol) for i, symbol in enumerate(symbols, 1))
        if strict:
            params.append(("strictmax", str(len(symbols))))
    else:
        numbers = parse_ids(" ".join(str(i) for i in ids))
        params.append(("id", ",".join(str(x) for x in numbers)))
    return API_URL + "?" + urllib.parse.urlencode(params)


def _num(value, integer: bool = False):
    if value is None or str(value).strip() in ("", "\\N", ".", "?", "None", "null"):
        return None
    s = str(value).strip()
    if s.startswith("+"):
        s = s[1:]
    # CIF numbers with standard uncertainty: 3.4412(5)
    match = re.fullmatch(r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)(?:\(\d+\))?", s)
    if not match:
        return None
    try:
        f = float(match.group(1))
        if not math.isfinite(f):
            return None
        return int(f) if integer and f.is_integer() else (None if integer else f)
    except (ValueError, OverflowError):
        return None


def _normalize_field(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).lower().lstrip("\ufeff"))


def _schema_candidates(values: dict) -> dict:
    # The COD CSV API is documented as a relational data export. Also accept
    # commonly occurring alternative header labels to preserve compatibility.
    lookup = {_normalize_field(key): val for key, val in values.items() if key is not None}
    alias = {
        "file": ["file", "id", "codid"], "sgNumber": ["sgNumber", "space_group_number", "sg_number"],
        "sg": ["sg", "spacegroup", "space_group"], "a": ["a"], "b": ["b"], "c": ["c"],
        "alpha": ["alpha", "alp"], "beta": ["beta", "bet"], "gamma": ["gamma", "gam"],
        "formula": ["formula", "chemical_formula_sum"], "nel": ["nel"],
    }
    out = {}
    for key in DB_COLS:
        candidates = alias.get(key, [key])
        raw = None
        for cand in candidates:
            norm = _normalize_field(cand)
            if norm in lookup:
                raw = lookup[norm]
                break
        if key in INTEGER_FIELDS:
            out[key] = _num(raw, integer=True)
        elif key in REAL_FIELDS:
            out[key] = _num(raw)
        else:
            out[key] = None if raw is None or str(raw).strip() in ("", ".", "?", "\\N") else str(raw).strip()[:4000]
    if out["nel"] is None:
        symbols = set(re.findall(r"[A-Z][a-z]?", out.get("formula") or "")) & _SYMBOLS
        out["nel"] = len(symbols) if symbols else None
    return out


def parse_cod_csv(data: bytes, *, max_records: int = MAX_RECORDS_PER_QUERY) -> list[dict]:
    if len(data) > MAX_RESPONSE_BYTES:
        raise CODSyncError("Odpowiedź COD przekracza limit bezpieczeństwa; zawęź filtr.")
    try:
        txt = data.decode("utf-8-sig")
    except UnicodeError:
        txt = data.decode("latin-1")
    if not txt.strip():
        raise CODSyncError("COD zwrócił pustą odpowiedź (nie zmieniono bazy).")
    if txt.lstrip()[:50].lower().startswith(("<!doctype", "<html", "<head", "<body")):
        raise CODSyncError("Serwer COD zwrócił HTML zamiast CSV; sprawdź filtry i spróbuj później.")
    snippet = txt[:16000]
    try:
        sep = csv.Sniffer().sniff(snippet, delimiters=",;\t").delimiter
    except csv.Error:
        sep = ","
    reader = csv.DictReader(io.StringIO(txt), delimiter=sep)
    headers = {_normalize_field(f) for f in (reader.fieldnames or []) if f}
    if "file" not in headers and "codid" not in headers and "id" not in headers:
        raise CODSyncError("COD CSV nie ma nagłówka identyfikatora (file/ID); import wstrzymany.")
    required = {"a", "b", "c", "alpha", "beta", "gamma"}
    if not required.issubset(headers):
        raise CODSyncError("COD CSV nie zawiera wymaganych parametrów sieci. Nie zapisano danych.")
    out = []
    for i, row in enumerate(reader):
        if i >= max_records:
            raise CODSyncError(f"Wynik przekracza {max_records:,} rekordów; zawęź filtrowanie.")
        if not row or all(not v for v in row.values()):
            continue
        values = _schema_candidates(row)
        if values["file"] is None or values["file"] < 1000000 or values["file"] > 9999999:
            continue
        if not all(values[k] is not None and values[k] > 0 for k in ("a", "b", "c")):
            continue
        if not all(values[k] is not None and 0 < values[k] <= 180 for k in ("alpha", "beta", "gamma")):
            continue
        out.append(values)
    return out


def _db_schema_sql() -> str:
    types = {}
    for col in DB_COLS:
        types[col] = "INTEGER" if col in INTEGER_FIELDS else "REAL" if col in REAL_FIELDS else "TEXT"
    return "CREATE TABLE IF NOT EXISTS data (" + ", ".join(
        f'"{col}" {types[col]}' + (" PRIMARY KEY" if col == "file" else "") for col in DB_COLS
    ) + ")"


def _ensure_db(conn: sqlite3.Connection) -> None:
    conn.execute(_db_schema_sql())
    have = {x[1] for x in conn.execute("PRAGMA table_info(data)")}
    if not REQUIRED_COD_FIELDS.issubset(have):
        raise CODSyncError("Istniejąca baza ma niezgodny schemat COD; wybierz nową ścieżkę.")
    conn.execute("CREATE TABLE IF NOT EXISTS cod_sync_log (id INTEGER PRIMARY KEY, kind TEXT, query TEXT, inserted INTEGER, skipped INTEGER, timestamp_utc TEXT)")
    conn.execute("CREATE INDEX IF NOT EXISTS cod_formula_idx ON data(formula)")
    conn.commit()


def _write_batch(conn: sqlite3.Connection, rows: list[dict]):
    if not rows:
        return
    cols = ",".join('"'+c+'"' for c in DB_COLS)
    marks = ",".join("?" for _ in DB_COLS)
    changes = ",".join(f'"{c}"=excluded."{c}"' for c in DB_COLS if c != "file")
    sql = f"INSERT INTO data ({cols}) VALUES ({marks}) ON CONFLICT(file) DO UPDATE SET {changes}"
    conn.executemany(sql, [tuple(row.get(c) for c in DB_COLS) for row in rows])


def _db_connect(path: Path):
    path = Path(path)
    if path.exists() and not path.is_file():
        raise CODSyncError("Docelowa ścieżka COD nie jest plikiem.")
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(str(path), timeout=45)
    try:
        connection.execute("PRAGMA busy_timeout=45000")
        _ensure_db(connection)
        return connection
    except Exception:
        connection.close()
        raise


def sync_api(db_path: Path, *, elements=None, ids=None, strict=False,
             opener: Callable | None = None, cancelled: Callable[[], bool] | None = None,
             progress: Callable[[int, str], None] | None = None) -> dict:
    """Synchronize a bounded subset. Rows replace matching COD IDs, leaving others intact."""
    url = make_api_url(elements=elements, ids=ids, strict=strict)
    if cancelled and cancelled():
        raise CODSyncError("Anulowano przed pobraniem.")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "text/csv"})
    getter = opener or urllib.request.urlopen
    if progress:
        progress(0, "Pobieranie metadanych z oficjalnego API COD...")
    try:
        response = getter(req, timeout=45)
        with response as remote:
            announced = remote.headers.get("Content-Length") if hasattr(remote, "headers") else None
            if announced and int(announced) > MAX_RESPONSE_BYTES:
                raise CODSyncError("Zbyt duża odpowiedź COD; użyj bardziej szczegółowego filtra.")
            chunks=[]; received=0
            while True:
                if cancelled and cancelled():
                    raise CODSyncError("Anulowano pobieranie; baza pozostaje bez zmian.")
                chunk=remote.read(131072)
                if not chunk: break
                received += len(chunk)
                if received > MAX_RESPONSE_BYTES:
                    raise CODSyncError("Przekroczono limit odpowiedzi COD (48 MB).")
                chunks.append(chunk)
                if progress and received % 1048576 < 131072:
                    progress(-1, f"Pobrano {received/1048576:.1f} MB...")
            data=b"".join(chunks)
    except CODSyncError:
        raise
    except Exception as exc:
        raise CODSyncError(f"Nie udało się pobrać COD: {type(exc).__name__}: {exc}") from exc
    rows=parse_cod_csv(data)
    if cancelled and cancelled():
        raise CODSyncError("Anulowano; baza nie została zmieniona.")
    if progress:
        progress(70, f"Wczytano {len(rows)} rekordów, zapis SQLite...")
    connection=_db_connect(db_path)
    try:
        with connection:
            _write_batch(connection, rows)
            connection.execute(
                "INSERT INTO cod_sync_log(kind,query,inserted,skipped,timestamp_utc) VALUES(?,?,?,?,?)",
                ("api", url, len(rows), 0, datetime.now(timezone.utc).isoformat())
            )
            total=connection.execute("SELECT COUNT(*) FROM data").fetchone()[0]
    finally:
        connection.close()
    if progress:
        progress(100, "Synchronizacja zakończona.")
    return {"imported": len(rows), "total": total, "path": str(Path(db_path).resolve()), "source": url}


CIF_TAGS = {
    "formula": ("_chemical_formula_sum",),
    "calcformula": ("_chemical_formula_structural",),
    "chemname": ("_chemical_name_systematic",),
    "commonname": ("_chemical_name_common",),
    "mineral": ("_chemical_name_mineral",),
    "sg": ("_space_group_name_H-M_alt", "_symmetry_space_group_name_H-M"),
    "sgNumber": ("_space_group_IT_number", "_symmetry_Int_Tables_number"),
    "sgHall": ("_space_group_name_Hall", "_symmetry_space_group_name_Hall"),
    "a": ("_cell_length_a",), "b": ("_cell_length_b",), "c": ("_cell_length_c",),
    "alpha": ("_cell_angle_alpha",), "beta": ("_cell_angle_beta",), "gamma": ("_cell_angle_gamma",),
    "vol": ("_cell_volume",), "doi": ("_journal_paper_doi",),
    "journal": ("_journal_name_full", "_journal_name_abbrev"),
    "year": ("_journal_year",), "title": ("_publ_section_title",),
    "method": ("_diffrn_measurement_method",),
    "celltemp": ("_cell_measurement_temperature",),
    "diffrtemp": ("_diffrn_ambient_temperature",),
}


def cif_to_row(content: bytes, cod_id: int) -> dict:
    """Parse selected CIF metadata via gemmi, using no atomic coordinates or diffraction intensities."""
    try:
        import gemmi
    except ImportError as exc:
        raise CODSyncError("Import CIF wymaga biblioteki gemmi (pip install gemmi).") from exc
    if len(content) > MAX_CIF_BYTES:
        raise CODSyncError(f"Pominięto zbyt duży plik CIF: COD {cod_id}.")
    try:
        document=gemmi.cif.read_string(content.decode("utf-8", errors="replace"))
        block=document.sole_block() if len(document)==1 else document[0]
        values={"file":cod_id}
        for key,tags in CIF_TAGS.items():
            for tag in tags:
                val=block.find_value(tag)
                if val and val not in ("?", "."):
                    # gemmi represents quoted strings with quotes; as_string removes them
                    values[key]=val.strip("'\"") if val.startswith(("'", '"')) else val
                    break
        return _schema_candidates(values)
    except CODSyncError:
        raise
    except Exception as exc:
        raise CODSyncError(f"Niepoprawny CIF COD {cod_id}: {type(exc).__name__}: {exc}") from exc


def import_zip_archive(db_path: Path, archive_path: Path,
                       cancelled: Callable[[], bool] | None = None,
                       progress: Callable[[int, str], None] | None = None,
                       batch_size: int = 250) -> dict:
    """Import the COD CIF tree from an official local .zip archive without extracting.

    Large bulk archives remain local. Partial import on cancel is possible and reported,
    and previous records are preserved. Per-file memory capped at 6 MB.
    """
    archive_path=Path(archive_path)
    if not archive_path.is_file() or not zipfile.is_zipfile(archive_path):
        raise CODSyncError("Wybierz poprawne archiwum ZIP zawierające pliki COD CIF.")
    conn=_db_connect(db_path)
    imported=0; skipped=0
    try:
        with zipfile.ZipFile(archive_path) as z:
            members = [x for x in z.infolist() if re.search(r"(?i)(?:^|/)(\d{7})\.cif$", x.filename.replace("\\", "/")) and not x.is_dir()]
            if not members:
                raise CODSyncError("W archiwum brak plików CIF o nazwach ID COD (1234567.cif).")
            if progress: progress(0, f"Import archiwum COD: {len(members):,} plików CIF")
            batch=[]
            for index, member in enumerate(members):
                if cancelled and cancelled():
                    raise CODSyncError(f"Anulowano: zapisano {imported} rekordów z dotychczasowych partii. Można wznowić import.")
                try:
                    cod_id = int(re.search(r"(\d{7})\.cif$", member.filename, re.I).group(1))
                    if member.file_size > MAX_CIF_BYTES:
                        skipped += 1; continue
                    with z.open(member, "r") as stream:
                        data=stream.read(MAX_CIF_BYTES + 1)
                    if len(data) > MAX_CIF_BYTES:
                        skipped += 1; continue
                    row=cif_to_row(data,cod_id)
                    if not all(row[k] is not None and row[k]>0 for k in ("a", "b", "c", "alpha", "beta", "gamma")):
                        skipped+=1; continue
                    batch.append(row)
                except (CODSyncError, zipfile.BadZipFile, OSError, KeyError):
                    skipped+=1; continue
                if len(batch) >= batch_size:
                    with conn:
                        _write_batch(conn,batch)
                    imported+=len(batch);batch.clear()
                if progress and (index % 250 == 0 or index==len(members)-1):
                    progress(int(100*(index+1)/len(members)), f"Pliki: {index+1:,}/{len(members):,}; zapisano: {imported:,}; pominięto: {skipped:,}")
            if batch:
                with conn: _write_batch(conn,batch)
                imported+=len(batch)
            with conn:
                conn.execute("INSERT INTO cod_sync_log(kind,query,inserted,skipped,timestamp_utc) VALUES(?,?,?,?,?)",
                             ("archive", archive_path.name, imported, skipped, datetime.now(timezone.utc).isoformat()))
            total=conn.execute("SELECT COUNT(*) FROM data").fetchone()[0]
    finally:
        conn.close()
    if progress: progress(100, "Import z archiwum zakończony.")
    return {"imported": imported, "skipped": skipped, "total": total,
            "path": str(Path(db_path).resolve()), "source": str(archive_path.resolve())}


def download_full_cod_zip(dest: Path, *, cancelled=None, progress=None,
                          opener=None, max_bytes: int = 20*1024*1024*1024) -> Path:
    """Explicit bulk opt-in only; staged .part download and atomic final rename.

    No zip extraction; see import_zip_archive. Do not poll COD in a loop.
    """
    dest=Path(dest)
    dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():
        raise CODSyncError("Plik docelowy już istnieje. Nie zostanie nadpisany.")
    tmp=dest.with_name(dest.name + ".part")
    if tmp.exists():
        raise CODSyncError("Istnieje wcześniejszy plik .part; usuń go po sprawdzeniu, aby zacząć od nowa.")
    request=urllib.request.Request(FULL_ARCHIVE_URL,headers={"User-Agent":USER_AGENT})
    getter=opener or urllib.request.urlopen
    size=0
    try:
        with getter(request,timeout=90) as src, tmp.open("xb") as out:
            while True:
                if cancelled and cancelled():
                    raise CODSyncError("Pobieranie anulowane.")
                chunk=src.read(1024*1024)
                if not chunk: break
                size += len(chunk)
                if size>max_bytes:
                    raise CODSyncError("Archiwum przekracza limit pobierania; przerwano.")
                out.write(chunk)
                if progress: progress(-1,f"Pobieranie archiwum COD: {size/1048576:.0f} MB")
        if not zipfile.is_zipfile(tmp):
            raise CODSyncError("Pobrany plik nie jest poprawnym archiwum ZIP.")
        os.replace(tmp,dest)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    return dest
