"""Command-line management of COD Offline, alongside the XRD Analyzer GUI.

Examples:
  python COD_OFFLINE_CLI.py elements --elements "Fe Cr Ni" --strict
  python COD_OFFLINE_CLI.py ids --ids "1000000 2002926"
  python COD_OFFLINE_CLI.py archive --file "D:\\COD\\cod-cifs-mysql.zip"
  python COD_OFFLINE_CLI.py full --confirm-large-download
"""
import argparse
import sys
from pathlib import Path
from datetime import datetime

from xrd_cod_offline import (
    CODSyncError, parse_elements, parse_ids, default_cod_path,
    sync_api, import_zip_archive, download_full_cod_zip,
)


def show_progress(value, message):
    print(f"{value:3d}% {message}" if value >= 0 else message, flush=True)


def main(argv=None) -> int:
    parser=argparse.ArgumentParser(description="XRD Analyzer COD Offline – 20.3.5 (CC0 metadata)")
    parser.add_argument('--db',type=Path,default=default_cod_path(),help='Ścieżka lokalnej bazy SQLite (COD .db3)')
    sub=parser.add_subparsers(dest='mode',required=True)
    p1=sub.add_parser('elements',help='Pobierz/odśwież struktury zawierające wybrane pierwiastki')
    p1.add_argument('--elements',required=True,help='np. Fe Cr Ni')
    p1.add_argument('--strict',action='store_true',help='Nie pobieraj rekordów o liczbie pierwiastków większej od liczby podanych')
    p2=sub.add_parser('ids',help='Pobierz/odśwież do 100 konkretnych rekordów COD')
    p2.add_argument('--ids',required=True,help='np. 1000000 2002926')
    p3=sub.add_parser('archive',help='Importuj lokalny ZIP COD zawierający pliki *.cif')
    p3.add_argument('--file',type=Path,required=True)
    p4=sub.add_parser('full',help='Ściągnij oficjalne archiwum COD i zaindeksuj offline; bardzo duży transfer')
    p4.add_argument('--confirm-large-download',action='store_true',help='Wymagane świadome potwierdzenie transferu wielu GB')
    p4.add_argument('--file',type=Path,help='Ścieżka zapisu pobranego ZIP')
    opts=parser.parse_args(argv)
    try:
        if opts.mode=='elements':
            result=sync_api(opts.db,elements=parse_elements(opts.elements),strict=opts.strict,progress=show_progress)
        elif opts.mode=='ids':
            result=sync_api(opts.db,ids=parse_ids(opts.ids),progress=show_progress)
        elif opts.mode=='archive':
            result=import_zip_archive(opts.db,opts.file,progress=show_progress)
        else:
            if not opts.confirm_large_download:
                parser.error('Pełne pobranie wymaga --confirm-large-download.')
            archive=opts.file or opts.db.parent / f"cod-cifs-mysql_{datetime.now():%Y%m%d}.zip"
            if not archive.exists():
                download_full_cod_zip(archive,progress=show_progress)
            result=import_zip_archive(opts.db,archive,progress=show_progress)
        print("ZAKOŃCZONO:")
        for k,v in result.items():print(f"  {k}: {v}")
        return 0
    except (CODSyncError,ValueError,OSError) as exc:
        print(f"BŁĄD COD Offline: {exc}",file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
