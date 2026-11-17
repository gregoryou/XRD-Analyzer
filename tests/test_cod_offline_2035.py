"""Offline regression tests for COD 20.3.5: no live network required."""
import io
import sqlite3
from pathlib import Path
import zipfile

import pytest

from xrd_cod_offline import (
    parse_elements, parse_ids, make_api_url, parse_cod_csv, sync_api,
    cif_to_row, import_zip_archive, download_full_cod_zip, CODSyncError,
)

DATA = ("file,formula,chemname,mineral,commonname,sg,sgNumber,a,b,c,alpha,beta,gamma,doi,journal,year,method,nel\n"
        '1000000,"Ni O","Nickel oxide",,,"F m -3 m",225,4.17,4.17,4.17,90,90,90,"10.1111/test","J. Test",2020,"powder",2\n'
        '2002926,"Fe Ni","Fe-Ni alloy",,,"F m -3 m",225,3.54,3.54,3.54,90,90,90,,,2019,,2\n').encode('utf-8')


class FakeResponse(io.BytesIO):
    def __init__(self, data):
        super().__init__(data)
        self.headers = {"Content-Length": str(len(data))}


def opener_with(data):
    def opener(request, timeout):
        assert request.full_url.startswith("https://www.crystallography.net/cod/result?")
        assert timeout > 0
        return FakeResponse(data)
    return opener


def test_element_and_id_validation():
    assert parse_elements('fe cr NI') == ['Fe', 'Cr', 'Ni']
    assert parse_ids('1000000,2002926') == [1000000, 2002926]
    with pytest.raises(ValueError): parse_elements('Fe Fe')
    with pytest.raises(ValueError): parse_elements('DROP TABLE')
    with pytest.raises(ValueError): parse_ids('1234')
    assert 'strictmax=3' in make_api_url(elements=['Fe','Cr','Ni'], strict=True)
    assert 'id=1000000%2C2002926' in make_api_url(ids=[1000000,2002926])


def test_csv_parser_safety_and_metadata():
    rows = parse_cod_csv(DATA)
    assert len(rows) == 2
    assert rows[0]['file'] == 1000000
    assert rows[0]['a'] == 4.17
    assert rows[0]['formula'] == 'Ni O'
    assert rows[0]['sgNumber'] == 225
    assert rows[1]['doi'] is None
    with pytest.raises(CODSyncError): parse_cod_csv(b'<html>Too many records</html>')
    with pytest.raises(CODSyncError): parse_cod_csv(b'file,text\n1234567,Fe\n')
    with pytest.raises(CODSyncError): parse_cod_csv(DATA, max_records=1)


def test_live_sync_with_fake_official_server(tmp_path):
    dest = tmp_path / 'cod.db3'
    result = sync_api(dest, elements=['Fe','Ni'], strict=True, opener=opener_with(DATA))
    assert result['total'] == 2
    assert result['imported'] == 2
    result2 = sync_api(dest, ids=[2002926], opener=opener_with(DATA))
    assert result2['total'] == 2  # UPSERT does not duplicate
    with sqlite3.connect(dest) as db:
        assert db.execute("SELECT COUNT(*) FROM data").fetchone()[0] == 2
        assert db.execute("SELECT formula, sgNumber, doi FROM data WHERE file=1000000").fetchone() == ('Ni O', 225, '10.1111/test')
        assert db.execute("SELECT COUNT(*) FROM cod_sync_log").fetchone()[0] == 2
        fields = {r[1] for r in db.execute('PRAGMA table_info(data)')}
        assert {'file','a','b','c','alpha','beta','gamma','sg','sgNumber','formula','nel','authors','title'}.issubset(fields)


def test_failed_download_never_updates_db(tmp_path):
    dest=tmp_path/'cod.db3'
    sync_api(dest, elements=['Ni'], opener=opener_with(DATA))
    original=dest.read_bytes()
    with pytest.raises(CODSyncError):
        sync_api(dest, elements=['Ni'], opener=opener_with(b'<html>Rate limit</html>'))
    assert dest.read_bytes() == original
    with pytest.raises(CODSyncError):
        sync_api(tmp_path/'missing.db3',elements=['Ni'],opener=opener_with(DATA),cancelled=lambda: True)
    assert not (tmp_path/'missing.db3').exists()


CIF = b'''data_1000000
_chemical_formula_sum 'Ni O'
_space_group_name_H-M_alt 'F m -3 m'
_space_group_IT_number 225
_cell_length_a 4.172(2)
_cell_length_b 4.172
_cell_length_c 4.172
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_journal_year 2020
_journal_paper_doi '10.1234/xrd.test'
'''


def test_cif_archive_indexing_without_extraction(tmp_path):
    row = cif_to_row(CIF,1000000)
    assert row['a'] == 4.172
    assert row['formula'] == 'Ni O'
    archive=tmp_path/'fixture.zip'
    with zipfile.ZipFile(archive,'w') as z:
        z.writestr('cod/cif/1/00/00/1000000.cif',CIF)
        z.writestr('cod/cif/1/00/00/1000001.cif', b'bad cif')
        z.writestr('unrelated/README.txt','test')
    result=import_zip_archive(tmp_path/'COD.db3',archive)
    assert result['imported'] == 1
    assert result['skipped'] == 1
    assert result['total'] == 1
    # Repeat import to verify safe upsert.
    again=import_zip_archive(tmp_path/'COD.db3',archive)
    assert again['total']==1


def test_download_archive_rejects_wrong_content(tmp_path):
    path=tmp_path/'full.zip'
    with pytest.raises(CODSyncError):
        download_full_cod_zip(path, opener=lambda req, timeout: FakeResponse(b'not zip'))
    assert not path.exists() and not (tmp_path/'full.zip.part').exists()


def test_existing_cod_connector_reads_generated_sqlite(tmp_path):
    pytest.importorskip('PyQt6')
    from XRD_Analyzer_GUI import CODDatabase
    path=tmp_path/'auto_cod.db3'
    sync_api(path,ids=[1000000],opener=opener_with(DATA))
    reader=CODDatabase(path)
    try:
        assert reader.is_open
        assert reader.record_count==2
        assert reader.get_entry(1000000)[0] == 1000000
        assert reader.get_entry_details(1000000)['sgNumber'] == 225
    finally:
        reader.close()
