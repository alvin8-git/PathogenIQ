import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "dl_validation",
    Path(__file__).resolve().parent.parent / "scripts" / "04_download_validation_data.py",
)
dl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dl)


SAMPLE = (
    "run_accession\tfastq_ftp\tfastq_md5\n"
    "ERR123\tftp.sra.ebi.ac.uk/vol1/ERR123_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/ERR123_2.fastq.gz\taaa;bbb\n"
    "ERR124\tftp.sra.ebi.ac.uk/vol1/ERR124.fastq.gz\tccc\n"
)


def test_parse_filereport_pairs_and_https():
    rows = dl._parse_filereport(SAMPLE)
    assert len(rows) == 3
    assert rows[0] == ("ERR123", "https://ftp.sra.ebi.ac.uk/vol1/ERR123_1.fastq.gz", "aaa")
    assert rows[1][2] == "bbb"   # second file of the paired run keeps its own md5
    assert rows[2] == ("ERR124", "https://ftp.sra.ebi.ac.uk/vol1/ERR124.fastq.gz", "ccc")
    assert all(u.startswith("https://") for _, u, _ in rows)


def test_parse_filereport_skips_header_only():
    assert dl._parse_filereport("run_accession\tfastq_ftp\tfastq_md5\n") == []


def test_filename_from_url_zenodo_content():
    # Zenodo: .../files/<name>/content -> the real name is before 'content'
    assert dl._filename_from_url(
        "https://zenodo.org/api/records/3632528/files/cami2_mouse_gut_gs.profile/content"
    ) == "cami2_mouse_gut_gs.profile"


def test_filename_from_url_plain():
    assert dl._filename_from_url(
        "https://ftp.sra.ebi.ac.uk/vol1/ERR588954_1.fastq.gz"
    ) == "ERR588954_1.fastq.gz"


def test_parse_filereport_keeps_existing_scheme():
    rows = dl._parse_filereport(
        "run_accession\tfastq_ftp\tfastq_md5\nERR9\thttps://example/x.fastq.gz\tzzz\n"
    )
    assert rows == [("ERR9", "https://example/x.fastq.gz", "zzz")]
