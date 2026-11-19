"""Map every planned UKB assay to the complete published AGES accession table.

Read-only workbook extraction. No imaging result or historical-target preference
enters the mapping. A published accession is not verified regional availability.
"""
import csv
import hashlib
import json
from pathlib import Path
import re
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

import openpyxl

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'data/ages_catalog'
URL = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-27850-z/MediaObjects/41467_2021_27850_MOESM18_ESM.xlsx'
EXPECTED_SHA = '4799deaa25e03bc3bc95123ef247830c78c7c745caa81cac071e6b7366afdfbd'
NS = {'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows, fields):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def main():
    source = OUT/'download/41467_2021_27850_MOESM18_ESM.xlsx'
    source.parent.mkdir(parents=True, exist_ok=True)
    if not source.exists():
        with urllib.request.urlopen(URL, timeout=40) as response:
            data = response.read()
        assert hashlib.sha256(data).hexdigest() == EXPECTED_SHA
        source.write_bytes(data)
    assert sha(source) == EXPECTED_SHA
    wb = openpyxl.load_workbook(source, read_only=True, data_only=False)
    assert wb.sheetnames == ['Data 16']
    sheet = wb['Data 16']
    assert (sheet.max_row, sheet.max_column) == (4785, 4)
    headers = list(next(sheet.iter_rows(min_row=3, max_row=3, values_only=True)))
    assert headers == ['Study tag', 'Study Accession', 'Reported trait', 'Summary statistics file']
    observed = {}
    rows = []
    for number, cells in enumerate(sheet.iter_rows(min_row=4), 4):
        values = []
        for cell in cells:
            assert cell.data_type != 'f' and isinstance(cell.value, str)
            observed[cell.coordinate] = cell.value
            values.append(cell.value)
        assay, accession, trait, filename = values
        assert re.fullmatch(r'\d+_\d+_\d+', assay)
        assert re.fullmatch(r'GCST\d{8}', accession)
        prefix = 'Serum levels of protein '
        assert trait.startswith(prefix)
        rows.append(dict(ages_assay=assay, accession=accession,
                         reported_trait=trait, ages_target=trait.removeprefix(prefix),
                         published_summary_filename=filename, source_sheet='Data 16',
                         source_row=number, source_url=URL,
                         study_url='https://www.ebi.ac.uk/gwas/studies/'+accession,
                         access_status='published_accession_only_not_file_verified'))
    wb.close()
    assert len(rows) == len({r['accession'] for r in rows}) == len({r['ages_assay'] for r in rows}) == 4782
    # Independent ZIP/XML reader protects string identifiers and exact cell values.
    with zipfile.ZipFile(source) as archive:
        strings = [''.join(x.itertext()) for x in ET.fromstring(archive.read('xl/sharedStrings.xml')).findall('s:si', NS)]
        xml = ET.fromstring(archive.read('xl/worksheets/sheet1.xml'))
        checked = 0
        for cell in xml.findall('.//s:sheetData/s:row/s:c', NS):
            address = cell.attrib['r']
            if address not in observed:
                continue
            assert cell.find('s:f', NS) is None
            value = cell.find('s:v', NS)
            assert value is not None and cell.attrib['t'] == 's'
            assert strings[int(value.text)] == observed[address]
            checked += 1
    assert checked == 4782*4
    catalog = OUT/'all_4782_assays.csv'
    write_csv(catalog, rows, list(rows[0]))
    coverage_path = ROOT/'data/ukb_ppp/extracted/assay_coverage.csv'
    selected_path = ROOT/'data/ukb_ppp/extracted/cis_primary_selected.csv'
    ukb = list(csv.DictReader(coverage_path.open()))
    selected = {r['exposure_id'] for r in csv.DictReader(selected_path.open())}
    assert len(ukb) == 2941 and len({r['exposure_id'] for r in ukb}) == 2941
    lookup = {}
    for r in rows:
        lookup.setdefault(r['ages_target'], []).append(r)
    result = []
    matches = set()
    for r in ukb:
        found = lookup.get(r['gene'], [])
        base = dict(exposure_id=r['exposure_id'], gene=r['gene'], ukb_assay=r['assay'],
                    has_selected_cis_instrument=r['exposure_id'] in selected,
                    match_status='exact_reported_target_match' if found else 'no_exact_target_match',
                    matched_ages_assays=len(found), identity_scope='gene-label-only; cross-platform assay specificity unverified')
        if found:
            matches.add(r['exposure_id'])
            result.extend(dict(**base, **hit) for hit in found)
        else:
            result.append(dict(**base, **{k:'' for k in rows[0]}))
    mapping = OUT/'ukb_all_assays_to_ages.csv'
    write_csv(mapping, result, list(result[0]))
    assert {r['exposure_id'] for r in result} == {r['exposure_id'] for r in ukb}
    record = dict(complete=True, source_url=URL, source_sha256=EXPECTED_SHA,
                  source_sheet='Data 16', source_range='A4:D4785', source_rows=4782,
                  independent_cell_checks=checked, formula_cells=0,
                  planned_ukb_assays=2941, matched_ukb_assays=len(matches),
                  unmatched_ukb_assays=2941-len(matches),
                  matched_selected_cis_assays=len(matches & selected),
                  mapping_rows=len(result), imaging_results_used=False,
                  historical_gene_preference=False, input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [source,coverage_path,selected_path]},
                  output_sha256={str(p.relative_to(ROOT)):sha(p) for p in [catalog,mapping]},
                  code_sha256=sha(Path(__file__)), reader_version=openpyxl.__version__,
                  interpretation='Exact published gene labels identify potential external datasets, not equivalence of assays, instrument strength, verified data access or evidence for colocalization. Unmatched labels are retained; no fuzzy gene aliases or complex-target splitting.')
    (OUT/'MANIFEST.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k not in ['input_sha256','output_sha256','interpretation']},indent=2))


if __name__ == '__main__':
    main()
