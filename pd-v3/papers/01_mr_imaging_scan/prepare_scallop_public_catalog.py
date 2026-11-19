"""Source identity inventory for the original public regional follow-up.

No outcome results, selection, MR estimates or full pQTL downloads are used.
All UKB assays and all published protein identities remain in the inventory.
"""
import csv
from collections import Counter, defaultdict
import hashlib
import io
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile

import openpyxl

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'data/scallop_catalog'
SHEET = '1 Protein_info'
URL = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs42255-020-00287-2/MediaObjects/42255_2020_287_MOESM3_ESM.xlsx'
INPUTS = {
    'data/scallop_catalog/Supplementary_tables.xlsx': '928e2526e972469f52cfea931e1bcec1be9d699f3305db4a197aad0ffeb64d4d',
    'data/scallop_cd40/ZENODO_RECORD.json': '00d8a93e707fab32cd1c321e9fbc85b708f925490de6aa91c3347ce303e0cdb7',
    'data/ukb_ppp/extracted/assay_coverage.csv': '4f5a072ae5051be45f213e63d7378ffd4a0c149b3d1fa6afed651b5a695acba2',
}
NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
FIELDS = ['uniprot_reported', 'gene_reported', 'olink_assay_reported',
          'short_annotation', 'protein_name', 'cvd1_included']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_exact(path, content):
    payload = content.encode()
    if path.exists():
        if path.read_bytes() != payload:
            raise RuntimeError(f'Refusing to replace a different existing artifact: {path}')
    else:
        path.write_bytes(payload)


def write_csv(name, rows, fields):
    buffer = io.StringIO(newline='')
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    save_exact(OUT / name, buffer.getvalue())


def identity_relation(ukb, source):
    gene_equal = ukb['gene'] == source['gene_reported']
    raw = source['uniprot_reported']
    source_ids = set(raw.split('.')) if raw != 'NA' else set()
    ukb_ids = set(ukb['uniprot'].split(';'))
    if source_ids and source_ids == ukb_ids:
        if len(source_ids) > 1:
            return 'same_component_set_complex_assay_review_required'
        return 'same_reported_gene_and_uniprot' if gene_equal else 'same_uniprot_gene_review_required'
    if gene_equal:
        return 'same_reported_gene_only_identity_review_required'
    if source_ids & ukb_ids:
        return 'partial_complex_component_overlap_review_required'
    return None


def main():
    for name, expected in INPUTS.items():
        assert sha(ROOT / name) == expected, name
    source = OUT / 'Supplementary_tables.xlsx'
    workbook = openpyxl.load_workbook(source, read_only=True, data_only=False)
    assert workbook.sheetnames[0] == SHEET
    sheet = workbook[SHEET]
    assert [x.value for x in next(sheet.iter_rows(min_row=3, max_row=3, max_col=6))] == [
        'UniProt ID', 'Gene', 'Olink_assay', 'Short_annotation', 'Long name', 'CVDI']
    observed = {}; authors = []
    for number, cells in enumerate(sheet.iter_rows(min_row=4, max_row=95, max_col=6), 4):
        values = []
        for cell in cells:
            assert cell.data_type != 'f' and isinstance(cell.value, str)
            observed[cell.coordinate] = cell.value; values.append(cell.value)
        row = dict(zip(FIELDS, values))
        row.update(source_sheet=SHEET, source_row=number, source_url=URL)
        assert row['cvd1_included'] == 'Y'
        assert re.fullmatch(r'\d+_.+', row['olink_assay_reported'])
        authors.append(row)
    workbook.close()
    assert len(authors) == len({x['short_annotation'] for x in authors}) == 92
    # Independent raw worksheet XML check of every copied identity cell.
    with zipfile.ZipFile(source) as archive:
        strings = [''.join(x.itertext()) for x in ET.fromstring(archive.read('xl/sharedStrings.xml')).findall('s:si', NS)]
        tree = ET.fromstring(archive.read('xl/worksheets/sheet1.xml'))
        checked = set()
        for cell in tree.findall('.//s:sheetData/s:row/s:c', NS):
            address = cell.attrib['r']
            if address not in observed:
                continue
            assert cell.find('s:f', NS) is None
            assert cell.attrib['t'] == 's'
            assert strings[int(cell.find('s:v', NS).text)] == observed[address]
            checked.add(address)
        assert len(checked) == 92 * 6
    record = json.loads((ROOT / 'data/scallop_cd40/ZENODO_RECORD.json').read_text())
    assert record['id'] == 2615265 and len(record['files']) == 90
    by_name = {x['short_annotation']: x for x in authors}
    files = []
    for item in sorted(record['files'], key=lambda x: x['key']):
        assert item['key'].endswith('.txt.gz')
        token = item['key'][:-7]
        assert token in by_name
        assert re.fullmatch(r'md5:[0-9a-f]{32}', item['checksum'])
        assert item['size'] > 0
        url = item['links']['self']
        assert url == 'https://zenodo.org/api/records/2615265/files/' + item['key'] + '/content'
        files.append(dict(by_name[token], filename=item['key'], size_bytes=item['size'],
                          published_checksum=item['checksum'], download_url=url,
                          file_identity_match='exact_author_short_annotation',
                          file_status='published_file_inventory_not_newly_download_verified'))
    file_by_name = {x['short_annotation']: x for x in files}
    assert len(file_by_name) == 90
    missing = sorted(set(by_name) - set(file_by_name))
    assert missing == ['BNP', 'IL-4']
    # Author labels, not filename/gene similarity, establish these identities.
    assert file_by_name['TF']['gene_reported'] == 'F3'
    assert file_by_name['CD40-L']['gene_reported'] == 'CD40LG'
    assert file_by_name['CD40']['gene_reported'] == 'CD40'
    assert file_by_name['NT-pro_BNP']['uniprot_reported'] == 'NA'
    ukb = list(csv.DictReader((ROOT / 'data/ukb_ppp/extracted/assay_coverage.csv').open()))
    assert len(ukb) == len({x['exposure_id'] for x in ukb}) == 2941
    edges = []; coverage = []; counts = Counter()
    for u in ukb:
        matches = []
        for a in authors:
            relation = identity_relation(u, a)
            if relation is None:
                continue
            public = file_by_name.get(a['short_annotation'])
            edge = dict(exposure_id=u['exposure_id'], ukb_gene=u['gene'], ukb_uniprot=u['uniprot'],
                        ukb_assay=u['assay'], ukb_st3_original_row=u['st3_original_row'],
                        scallop_short_annotation=a['short_annotation'], scallop_gene=a['gene_reported'],
                        scallop_uniprot=a['uniprot_reported'], source_row=a['source_row'],
                        identity_relation=relation, public_file_listed=bool(public),
                        filename=public['filename'] if public else '',
                        assay_equivalence='not_established',
                        regional_readiness='not_assessed')
            edges.append(edge); matches.append(edge)
        strict = [x for x in matches if x['public_file_listed'] and x['identity_relation'] == 'same_reported_gene_and_uniprot']
        public_matches = [x for x in matches if x['public_file_listed']]
        status = ('same_gene_uniprot_public_file_listed' if strict else
                  'public_file_identity_review_required' if public_matches else
                  'author_identity_present_but_no_public_file' if matches else
                  'no_identity_match_in_this_source')
        counts[status] += 1
        coverage.append(dict(exposure_id=u['exposure_id'], gene=u['gene'], uniprot=u['uniprot'],
                             assay=u['assay'], has_cis_primary=u['has_cis_primary'],
                             status=status, n_author_matches=len(matches), n_public_file_links=len(public_matches),
                             n_same_gene_uniprot_public_files=len(strict),
                             matching_author_labels=json.dumps([x['scallop_short_annotation'] for x in matches]),
                             public_filenames=json.dumps([x['filename'] for x in public_matches]),
                             identity_relations=json.dumps([x['identity_relation'] for x in matches]),
                             assay_equivalence='not_established', regional_readiness='not_assessed'))
    write_csv('author_92_protein_identities.csv', authors, list(authors[0]))
    write_csv('all_90_public_files.csv', files, list(files[0]))
    write_csv('ukb_scallop_identity_links.csv', edges, list(edges[0]))
    write_csv('ukb_all_2941_assays.csv', coverage, list(coverage[0]))
    # Read back the complete output denominator and exact link references.
    readback = list(csv.DictReader((OUT / 'ukb_all_2941_assays.csv').open()))
    assert [x['exposure_id'] for x in readback] == [x['exposure_id'] for x in ukb]
    assert sum(int(x['n_author_matches']) for x in readback) == len(edges)
    assert sum(int(x['n_public_file_links']) for x in readback) == sum(x['public_file_listed'] for x in edges)
    for name, expected in INPUTS.items():
        assert sha(ROOT / name) == expected, name
    outputs = {name: sha(OUT / name) for name in [
        'author_92_protein_identities.csv', 'all_90_public_files.csv',
        'ukb_scallop_identity_links.csv', 'ukb_all_2941_assays.csv']}
    report = dict(complete=True, scope='public_source_identity_catalog_only',
                  author_identity_rows=len(authors), source_identity_cells_independently_checked=len(checked),
                  public_files=len(files), total_advertised_bytes=sum(x['size_bytes'] for x in files),
                  author_labels_without_public_file=missing, ukb_assays=len(ukb), identity_links=len(edges),
                  public_identity_links=sum(x['public_file_listed'] for x in edges),
                  coverage_status_counts=dict(counts),
                  public_source_targets_with_any_identity_link=len({x['scallop_short_annotation'] for x in edges if x['public_file_listed']}),
                  imaging_results_read=False, new_MR=False, candidate_nomination=False,
                  new_full_pqtl_files_downloaded=0, regional_colocalization_complete=False,
                  assay_equivalence_established=False, per_variant_N_or_sdY_assigned=False,
                  inputs_sha256=INPUTS, outputs_sha256=outputs,
                  source_code_sha256=sha(Path(__file__)), source_workbook_url=URL,
                  source_record_url='https://zenodo.org/records/2615265')
    save_exact(OUT / 'VERIFICATION.json', json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({k: v for k, v in report.items() if not k.endswith('sha256')}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
