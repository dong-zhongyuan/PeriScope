"""Freeze the complete official BIG40 catalog, before extracting associations."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pdv3.tabix_http import atomic, save_json


def main():
    source = ROOT / 'data/big40/download/IDPs_v4.html'
    spec = importlib.util.spec_from_file_location('big40_source', Path(__file__).with_name('s07_acquire_big40.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    catalog = module.parse_index(source)
    if [int(r['idp']) for r in catalog] != list(range(1, 3936)):
        raise ValueError('Official catalog no longer matches the original 3935-study release')
    dest = ROOT / 'data/big40_atlas'
    dest.mkdir(parents=True, exist_ok=True)
    for row in catalog:
        # Candidate mapping from the publication's contiguous accession range;
        # each trait name, assembly and cohort must be independently verified
        # against EBI metadata by s16 before it can produce statistics.
        number = 90002425 + int(row['idp'])
        accession = f'GCST{number}'
        low = (number - 1) // 1000 * 1000 + 1
        folder = f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{low}-GCST{low+999}/{accession}/harmonised/'
        row.update(accession=accession, phenotype=f'BIG40:release2:discovery:{row["idp"]}',
                   family=row['category'], is_qc=str(row['category_id'] == '0').lower(),
                   url=folder + accession + '.h.tsv.gz',
                   metadata_url=folder + accession + '.h.tsv.gz-meta.yaml',
                   accession_mapping_status='requires_exact_trait_metadata_check',
                   cohort='discovery', units_y='SD of deconfounded phenotype')
    import io
    buf = io.StringIO(newline='')
    writer = csv.DictWriter(buf, fieldnames=list(catalog[0]), lineterminator='\n')
    writer.writeheader(); writer.writerows(catalog)
    content = buf.getvalue().encode()
    atomic(dest / 'phenotypes_3935.csv', content)
    save_json(dest / 'PLAN.json', {
        'scope': 'Original first-paper whole BIG40 imaging catalog; no outcome-based selection',
        'phenotypes': len(catalog), 'biological_idps': sum(r['is_qc']=='false' for r in catalog),
        'qc_measures': sum(r['is_qc']=='true' for r in catalog),
        'assays': 2941, 'selected_cis_assays': 1954,
        'planned_assay_idp_pairs': 2941 * len(catalog),
        'planned_selected_cis_pairs': 1954 * len(catalog),
        'qsm_extension': '28 additional Oxford phenotypes already acquired; kept distinct from BIG40 IDs',
        'full_combined_catalog_size': len(catalog) + 28,
        'full_combined_planned_pairs': 2941 * (len(catalog) + 28),
        'source_catalog': str(source.relative_to(ROOT)),
        'source_catalog_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'plan_sha256': hashlib.sha256(content).hexdigest(),
        'publication': 'https://doi.org/10.1038/s41593-021-00826-4',
        'catalog_url': module.INDEX_URL,
        'complete_acquisition': False,
        'scope_note': 'All 3935 catalog entries retained; 16 QC measures explicitly marked, not called brain biological endpoints.',
        'independent_replication': False,
        'metadata_gate': 'Every accession must match exact normalized Oxford short name and GRCh38 harmonized metadata before extraction.'})
    print(json.dumps({'catalog_entries': len(catalog), 'planned_pairs': 2941 * len(catalog),
                      'plan': str(dest/'PLAN.json')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
