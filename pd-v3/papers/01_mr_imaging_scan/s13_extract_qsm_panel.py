"""Stream all 28 officially mapped QSM/new T2* IDPs; retain all selected cis loci."""
import gzip
import io
import json
from pathlib import Path
import sys
import tarfile
import numpy as np
import pandas as pd
from common import ROOT
from s07_acquire_big40 import atomic_json, sha256, utc
from s08_acquire_extract_panel import chrom

DATA = ROOT / 'data/oxford_qsm'


class SequentialReader(io.RawIOBase):
    """Adapt tar's stream member without invoking its unsupported seek check."""
    def __init__(self, binary):
        super().__init__()
        self.binary = binary

    def readable(self):
        return True

    def seekable(self):
        return False

    def readinto(self, buffer):
        value = self.binary.read(len(buffer))
        buffer[:len(value)] = value
        return len(value)


def parse_member(binary, targets, meta, output):
    positions = {p for c, p in targets}
    expected = ['chr', 'rsid', 'pos', 'a1', 'a2', 'beta', 'se', 'pval(-log10)']
    rows, pieces, errors = 0, [], {}
    with io.TextIOWrapper(io.BufferedReader(SequentialReader(binary)), encoding='ascii') as handle:
        header = handle.readline().split()
        if header != expected:
            raise ValueError(f'Unexpected header: {meta["member_basename"]}: {header}')
        for chunk in pd.read_csv(handle, sep=r'\s+', names=expected, chunksize=300000,
                                 dtype={'chr': str, 'rsid': str, 'a1': str, 'a2': str}):
            ch = chrom(chunk['chr'])
            for key, mask in {
                'bad_coordinate': ~ch.isin([str(i) for i in range(1, 23)] + ['X']) | (chunk.pos <= 0),
                'bad_beta': ~np.isfinite(chunk.beta),
                'bad_se': ~np.isfinite(chunk.se) | (chunk.se <= 0),
                'bad_logp': ~np.isfinite(chunk['pval(-log10)']) | (chunk['pval(-log10)'] < 0),
                'bad_alleles': chunk.a1.isna() | chunk.a2.isna() | chunk.a1.eq(chunk.a2),
            }.items():
                errors[key] = errors.get(key, 0) + int(mask.sum())
            subset = chunk.loc[chunk.pos.isin(positions)].copy()
            subset['chrom_y'] = ch.loc[subset.index]
            subset = subset.loc[[p in targets for p in zip(subset.chrom_y, subset.pos)]]
            subset['source_data_row'] = subset.index.to_numpy() + 1
            pieces.append(subset)
            rows += len(chunk)
    if rows == 0 or any(errors.values()):
        raise ValueError(f'{meta["member_basename"]}: invalid rows {errors}')
    matched = pd.concat(pieces, ignore_index=True)
    out = pd.DataFrame({
        'phenotype': meta['phenotype'], 'trait': meta['trait'], 'family': meta['family'],
        'variant_y': matched.rsid, 'chrom_y': matched.chrom_y, 'pos_y': matched.pos,
        'ea_y': matched.a2, 'oa_y': matched.a1, 'beta_y': matched.beta, 'se_y': matched.se,
        'p_y': np.power(10., -matched['pval(-log10)']), 'neg_log10_p_y': matched['pval(-log10)'],
        'eaf_y': np.nan, 'n_y': np.nan, 'n_y_status': 'exact_IDP_sample_size_not_verified',
        'study_cohort_n_discovery': 19720, 'build_y': 'GRCh37',
        'build_evidence': 'inherited_official_BIG40_convention',
        'units_y': 'SD of Gaussianized and deconfounded imaging phenotype',
        'source_y': meta['source_y'], 'source_member': meta['member_basename'],
        'source_data_row': matched.source_data_row,
    })
    out.to_csv(output, index=False)
    return {'phenotype': meta['phenotype'], 'member': meta['member_basename'],
            'source_rows': rows, 'source_error_counts': errors, 'matched_rows': len(out),
            'matched_loci': len(set(zip(out.chrom_y, out.pos_y))), 'requested_loci': len(targets),
            'output': str(output.relative_to(ROOT)), 'output_sha256': sha256(output),
            'outcome_p_filter': None, 'allele_filter_at_extraction': None}


def main():
    manifest = json.loads((DATA / 'MANIFEST.json').read_text())
    archive = DATA / manifest['path']
    if archive.stat().st_size != manifest['size_bytes'] or sha256(archive) != manifest['sha256']:
        raise SystemExit('Archive does not match completed acquisition manifest')
    mapping_path = DATA / 'extracted/official_28_idp_mapping.csv'
    if sha256(mapping_path) != manifest['official_mapping_sha256']:
        raise SystemExit('Official phenotype mapping changed after acquisition')
    mapping = pd.read_csv(mapping_path)
    if len(mapping) != 28 or mapping.member_basename.duplicated().any():
        raise SystemExit('Expected 28 unique official phenotype members')
    instruments_path = ROOT / 'data/ukb_ppp/extracted/cis_primary_selected.csv'
    inst = pd.read_csv(instruments_path)
    targets = set(zip(chrom(inst.chrom), inst.pos37.astype(int)))
    destination = DATA / 'extracted/cis_loci_28'
    destination.mkdir(parents=True, exist_ok=True)
    metas = []
    for row in mapping.to_dict('records'):
        idp = row['member_basename'].removesuffix('.txt')
        metas.append(dict(phenotype=f'OxfordQSM2022:discovery:{idp}', trait=row['short_name'],
                          family=row['family'], member_basename=row['member_basename'], build_y='GRCh37',
                          n_y=np.nan, study_cohort_n_discovery=19720,
                          source_y=manifest['http']['url'],
                          units_y='SD of Gaussianized and deconfounded imaging phenotype',
                          metadata_source='https://www.fmrib.ox.ac.uk/ukbiobank/gwas_resources/index.html',
                          methods_source='https://pmc.ncbi.nlm.nih.gov/articles/PMC9174052/#Sec31'))
    plan = DATA / 'qsm_28_phenotypes.csv'
    pd.DataFrame(metas).to_csv(plan, index=False)
    meta_by_name = {m['member_basename']: m for m in metas}
    outputs, inventory = [], []
    started = utc()
    with gzip.open(archive, 'rb') as uncompressed:
        with tarfile.open(fileobj=uncompressed, mode='r|') as tar:
            for member in tar:
                if not member.isfile():
                    inventory.append({'name': member.name, 'bytes': member.size, 'type': 'nonfile'})
                    continue
                name = Path(member.name).name
                inventory.append({'name': member.name, 'bytes': member.size, 'officially_mapped': name in meta_by_name})
                stream = tar.extractfile(member)
                if name in meta_by_name:
                    if any(r['member'] == name for r in outputs):
                        raise ValueError(f'Duplicate mapped archive member {name}')
                    print(f'Extracting {name} {meta_by_name[name]["trait"]}...', flush=True)
                    output = destination / f'{name.removesuffix(".txt")}.cis_loci.csv'
                    outputs.append(parse_member(stream, targets, meta_by_name[name], output))
                    print(f'Complete {name}: {outputs[-1]["matched_rows"]} rows', flush=True)
                else:
                    while stream.read(1 << 20):
                        pass
        while uncompressed.read(1 << 20):
            pass
    if {r['member'] for r in outputs} != set(meta_by_name):
        raise SystemExit('Archive does not contain exactly the mapped requested family')
    combined = pd.concat([pd.read_csv(ROOT / r['output']) for r in outputs], ignore_index=True)
    output = DATA / 'extracted/qsm_28_outcomes.csv'
    combined.to_csv(output, index=False)
    record = dict(complete=True, started_at_utc=started, finished_at_utc=utc(),
                  n_phenotypes=len(outputs), source_archive_sha256=manifest['sha256'],
                  gzip_crc_eof_verified=True, outputs=outputs, inventory=inventory,
                  input_instruments_sha256=sha256(instruments_path), plan_sha256=sha256(plan),
                  output=str(output.relative_to(ROOT)), output_sha256=sha256(output),
                  script_sha256=sha256(__file__), command=sys.argv)
    atomic_json(DATA / 'QSM_28_EXTRACTION_MANIFEST.json', record)
    big40 = ROOT / 'data/big40'
    bm = json.loads((big40 / 'PHASE1_17_MANIFEST.json').read_text())
    bp, by = big40 / 'phase1_17_phenotypes.csv', big40 / 'extracted/phase1_17_outcomes.csv'
    if bm.get('complete') is not True or sha256(by) != bm['output_sha256']:
        raise SystemExit('BIG40 17-IDP extraction is incomplete or changed')
    both = ROOT / 'data/imaging_combined'
    (both / 'extracted').mkdir(parents=True, exist_ok=True)
    pd.concat([pd.read_csv(bp), pd.read_csv(plan)], ignore_index=True).to_csv(both / 'phenotypes_45.csv', index=False)
    all_y = both / 'extracted/outcomes_45.csv'
    pd.concat([pd.read_csv(by), combined], ignore_index=True).to_csv(all_y, index=False)
    atomic_json(both / 'MANIFEST.json', dict(complete=True, n_phenotypes=45,
        written_at_utc=utc(), output=str(all_y.relative_to(ROOT)), output_sha256=sha256(all_y),
        input_manifests={str(p.relative_to(ROOT)): sha256(p) for p in
                         [big40 / 'PHASE1_17_MANIFEST.json', DATA / 'QSM_28_EXTRACTION_MANIFEST.json']},
        phenotype_plan_sha256=sha256(both / 'phenotypes_45.csv'),
        note='17 BIG40 and 28 official QSM/new-T2* discovery IDPs; overlapping UKB sources'))
    print('Complete: 28 QSM/new-T2* and combined 45-IDP input plans', flush=True)


if __name__ == '__main__':
    main()
