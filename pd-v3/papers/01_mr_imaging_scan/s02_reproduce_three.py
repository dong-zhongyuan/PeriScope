"""Independently recompute all 576 legacy estimates, then inspect three proteins."""
import argparse
import json
import subprocess
import sys
import numpy as np
import pandas as pd
from common import ROOT, SNAPSHOT, RESULTS, provenance, save_csv
from pdv3.legacy import prepare_legacy
from pdv3.scan import scan_pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', type=type(ROOT), default=SNAPSHOT)
    parser.add_argument('--selftest', action='store_true')
    args = parser.parse_args()
    if args.selftest:
        subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(ROOT / 'tests'), '-v'],
                       check=True, cwd=ROOT)
        return
    inputs = [args.snapshot / 'results' / n for n in
              ('combined_exposure_instruments.csv', 'brain_sentinel_statistics.csv')]
    public_excerpt = ROOT / 'reference/sun2023_three_leads_verified_20261115.json'
    public = json.loads(public_excerpt.read_text())
    exposure = pd.read_csv(inputs[0]).set_index('assay')
    source_checks = []
    for record in public['records']:
        fields = record['raw_fields']
        row = exposure.loc[fields['UKBPPP ProteinID']]
        ident = fields['Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)'].split(':')
        exact = (row.rsid == fields['rsID'] and row.effect_allele == ident[3]
                 and row.other_allele == ident[2] and row.chrom == int(fields['CHROM'])
                 and row.pos == int(fields['GENPOS (hg38)']) and row.pos37 == int(ident[1]))
        errors = {key: float(row[key]) - float(fields[field]) for key, field in [
            ('beta_exposure', 'BETA (discovery, wrt. A1)'), ('se_exposure', 'SE (discovery)'),
            ('eaf_exposure', 'A1FREQ (discovery)')]}
        source_checks.append(dict(gene=record['protein'], exact_variant_alleles_coordinates=bool(exact),
                                  differences=errors, passed=bool(exact and all(abs(e) < 1e-12 for e in errors.values()))))
    output = RESULTS / 'public_st9_crosscheck.json'
    RESULTS.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({'source': public['source'], 'checks': source_checks,
                                 'portal_conflict': public['portal_conflict']}, indent=2) + '\n')
    provenance(output, [inputs[0], public_excerpt])
    if not all(c['passed'] for c in source_checks):
        raise SystemExit('Public ST9 versus legacy exposure discrepancy; see public_st9_crosscheck.json')
    pairs, input_audit = prepare_legacy(args.snapshot)
    full, audit = scan_pairs(pairs)
    full = full.merge(pairs[['gene', 'phenotype', 'trait', 'PAV_linked']].drop_duplicates(),
                      on=['gene', 'phenotype'], validate='one_to_one')
    save_csv(pairs, RESULTS / 'legacy_planned_pairs.csv', inputs)
    save_csv(input_audit, RESULTS / 'legacy_input_audit.csv', inputs)
    save_csv(full, RESULTS / 'legacy_full_family.csv', inputs,
             {'scope': 'historical 36 proteins x 16 IDPs; not proteome-wide scan'})
    save_csv(audit, RESULTS / 'legacy_harmonization_audit.csv', inputs)
    save_csv(full[full.gene.isin(['CD22', 'CD40', 'TFRC'])], RESULTS / 'three_proteins_all_48.csv', inputs,
             {'full_family': 'legacy_full_family.csv'})
    # Saved results enter only here, after independent computation and output.
    expected_path = args.snapshot / 'results/all_brain_MR_estimates.csv'
    expected = pd.read_csv(expected_path, float_precision='round_trip')
    comparison = full.merge(expected[['gene', 'outcome_accession', 'beta_MR', 'se_MR', 'p_MR', 'q_BH']],
                            left_on=['gene', 'phenotype'], right_on=['gene', 'outcome_accession'],
                            how='outer', validate='one_to_one', indicator=True)
    errors = {}
    for new, old in [('beta', 'beta_MR'), ('se', 'se_MR'), ('p', 'p_MR'), ('q_bh_global', 'q_BH')]:
        comparison[f'difference_{new}'] = comparison[new] - comparison[old]
        errors[new] = float(comparison[f'difference_{new}'].abs().max())
    complete = len(full) == 576 and full.gene.nunique() == 36 and full.phenotype.nunique() == 16
    passed = complete and comparison['_merge'].eq('both').all() and full.status.eq('estimated').all()
    passed = bool(passed and all(np.isfinite(v) and v < 1e-10 for v in errors.values()))
    save_csv(comparison, RESULTS / 'legacy_comparison.csv', inputs + [expected_path])
    summary = dict(status='passed' if passed else 'failed', tests=len(full), proteins=full.gene.nunique(),
                   phenotypes=full.phenotype.nunique(), max_absolute_error=errors,
                   coordinate_and_allele_checks=True,
                   scope='Independent arithmetic using preserved source extractions; not independent biological replication',
                   source_plan_pairs=len(input_audit), source_acquired_pairs=int(input_audit.acquisition_status.eq('acquired').sum()),
                   input_pair_status_counts=input_audit.status.value_counts().to_dict(),
                   covariance_assumption='First-order and delta SE assume zero exposure/outcome covariance; UKB overlap is possible')
    output = RESULTS / 'reproduction_summary.json'
    output.write_text(json.dumps(summary, indent=2) + '\n')
    provenance(output, inputs + [expected_path])
    print(json.dumps(summary, indent=2))
    if not passed:
        raise SystemExit('Legacy comparison failed; inspect comparison and audit before use')


if __name__ == '__main__':
    main()
