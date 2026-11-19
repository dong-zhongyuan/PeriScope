"""Run a complete supplied pair plan, retaining unestimable hypotheses."""
import argparse
import pandas as pd
from common import ROOT, RESULTS, save_csv
from pdv3.scan import scan_pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pairs', type=type(ROOT), required=True,
                        help='CSV of all planned instrument/outcome pairs, including missing outcomes')
    parser.add_argument('--output', type=type(ROOT), default=RESULTS / 'scan')
    parser.add_argument('--independent-instruments', action='store_true',
                        help='Attest that distinct instruments within each exposure are LD-independent')
    args = parser.parse_args()
    pairs = pd.read_csv(args.pairs, float_precision='round_trip')
    full, audit = scan_pairs(pairs, independent=args.independent_instruments)
    save_csv(full, args.output / 'full_family.csv', [args.pairs],
             {'independence_attested': args.independent_instruments})
    save_csv(audit, args.output / 'instrument_audit.csv', [args.pairs])
    print(f'{len(full)} planned hypotheses; {full.status.eq("estimated").sum()} estimated')


if __name__ == '__main__':
    main()
