"""Construct the historical plan from source exposure/outcome extractions."""
import argparse
from common import ROOT, SNAPSHOT, RESULTS, save_csv
from pdv3.legacy import prepare_legacy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', type=type(ROOT), default=SNAPSHOT)
    args = parser.parse_args()
    inputs = [args.snapshot / 'results' / n for n in
              ('combined_exposure_instruments.csv', 'brain_sentinel_statistics.csv')]
    pairs, audit = prepare_legacy(args.snapshot)
    save_csv(pairs, RESULTS / 'legacy_planned_pairs.csv', inputs)
    save_csv(audit, RESULTS / 'legacy_input_audit.csv', inputs)
    print(f'{len(audit)} planned source-assay/outcome pairs; {len(pairs)} historical family pairs')


if __name__ == '__main__':
    main()
