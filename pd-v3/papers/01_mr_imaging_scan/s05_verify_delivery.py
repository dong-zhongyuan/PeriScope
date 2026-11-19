"""Verify input snapshot and every recorded result/figure hash, without recomputing."""
import json
from common import ROOT, SNAPSHOT, PAPER
from pdv3.provenance import sha256


def main():
    manifest = json.loads((ROOT / 'reference/periscope_mr_20261115_manifest.json').read_text())
    failures = []
    for name, expected in manifest['sha256'].items():
        path = SNAPSHOT / name
        if not path.is_file() or sha256(path) != expected:
            failures.append(str(path))
    records = list(PAPER.rglob('*_provenance.json'))
    for record in records:
        data = json.loads(record.read_text())
        target = record.with_name(record.name.removesuffix('_provenance.json'))
        if not target.is_file() or sha256(target) != data['result_sha256']:
            failures.append(str(target))
        for group in ['code', 'inputs']:
            for rel, expected in data[group].items():
                path = ROOT / rel
                if not path.is_file() or sha256(path) != expected:
                    failures.append(f'{record.name}: changed {group} {rel}')
    if failures:
        raise SystemExit('\n'.join(failures))
    if not records:
        raise SystemExit('No result provenance found')
    print(f'PASS: {len(manifest["sha256"])} untouched snapshot files; {len(records)} artifact records and current code/input hashes')


if __name__ == '__main__':
    main()
