"""Stream-verify every NPY member of all six rebuilt inputs against the frozen run."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile


def verify(directory):
    expected = json.loads((Path(__file__).with_name('input_sources.json')).read_text())['expected_arrays']
    result = {}
    errors = []
    for name, entries in expected.items():
        path = Path(directory) / name
        if not path.is_file():
            errors.append(f'Missing {name}')
            continue
        checks = {}
        with zipfile.ZipFile(path) as archive:
            if sorted(archive.namelist()) != sorted(entries):
                errors.append(f'{name}: array inventory differs')
            for key, spec in entries.items():
                if key not in archive.namelist():
                    continue
                h = hashlib.sha256()
                with archive.open(key) as f:
                    for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
                        h.update(chunk)
                passed = archive.getinfo(key).file_size == spec['bytes'] and h.hexdigest() == spec['sha256']
                checks[key] = passed
                if not passed:
                    errors.append(f'{name}:{key} differs')
        result[name] = checks
        print(f'{name}: {sum(checks.values())}/{len(entries)} array hashes match', flush=True)
    if errors:
        raise RuntimeError('\n'.join(errors))
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inputs', type=Path, required=True)
    ap.add_argument('--report', type=Path)
    args = ap.parse_args()
    result = verify(args.inputs)
    if args.report:
        args.report.write_text(json.dumps({'all_arrays_match': True, 'files': result}, indent=2) + '\n')
