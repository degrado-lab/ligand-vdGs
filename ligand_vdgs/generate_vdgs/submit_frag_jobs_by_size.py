"""Submit fragment scripts in the order recorded by submission_order.tsv."""
import os
import re
import csv
import argparse
import subprocess

SLOTS_RE = re.compile(r'^#\$\s*-pe\s+smp\s+(\d+)', re.MULTILINE)

HEADER = ['order', 'script', 'fragment', 'slots', 'h_rt', 'tier_count']

def submission_rows(submit_dir):
    manifest = os.path.join(submit_dir, 'submission_order.tsv')
    with open(manifest, newline='') as handle:
        if not handle.readline().startswith('#'):
            raise SystemExit(f'[ERROR] missing manifest comment in {manifest}')
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != HEADER:
            raise SystemExit(f'[ERROR] unexpected manifest header in {manifest}: {reader.fieldnames}')
        rows = list(reader)
    if not rows:
        raise SystemExit(f'[ERROR] no jobs in {manifest}')
    seen = set()
    previous = None
    for index, row in enumerate(rows):
        try:
            order, slots, tier_count = (int(row[key]) for key in ('order', 'slots', 'tier_count'))
        except (ValueError, TypeError, KeyError) as exc:
            raise SystemExit(f'[ERROR] invalid manifest row {index} in {manifest}: {exc}') from exc
        path = row['script']
        if order != index or path in seen or not path.endswith('.sh'):
            raise SystemExit(f'[ERROR] invalid order or script at row {index} in {manifest}')
        seen.add(path)
        if previous is not None and (slots, tier_count) > previous:
            raise SystemExit(f'[ERROR] resource order increases at row {index} in {manifest}')
        previous = slots, tier_count
        try:
            with open(path) as script:
                match = SLOTS_RE.search(script.read())
        except OSError as exc:
            raise SystemExit(f'[ERROR] cannot read {path}: {exc}') from exc
        if not match or int(match.group(1)) != slots:
            raise SystemExit(f'[ERROR] slot request disagrees with manifest: {path}')
        yield slots, path

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('submit_dir', help='Directory containing submission_order.tsv.')
    parser.add_argument('--dry-run', action='store_true',
                        help="Print the submission order without calling qsub.")
    parser.add_argument('--print-order', action='store_true',
                        help='Print slot and script TSV for the production runner.')
    return parser.parse_args()

def main():
    args = parse_args()
    for slots, path in list(submission_rows(args.submit_dir)):
        if args.print_order:
            print(f'{slots}\t{path}')
        else:
            print(f'{"[DRY RUN] " if args.dry_run else ""}slots={slots} qsub {path}', flush=True)
        if not args.dry_run and not args.print_order:
            subprocess.run(['qsub', path], check=True)

if __name__ == '__main__':
    main()
