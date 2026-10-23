#!/usr/bin/env python3
"""Check the package inventory, exact witnesses, lower bounds and rejection tests."""
import sys
if not __debug__:
    raise SystemExit('Run without -O or PYTHONOPTIMIZE; checks use assert')

import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import re
import resource
import subprocess
import time

import exact_check as research
import face_proof
import independent_fraction as independent

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from regenerate import scaled_certificate, encoded


def read(path):
    return json.loads((ROOT/path).read_text())


def check_manifest():
    manifest = read('manifest.json')
    assert manifest['schema'] == 'minimal-terms.manifest.v1'
    assert manifest['version'] == '1.0.0'
    seen = set()
    total = 0
    for row in manifest['files']:
        path = Path(row['path'])
        assert not path.is_absolute() and '..' not in path.parts
        assert row['path'] not in seen
        seen.add(row['path'])
        target = ROOT/path
        assert target.is_file() and not target.is_symlink(), str(path)
        data = target.read_bytes()
        assert len(data) == row['bytes'], f'Byte count: {path}'
        assert hashlib.sha256(data).hexdigest() == row['sha256'], f'SHA-256: {path}'
        total += len(data)
    required = {'README.md', 'claims.json', 'CITATION.cff', '.zenodo.json', 'LICENSE',
                'Apache-2.0.txt', 'CC-BY-4.0.txt', 'AGENTS.md', '.github/workflows/check.yml',
                'paper/paper-v1.pdf', 'paper/paper-v1.tex', 'VERIFICATION-RECORD.md',
                'certificates/index.json', 'certificates/lower-bound.json',
                'certificates/independent-witnesses.json', 'checkers/exact_check.py',
                'checkers/face_proof.py', 'checkers/independent_fraction.py', 'checkers/verify.py',
                'tools/regenerate.py', 'tools/update_manifest.py', 'PROVENANCE.json'}
    for i in range(1, 10):
        required |= {f'certificates/{folder}/{i:02d}.json' for folder in ('witnesses','source-witnesses')}
    assert required <= seen
    assert total == manifest['total_bytes_without_manifest']
    assert total + (ROOT/'manifest.json').stat().st_size < 50_000_000
    return {'status': 'PASS', 'files': len(seen), 'bytes_without_manifest': total}


def check_paper_table(rows):
    tex = (ROOT/'paper/paper-v1.tex').read_text()
    entries = {}
    for line in tex.splitlines():
        # Only the nine-row parameter table has nine ampersand-separated columns.
        cells = [part.strip() for part in line.split('&')]
        if len(cells) == 9 and cells[0].isdigit() and 1 <= int(cells[0]) <= 9:
            values = [int(re.sub(r'\\+$', '', x).strip()) for x in cells[2:]]
            assert int(cells[0]) not in entries
            entries[int(cells[0])] = values
    assert len(entries) == 9
    for row in rows:
        assert entries[row['id']] == row['parameters'] + [row['minimum']]
    return {'status': 'PASS', 'parameter_rows': 9}


def check_independent():
    data = read('certificates/independent-witnesses.json')
    basis = independent.derive_basis()
    output = {}
    for key in ('witnesses', 'paper_controls'):
        output[key] = []
        for entry in data[key]:
            coeff = (independent.from_parameters(entry['parameters'], basis) if 'parameters' in entry
                     else independent.from_homogeneous(entry['ABCDEF']))
            assert coeff == independent.from_parameters([coeff[j] for j in independent.FREE], basis)
            output[key].append({'id': entry['id'],
                                **independent.verify_coefficients(coeff, entry['type'], entry['M'])})
    assert len(output['witnesses']) == 9 and len(output['paper_controls']) == 11
    output['excluded_controls'] = []
    for entry in data['excluded_controls']:
        try:
            independent.verify_coefficients(independent.from_homogeneous(entry['ABCDEF']))
        except AssertionError as err:
            assert str(err) == 'homogeneous common factor'
            output['excluded_controls'].append(entry['id'])
        else:
            raise AssertionError(f'Excluded control accepted: {entry["id"]}')
    assert len(output['excluded_controls']) == 2
    output['lower_bound'] = independent.lower_bound(basis)
    output['mutation_tests'] = independent.mutation_tests(basis, data['witnesses'])
    return output, basis


def rejection_tests(certs, lower):
    rejected = []

    def expect_rejected(name, obj, change, checker):
        bad = copy.deepcopy(obj)
        change(bad)
        try:
            checker(bad)
        except (AssertionError, ValueError, ZeroDivisionError):
            rejected.append(name)
        else:
            raise AssertionError(f'Corruption accepted: {name}')

    def set_key(key, value):
        return lambda obj: obj.__setitem__(key, value)

    expect_rejected('wrong-resultant', certs[0], set_key('resultant', '10'), research.witness)
    expect_rejected('wrong-term-count', certs[0], set_key('terms', 7), research.witness)
    expect_rejected('broken-equilibrium', certs[0], lambda c: c['a'].__setitem__(0, '1'), research.witness)
    expect_rejected('missing-real-root', certs[0], set_key('roots', []), research.witness)
    expect_rejected('wrong-multiplicity', certs[0], lambda c: c['roots'][0].__setitem__('multiplicity', 2), research.witness)
    expect_rejected('wrong-root-sign', certs[0], lambda c: c['roots'][0]['signs'].__setitem__('n', 1), research.witness)
    expect_rejected('wrong-reaction-rate', certs[0], lambda c: c['reactions'][0].__setitem__('rate', '4'), research.witness)
    expect_rejected('wrong-type', certs[0], set_key('code', 'P2H'), research.witness)
    expect_rejected('missing-infinity', certs[1], set_key('infinity_multiplicity', 0), research.witness)
    expect_rejected('missing-support-face', lower, lambda c: c['base_lower_six']['faces'].pop(), face_proof.lower_seven)
    expect_rejected('wrong-lower-bound', lower, set_key('bound', 7), face_proof.lower_seven)
    expect_rejected('missing-small-face', lower, lambda c: c['low_faces'].pop(), face_proof.lower_seven)
    optimized = subprocess.run([sys.executable, '-B', '-O', str(Path(__file__).resolve()), '--help'],
                               text=True, capture_output=True, timeout=15)
    assert optimized.returncode != 0 and 'Run without -O' in optimized.stderr
    return {'rejected': rejected, 'optimized_python_refused': True}


def check_claims(rows):
    claims = read('claims.json')
    assert claims['version'] == '1.0.0'
    assert claims['results'] == [{'type': r['type'], 'minimum': r['minimum'],
                                 'best_count_in_ye_v1': r['best_count_in_ye_v1']} for r in rows]
    assert claims['scope']['resultant_nonzero'] is True
    for key in ('reaction_minimum', 'historical_domain_equivalence', 'priority_claimed',
                'formal_proof_assistant', 'human_peer_review'):
        assert claims['scope'][key] is False
    assert len({c['id'] for c in claims['claims']}) == len(claims['claims'])
    for claim in claims['claims']:
        assert claim['kind'] in ('exact-computation', 'independent-implementation', 'prose')
        assert claim['evidence']
        for path in claim['evidence']:
            assert (ROOT/path).is_file(), path
    citation = read('CITATION.cff')  # JSON is valid YAML 1.2.
    zenodo = read('.zenodo.json')
    assert citation['cff-version'] == '1.2.0' and citation['version'] == zenodo['version'] == '1.0.0'
    assert zenodo['upload_type'] == 'publication' and zenodo['publication_type'] == 'preprint'
    assert citation['authors'][0]['family-names'] == 'Takimoto'
    return {'status': 'PASS', 'claims': len(claims['claims'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Optional detailed JSON output; stdout is a short summary')
    parser.add_argument('--witness-dir', type=Path, help='Check nine regenerated certificates from this directory')
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU, (120, 125))
    start = time.monotonic()
    report = {'schema': 'minimal-terms.verification.v1', 'inventory': check_manifest()}
    index = read('certificates/index.json')
    rows = index['witnesses']
    assert len(rows) == 9 and [r['id'] for r in rows] == list(range(1, 10))
    assert {r['type'] for r in rows} == research.NINE
    assert [r['minimum'] for r in rows] == [8]*7+[6]*2
    independent_output, basis = check_independent()
    report['independent'] = independent_output
    certs = []
    report['paper_witnesses'] = []
    for row in rows:
        path = (args.witness_dir/Path(row['certificate']).name if args.witness_dir else ROOT/row['certificate'])
        cert = json.loads(path.read_text())
        source = read(row['source_certificate'])
        research.witness(source)
        assert path.read_bytes() == encoded(scaled_certificate(source, row['positive_time_scale']))
        checked = research.witness(cert)
        assert checked['code'] == row['type'] and checked['terms'] == row['minimum']
        assert checked['resultant'] == row['resultant']
        coeff = independent.from_parameters(row['parameters'], basis)
        assert list(map(F, cert['a']+cert['b'])) == coeff
        assert [F(x)*row['positive_time_scale'] for x in source['theta_normalized']] == row['parameters']
        second = independent.verify_coefficients(coeff, row['ye_table_row'], row['minimum'])
        assert second['R'] == checked['resultant']
        assert independent_output['lower_bound']['minima'][str(second['type'])] == row['minimum']
        certs.append(cert)
        report['paper_witnesses'].append({'id': row['id'], 'research': checked, 'independent': second})
    lower = read('certificates/lower-bound.json')
    assert len(lower['low_faces']) == 8
    report['research_lower_bound'] = face_proof.lower_seven(lower)
    second_small = independent_output['lower_bound']['small_cells']
    independent_masks = set()
    for row in second_small:
        bits = [int(name in row['active']) for name in independent.NAMES]
        independent_masks.add((sum(bits[i]*2**i for i in range(3)),
                               sum(bits[i+3]*2**i for i in range(3))))
    assert independent_masks == {(r['xmask'],r['ymask']) for r in lower['low_faces']}
    controls = {entry['id']: entry for entry in independent_output['paper_controls']}
    assert controls['discussion-six']['type'] == 15 and controls['discussion-six']['M'] == 6
    for row in rows:
        previous = controls[f'table-{row["ye_table_row"]}']
        best = previous['M']
        if previous['type'] == controls['discussion-six']['type']:
            best = min(best, controls['discussion-six']['M'])
        assert best == row['best_count_in_ye_v1']
    assert sum(r['minimum'] < r['best_count_in_ye_v1'] for r in rows) == 5
    report['paper_transcription'] = check_paper_table(rows)
    report['claims'] = check_claims(rows)
    report['rejection_tests'] = rejection_tests(certs, lower)
    report['environment'] = {'python': platform.python_version(), 'platform': platform.platform()}
    report['elapsed_seconds'] = time.monotonic()-start
    report['peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
    report['status'] = 'PASS'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'paper_witnesses': 9, 'independent_witnesses_and_controls': 20,
                      'excluded_controls': 2, 'lower_bound_supports': {'research':49, 'independent':64, 'below_eight':8},
                      'corruptions_rejected': len(report['rejection_tests']['rejected']) +
                                              len(independent_output['mutation_tests']['rejected']),
                      'optimized_python_refused': True, 'elapsed_seconds': report['elapsed_seconds'],
                      'peak_rss_bytes': report['peak_rss_bytes']}))


if __name__ == '__main__':
    main()
