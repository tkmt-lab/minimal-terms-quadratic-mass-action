#!/usr/bin/env python3
"""Maintainer command: record package bytes after an intentional edit."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED = {'.git', '.validation', 'work', '__pycache__'}
files = []
for path in sorted(ROOT.rglob('*')):
    rel = path.relative_to(ROOT)
    if any(part in EXCLUDED for part in rel.parts) or rel.as_posix() == 'manifest.json':
        continue
    if path.is_file():
        if path.is_symlink():
            raise ValueError(f'Symlink is not a package input: {rel}')
        data = path.read_bytes()
        files.append({'path': rel.as_posix(), 'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)})
out = {'schema': 'minimal-terms.manifest.v1', 'version': '1.0.0', 'files': files,
       'total_bytes_without_manifest': sum(row['bytes'] for row in files)}
(ROOT/'manifest.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({'files': len(files), 'bytes_without_manifest': out['total_bytes_without_manifest']}))
