#!/usr/bin/env python3
"""Verify immutable SHA manifest, replay 3 offline tests, verify unchanged outputs."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'E34_E36_SHA256_MANIFEST.json').read_text(encoding='utf-8'))

def check_hashes():
    for name,info in manifest['files'].items():
        p=ROOT/name
        if not p.is_file():
            raise RuntimeError('MISSING_PARENT_OR_OUTPUT: '+name)
        if p.stat().st_size!=info['bytes']:
            raise RuntimeError('BYTE_COUNT_MISMATCH: '+name)
        if hashlib.sha256(p.read_bytes()).hexdigest()!=info['sha256']:
            raise RuntimeError('SHA256_MISMATCH: '+name)
    print('E34_E36_MANIFEST_SHA_PASS files='+str(len(manifest['files'])),flush=True)

check_hashes()
for script in ('e34_offline_causal_history_qa.py','e35_offline_lowk_profile_bound_qa.py',
               'e36_offline_hidden_kinetic_initial_data_qa.py'):
    subprocess.run([sys.executable,str(ROOT/script)],cwd=ROOT,check=True)
check_hashes()
print('E34_E36_ALL_OFFLINE_REPLAY_PASS',flush=True)
