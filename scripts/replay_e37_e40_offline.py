import subprocess, sys
scripts=[
"e37_offline_lipschitz_history_qa.py",
"e38_offline_combined_certificate_qa.py",
"e39_offline_e28_highk_budget_qa.py",
"e40_offline_identifiability_qa.py",
]
for s in scripts:
    print("RUN", s, flush=True)
    subprocess.run([sys.executable,s],check=True)
print("E37_E40_REPLAY_PASS")
