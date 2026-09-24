# Local eBOSS nine-mock fine weighted random n(z) (blinded)

This command sequence is for WSL in the existing
\`~/stress-energy-closure\` clone. It reads **only public observed and
realistic mock RANDOM FITS**. It never opens galaxy data, real odd
multipoles or any physical wake fit.

## What is fixed

The preceding protocol is
\`source_data/eboss_dr16_9mock_fine_weighted_nz_protocol_2026-09-24.json\`.
The redshift grid is forty fixed bins of \`0.01\` on \`[0.6,1.0)\`.
Every retained row uses the same previously validated product
\`WEIGHT_SYSTOT * WEIGHT_CP * WEIGHT_NOZ * WEIGHT_FKP\`.
Rows with numerical-zero systematic weight, defined as
\`WEIGHT_SYSTOT <= 1e-12\`, are excluded consistently. All four
observed-random FITS SHA256s and all 36 mock-random compressed FITS
SHA256s must agree with the retained source audits. Each ELG
\`chunk\` is compared **only** to the identically named chunk, and all
nonmatching labels are reported.

The first local run is intentionally a one-cap, one-realization
sanity check, not the fixed nine-mock result. Do not choose the final
mock set or redshift ranges from the observed odd sector.

## Install and retrieve the existing compact mock provenance

\`\`\`bash
cd ~/stress-energy-closure
git pull --ff-only origin main
source .venv/bin/activate
python -m pip install -r requirements.txt

python scripts/audit_eboss_dr16_9mock_fine_weighted_nz.py --self-test

mkdir -p eboss_workspace/local_nz/source
gh run download 36017670812 \
  --repo dvlahek/stress-energy-closure \
  --name eboss-dr16-nine-mock-window-ensemble \
  --dir eboss_workspace/local_nz/source
test -s eboss_workspace/local_nz/source/ninemock_window_ensemble.json
\`\`\`

The \`gh run download\` operation fetches the previously **completed
random-only** ensemble (small JSON and NPZ). The new local script
only needs its JSON to pin the mock SHA256s. If GitHub CLI is not
available, download that run's artifact through the GitHub website
and extract \`ninemock_window_ensemble.json\` to the indicated
\`source\` directory.

## First run: the SGC mock 0001

This reuses any existing SHA-verified observed SGC random FITS under
\`eboss_workspace/local_rr/fits\`. Missing random files are fetched
from their pinned public release URLs. The compressed mock FITS are
deleted after each complete verified input audit unless
\`--keep-mock-cache\` is explicitly requested.

\`\`\`bash
cd ~/stress-energy-closure
source .venv/bin/activate
set -o pipefail

python scripts/audit_eboss_dr16_9mock_fine_weighted_nz.py \
  --ensemble-json eboss_workspace/local_nz/source/ninemock_window_ensemble.json \
  --mock-ids 1 --caps SGC \
  --observed-cache-dir eboss_workspace/local_rr/fits \
  --mock-cache-dir eboss_workspace/local_nz/mock_fits \
  --out-dir eboss_workspace/local_nz \
  2>&1 | tee eboss_workspace/local_nz/pilot_SGC_0001.log
\`\`\`

Successful completion writes
\`eboss_workspace/local_nz/fine_weighted_nz_summary.json\` with status
\`local_fine_weighted_nz_shards_complete\`. Individual validated
catalogue summaries are kept as \`observed_SGC_LRG.json\`,
\`observed_SGC_ELG.json\`, \`mock_0001_SGC_LRG.json\`, and
\`mock_0001_SGC_ELG.json\`. No 36-file mock download is needed for
this pilot. The script checkpoints after each complete cap-ID pair.
Its \`--self-test\` also writes and reads its own synthetic ELG
FITS with two distinct chunk labels and an excluded tiny systematic
weight, without downloading any catalogue.

## Complete the already predeclared cohort

Once the pilot passed in the same environment:

\`\`\`bash
python scripts/audit_eboss_dr16_9mock_fine_weighted_nz.py \
  --ensemble-json eboss_workspace/local_nz/source/ninemock_window_ensemble.json \
  --mock-ids all --caps all \
  --observed-cache-dir eboss_workspace/local_rr/fits \
  --mock-cache-dir eboss_workspace/local_nz/mock_fits \
  --out-dir eboss_workspace/local_nz \
  2>&1 | tee eboss_workspace/local_nz/full_ninemock.log
\`\`\`

Previously checked input JSON records are reused only if their pinned
SHA256, cap, tracer, fixed bin edges, blinded status, and protocol
commit match. The nine mock IDs are unchanged:
\`0001, 0125, 0250, 0375, 0500, 0625, 0750, 0875, 1000\`.
With all 18 cap-ID cases, the status becomes
\`nine_mock_fine_weighted_nz_control_complete\`. This is a **random
selection audit**, not an exact LRG×ELG joint mask, a validated
physical window or the significance of an observed odd signal.

## Resource and error behavior

Before starting an expensive run, check \`free -h\` and \`df -h .\`.
Compressed FITS input can require substantially more RAM when Astropy
decodes its BINTABLE; the script processes one file at a time but does
not claim true streaming decompression of gzip FITS. A saved FITS whose
SHA256 does not match is rejected, not silently overwritten.
Failed inputs stop the run and retain previous completed JSON
records plus \`fine_weighted_nz_failure.json\` for inspection. The
final reports include unweighted and weighted forty-bin fractions,
per-tracer/cap total variation, and separate exactly matched ELG
chunk diagnostics, with both hemispheres and all nine mocks required
for the complete status.

The official LRG MANGLE and ELG DECaLS brick-mask selection is still
a separate gate:
\`docs/EBOSS_DR16_JOINT_MASK_SELECTION_GATE.md\`.
