# E17D2b3b7 — lokalni read-only raw id/N gate (29. 9. 2026.)

**STATUS:** Source-only i sintetički GitHub Actions CI PASS. Stvarni numerički ASDF stupci još nisu pročitani; naredba ispod jest odobreni prvi lokalni pristup samo uz exact postojeću datoteku. [Zaključani prereg](../source_data/eboss_dr16_a03_e17d2b3b7_first_real_L1_id_N_columns_prereg_2026-09-29.json) Git blob `296bc22e3d38a414017cc018e1f8ba8658899278`, [čitač](../scripts/audit_eboss_dr16_a03_e17d2b3b7_first_real_L1_id_N.py) Git blob `f3df40e44a87aa0784fe1c017ce431da089920d8` i [source-only CI acceptance](../source_data/eboss_dr16_a03_e17d2b3b7_pinned_source_only_synthetic_CI_acceptance_2026-09-29.json).

## Što je odobreno i zašto

Samo prethodno lokalno verificirana `halo_info_000.asdf` iz `AbacusSummit_base_c000_ph000/halos/z0.950/halo_info/`. Stvarni YAML-only B6 report mora imati SHA256 `70820fb64965849ec692a035d4c2d6a1dedf38a8de492f19c88a4082dd42cedd`, YAML 14 416 B, 70 raw stupaca, 11 676 687 redaka po descriptorima. B7 otvara SAMO ASDF blok 0 (`id`, uint64) i 7 (`N`, uint32), unaprijed zaključane indekse 0, 5838343, 11676686 bez value-dependent selekcije, te punu N kolonu za deklarirane kontrole minimuma, maksimuma, zbroja i broja nula. Samo `N * ParticleMassHMsun` jest L1 **assigned catalog mass** u Msun/h. To NIJE pojedinačni M200c niti povijest istog haloa. Ovo NIJE fizički neutrinski wake niti statistički opservacijski rezultat. Bez drugog stvarnog haloa/PID/RV/cleaned/tree transfera, bez opaženog eBOSS odd i bez promjene E8/E16.

Sintetički [CI 36528381102](https://github.com/dvlahek/stress-energy-closure/actions/runs/36528381102) i [repo-integrity CI 36528385130](https://github.com/dvlahek/stress-energy-closure/actions/runs/36528385130) PASS na commitu `d2ba150cb48105402b2477a62dc9912aadc68882`. GitHub NIJE izvršio numeričko čitanje stvarnog ASDF-a.

## Jedan lokalni WSL blok

Radi samo u postojećem WSL Ubuntu terminalu, uz prethodno postojeće datoteke i .venv. Ispis je nebufferiran s `pipefail`. Ako ne postoji prethodno korišten lokalni Git-pinned ASDF dekoder, preuzima se *samo mala službena Python izvorna datoteka asdf.py* s točno zaključanoga commita i provjerava njezin Git blob. To NIJE transfer Abacus halo podataka.

```bash
set -euo pipefail
cd "$HOME/stress-energy-closure"
BRANCH="audit/eboss-elg-bit8-ra-orientation-20260925"
git fetch origin "$BRANCH"
git switch "$BRANCH"
git pull --ff-only origin "$BRANCH"

SCRIPT="scripts/audit_eboss_dr16_a03_e17d2b3b7_first_real_L1_id_N.py"
test "$(git hash-object "$SCRIPT")" = "f3df40e44a87aa0784fe1c017ce431da089920d8" || { echo 'STOP: B7 script drift'; exit 1; }
test "$(git hash-object source_data/eboss_dr16_a03_e17d2b3b7_first_real_L1_id_N_columns_prereg_2026-09-29.json)" = "296bc22e3d38a414017cc018e1f8ba8658899278" || { echo 'STOP: prereg drift'; exit 1; }
source .venv/bin/activate
DATA_DIR="$PWD/eboss_workspace/AbacusSummit_base_c000_ph000/halos/z0.950/halo_info"
B6_REPORT="$PWD/eboss_workspace/e17d2b3b6_reports/halo_info_000_schema_only.json"
OUT_DIR="$PWD/eboss_workspace/e17d2b3b7_reports"
mkdir -p "$OUT_DIR"

PINNED_DECODER_BLOB="8c5e5135736409bb0fab77bff06ad1259951189d"
EXT=""
while IFS= read -r candidate; do
  if [ "$(git hash-object "$candidate")" = "$PINNED_DECODER_BLOB" ]; then EXT="$candidate"; break; fi
done < <(find "$PWD" -type f -path '*/abacusnbody/data/asdf.py' -print)
if [ -z "$EXT" ]; then
  EXT="$PWD/eboss_workspace/e17d2b3b7_pinned_source/abacusnbody/data/asdf.py"
  mkdir -p "$(dirname "$EXT")"
  curl --fail --location --silent --show-error --max-time 45 --max-filesize 131072 \
    "https://raw.githubusercontent.com/abacusorg/abacusutils/24ab0dda5fea9ae406b1afdacaf4bbf989de9bc6/abacusnbody/data/asdf.py" \
    -o "$EXT.tmp"
  test "$(git hash-object "$EXT.tmp")" = "$PINNED_DECODER_BLOB" || { echo 'STOP: upstream decoder source mismatch'; exit 1; }
  mv "$EXT.tmp" "$EXT"
fi
test "$(git hash-object "$EXT")" = "$PINNED_DECODER_BLOB" || { echo 'STOP: decoder drift'; exit 1; }

echo "E17D2B3B7_PINNED_DECODER_READY $EXT"
(
  ulimit -v 6291456
  python -u "$SCRIPT" --self-test
  python -u "$SCRIPT" \
    --file "$DATA_DIR/halo_info_000.asdf" \
    --checksum-manifest "$DATA_DIR/checksums.crc32" \
    --b6-report "$B6_REPORT" \
    --official-extension "$EXT" \
    --output "$OUT_DIR/halo_info_000_id_N_locked_rows.json"
) 2>&1 | tee "$OUT_DIR/e17d2b3b7_full.log"
```

Ako bilo koji STOP ili greška nastupi, **ne** ponavljati uz promjenu datoteke, stupaca, odabira ili tolerancije. Sačuvati stdout i JSON ako je nastao. Očekivani završni marker je `E17D2B3B7_LOCAL_REAL_ID_N_ONLY_PASS`. Pošalji samo mali JSON i log, ne `halo_info_000.asdf`.
