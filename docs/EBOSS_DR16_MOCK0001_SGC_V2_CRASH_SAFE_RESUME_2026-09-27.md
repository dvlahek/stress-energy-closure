# Mock0001 SGC V2: crash-safe nastavak nakon originalnog STOP-a

**Problem i status, 27. 9. 2026.** Izvorni korisnički `full_eligible_id0001_4800_48000_report.json` SHA256 `6731e0aa1f18a192b8067bcd44fda5aae8ffdb5164f64a86eab55d2b3e585f5c` ima 117014 bajtova. NGC s originalnih svih 7537 LRG i 17149 ELG galaksija završio je 4800R i 48000R. SGC se zaustavio na `Full-D sparse forward reverse four-term parity failed: R1D2`, bez per-term SGC histograma ili broja neusklađenih parova u trajnom checkpointu. Korisnik je potom ponovno pokrenuo **originalni** runner: priloženi stdout sadrži samo Python `ValueError: Original full-D output exists` (zaštita postojećeg JSON-a). Iz ta dva priložena artefakta **ne slijedi** da je Ubuntu kernel ili WSL VM ubijen zbog manjka memorije. WSL/OOM status zahtijeva zaseban kernel log ako se stvarno gasi cijela distribucija.

**NGC rezultat nije dokazana konvergencija:** `4800→48000R` na istim full-galaxy redcima daje relativni `L1(ξ)` 1,6749929559624868; originalni NGC `ξ` L1 pada sa 7,442272 na 3,802123. Najveća apsolutna šest-bin promjena `ℓ1` iznosi ~0,025465 i `ℓ3` ~0,039177. Iz spremljenih ponderiranih `6×24` parnih histograma moguće je zasebno provjeriti originalni DD i relativne normalizirane promjene DR, RD, RR: približno 0, 0,02434, 0,04173 i 0,04346. Puni DD NGC ima 735451 prihvaćenih parova, identičan za obje random gustoće. Taj rezultat služi isključivo numeričkoj dijagnostici jednog mocka; ne daje eBOSS survey-window certifikat niti EV wake amplitude.

[Originalni v1 runner](../scripts/run_eboss_dr16_full_eligible_mock0001_sparse.py), njegove originalne JSON/log bajtove i SGC failure **ne prepisivati**. [Zasebni post-failure v2 izvorni protokol](../source_data/eboss_dr16_full_eligible_mock0001_sgc_resume_v2_after_failure_protocol_2026-09-27.json) eksplicitno priznaje da uzrok originalne SGC reverse-parity razlike nije izmjeren. [V2 runnable](../scripts/resume_eboss_dr16_full_eligible_mock0001_sgc_v2.py) koristi originalni NGC kao nepromjenjiv SHA-pinned parent te računa samo SGC; čak i pri ponovnom WSL prekidu atomarno ostavlja svaki dovršeni SGC parni histogram. Originalni A-02 neovisni forward/reverse dense 600D/1200R weighted hist SHA i sparse↔dense preflight i dalje se obavezno reproduciraju. Produkcijski veliki SGC reverse član dobiva se istovjetnim geometrijskim zrcaljenjem već izračunatog forward člana: predznak `mu` se mijenja, `D1R2` i `R1D2` zamijene mjesta, a svaka parna normalizacija prati pripadni član. To je matematička definicija iste reverse geometrije, **ne neovisno ponovno brojanje**; v1 originalna neovisna anomalija ostaje otvorena.

## Jedna sigurna WSL naredba: samo preostali SGC 4800R

Novi runner odbija promijenjeni originalni v1 JSON i prije novih FITS redaka ponovno u cijelosti provjerava svih 72 izvorna gzip SHA (1860198719 bajtova). Zahtijeva istovjetan SGC source SHA, originalni 600D/1200R i E4 4800R hash. Prvi lokalni korak namjerno **ne pokreće skupi SGC 48000R**. Nakon svake SGC pair statistike automatski sprema zaseban v2 output. Drugi poziv iste naredbe legalno nastavlja, bez brisanja checkpointa, ako se prethodni WSL proces prekinuo.

```bash
set -Eeuo pipefail
cd "$HOME/stress-energy-closure"
BRANCH="audit/eboss-elg-bit8-ra-orientation-20260925"
git fetch origin "$BRANCH"
git switch "$BRANCH"
git pull --ff-only origin "$BRANCH"
source .venv/bin/activate
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MALLOC_ARENA_MAX=2 PYTHONUNBUFFERED=1 PYTHONFAULTHANDLER=1
mkdir -p eboss_workspace/a03_full_eligible_mock
python -u scripts/resume_eboss_dr16_full_eligible_mock0001_sgc_v2.py --self-test
python -u scripts/resume_eboss_dr16_full_eligible_mock0001_sgc_v2.py --run \
  2>&1 | tee -a eboss_workspace/a03_full_eligible_mock/full_eligible_id0001_4800_48000_resume_v2_stdout.log
```

Nakon ovog *staged* koraka v2 ispisuje `V2_SGC_4800_STAGED_STOP_SAVED`, a novi izvještaj je `eboss_workspace/a03_full_eligible_mock/full_eligible_id0001_4800_48000_resume_v2_report.json`. Zajedno s pripadnim `...resume_v2_stdout.log` može se provjeriti bez daljnjeg velikog runa. To nije završni znanstveni rezultat. Nemoj dodavati `--include-48000` dok je Ubuntu nestabilan i dok ne provjerimo SGC 4800 rezultat.

**Ako je nestala cijela Ubuntu/WSL sesija, a ne samo Python:** restartaj WSL iz Windows terminala i na Linux strani pogledaj `journalctl -k -b --no-pager | grep -Ei 'out of memory|oom|killed process|segfault'` ili `dmesg -T | grep -Ei 'out of memory|oom|killed process|segfault'` ako log postoji/dopušta pristup. Ti su ispisi dijagnostika OS-a, ne dokaz unaprijed da je uzrok memorija. Sačuvaj v2 JSON/log: bez njih ne smijemo unaprijed tvrditi da je stage prošao. Ne mijenjati `main`, originalne datoteke ili opaženi odd seal.

[Source-only/synthetic V2 CI `36303134361` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36303134361), uključujući stvarni poziv 4-term synthetic mirror check i provjeru rekonstruiranog checkpointa; CI ne čita lokalni WSL FITS niti dokazuje izvedbu SGC 4800R.
