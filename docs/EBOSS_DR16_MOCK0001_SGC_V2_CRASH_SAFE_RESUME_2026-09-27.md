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

## Stvarni stage 4800R PASS, zatim sigurno ograničeni nastavak 48000R (V2.1)

Korisnikov dostavljeni SGC stage 4800R **stvarno je dovršen**: report 193920 bajtova, SHA256 `d3f298a9550981c0ece7d6c9e7aee0d3fb5cebaf0d80cdd39805cda3bb08e238`, stdout 9692 bajta, SHA256 `bad556ba41190aa8414255e3e82f81b679a7e08f7d4c690d454128461a3014a9`; originalni v1 NGC canonical SHA nepromijenjen. NGC 4800+48000 ostaje samo izvorni parent. SGC 4800 checkpoint sadrži sva četiri izvorno forward-brojana člana: DD 1006690, DR 304094, RD 1064361, RR 323373, uz puni RR support 144/144 i originalni E4 4800R RR count+sum replay. [Kompaktni audit korisničkog izvještaja](../source_data/eboss_dr16_mock0001_sgc_v2_uploaded_4800_stage_source_only_audit_2026-09-27.json) reproducira sve 24 kombinacije forward/reverse weighted histogram SHA u tri već izračunata cap×R slučaja, točno tri LS `ξ` i 12 exact-`μ` six-bin multipole vektora. To je **source-only interna audit provjera dostavljenog izvještaja**, ne neovisno ponovno čitanje originalnog lokalnog FITS-a. Izvorni veliki SGC independent-reverse fail i dalje nije dijagnosticiran; V2 SGC reverse je algebarski.

[Pre-48000R inženjerski V2.1 protokol](../source_data/eboss_dr16_mock0001_sgc_48000_chunked_v21_engineering_protocol_2026-09-27.json) uvodi **dodatne checkpointove unutar svakog od triju 48000R parnih članova**: do 4000 uzastopnih redaka prvog kataloga po podbloku, s cijelim drugim katalogom, bez izmjene izvornog para, redshift prozora, survey geometrije, random redaka, weights, binova i fizikalne predikcije. Suma disjunktnih redaka u izvornom redoslijedu koristi normalizaciju računanu nad **cijelim izvornim prvim i drugim katalogom**, ne zasebno normalizirane podblokove. Razlikovati moguće zadnje bitove zbrajanja floating point histograma od fizikalne promjene estimatora. Svaki podblok dobiva zaseban SHA i atomarni zapis; restart ga ponovno provjeri i preskoči. SGC 4800R i originalni NGC se ne računaju ponovno. Izvorno brojani reverse veliki SGC parovi se ne računaju, niti se v1 anomaly prikriva.

**Samo nakon ovog dostavljenog i provjerenog 4800R stagea** pokreni drugi, eksplicitni nastavak:

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
python -u scripts/resume_eboss_dr16_full_eligible_mock0001_sgc_v2.py --self-test
python -u scripts/resume_eboss_dr16_full_eligible_mock0001_sgc_v2.py --run --include-48000 \
  2>&1 | tee -a eboss_workspace/a03_full_eligible_mock/full_eligible_id0001_4800_48000_resume_v2_stdout.log
```

Ako se WSL ponovno ugasi, **ne brisati ni originalni v1 ni v2 output**. Pokretanjem potpuno iste `--run --include-48000` naredbe provjeravaju se isti izvori i nastavlja prvi nedovršeni **4000-row podblok**, bez ponavljanja dovršenih SGC 48000R parova/podblokova. Nakon završetka spremi novi v2 JSON i isti append-only log; očekuje se `V2_SGC_COMPLETE_FROZEN_NGC_REUSED`, ali ga ne prijavljivati bez stvarnog outputa. [Sintetička source-only V2.1 CI `36303866242` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36303866242) provjerava 4 disjunktna podbloka, originalni nesjeckani parni kernel, ponovno učitavanje SHA-checkpointa i negativnu tamper kontrolu. CI **nema** korisnikove originalne FITS datoteke; SGC 48000R numerički rezultat još je PENDING.

**Ako je nestala cijela Ubuntu/WSL sesija, a ne samo Python:** restartaj WSL iz Windows terminala i na Linux strani pogledaj `journalctl -k -b --no-pager | grep -Ei 'out of memory|oom|killed process|segfault'` ili `dmesg -T | grep -Ei 'out of memory|oom|killed process|segfault'` ako log postoji/dopušta pristup. Ti su ispisi dijagnostika OS-a, ne dokaz unaprijed da je uzrok memorija. Sačuvaj v2 JSON/log: bez njih ne smijemo unaprijed tvrditi da je stage prošao. Ne mijenjati `main`, originalne datoteke ili opaženi odd seal.

[Source-only/synthetic V2 CI `36303134361` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36303134361), uključujući stvarni poziv 4-term synthetic mirror check i provjeru rekonstruiranog checkpointa; CI ne čita lokalni WSL FITS niti dokazuje izvedbu SGC 4800R.
