# eBOSS DR16: puni mock-galaxy uzorak umjesto 600D, sparse 4800→48000R — jedan lokalni WSL run

**27. 9. 2026. · Tehnički cilj:** upotrijebiti sve već autentificirane i dostupne galaksije originalnog realistic EZmock ID `0001` za NGC i SGC, uz isti stari estimator i potpuno iste source-matched tracer-specific randome. Ovo **nije** opservacijski eBOSS odd unblinding, fizički certifikat prozora, nezavisna mock kovarijanca, Einstein–Vlasov predikcija ni detekcija.

## Zašto se mijenja veličina uzorka, ne science cut

Iz [izvornog 18-case E4 JSON-a](../source_data/eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_report_2026-09-27.json), koji ostaje nepromijenjen, originalnih `600D` i `1200→2400→4800R` daju medijan NGC `D1D2=1718` odnosno **11,93 prihvaćenih DD parova po 144 s–μ ćelije**; SGC `D1D2=5248` odnosno **36,44 po ćeliji**. DD je isti na sve tri R razine. Pri `4800R` medijan NGC `R1R2=114758` i SGC `323198`, ali nova RR gustoća **ne mijenja nedovoljno velik originalni DD poduzorak**. E4 relativni pomak pune `6×24 ξ` mreže pri 2400→4800 ima medijan NGC 0,438741 i SGC 0,434716. Ne predstavljati RR podršku 144/144 kao konvergenciju niti iz ovog pomaka izdvajati DD/DR/RD/RR doprinos po ćeliji: povijesni E4 izvještaj čuva SHA histograma, ne njegove 144 pojedinačne vrijednosti.

Originalni A-02 izvorno potvrđuje da mock0001 nakon *istog već odabranog* `z∈[0.9,1.0)` i pozitivnog `WEIGHT_SYSTOT×WEIGHT_CP×WEIGHT_NOZ×WEIGHT_FKP` ima NGC `7537 LRG-D/17149 ELG-D` i SGC `4159 LRG-D/15860 ELG-D`. To nije nova selekcija. Puni mock izvori i randome već čuva WSL. [Zaseban prospektivni *engineering* protokol nakon E4](../source_data/eboss_dr16_full_eligible_mock_galaxy_4800_48000_random_pilot_protocol_2026-09-27.json) fiksira jedan originalni ID, dvije kape, izvornih 4800R i dodatnih 43 200 istog-source random redaka po traceru/kapi; **48 000R nije post-E4 fizički prihvatni prag**, nego deset puta gušći numerički dijagnostički pilot, bez novih preuzimanja.

## Kopiraj cijeli blok u WSL

Izvršavanje i *stderr i stdout* idu istodobno na terminal i u trajni log. `set -o pipefail` čuva stvarni Python exit status i pri `tee`. Nema prepisivanja starog JSON-a: ako originalni output već postoji, program odbija ponovno pokretanje.

```bash
set -Eeuo pipefail
cd "$HOME/stress-energy-closure"
BRANCH="audit/eboss-elg-bit8-ra-orientation-20260925"
git fetch origin "$BRANCH"
git switch "$BRANCH"
git pull --ff-only origin "$BRANCH"
source .venv/bin/activate
mkdir -p eboss_workspace/a03_full_eligible_mock
PYTHONUNBUFFERED=1 python -u scripts/run_eboss_dr16_full_eligible_mock0001_sparse.py --self-test
PYTHONUNBUFFERED=1 python -u scripts/run_eboss_dr16_full_eligible_mock0001_sparse.py --run \
  2>&1 | tee eboss_workspace/a03_full_eligible_mock/full_eligible_id0001_4800_48000_stdout.log
printf 'FULL_MOCK_RUN_EXIT_OK\n'
```

**Predaj točno dva lokalna izlaza:** `eboss_workspace/a03_full_eligible_mock/full_eligible_id0001_4800_48000_report.json` i `eboss_workspace/a03_full_eligible_mock/full_eligible_id0001_4800_48000_stdout.log`. Uspješan run završava porukom `EBOSS_FULL_ELIGIBLE_MOCK0001 ... DESCRIPTIVE_ONLY`, sa `0001/NGC` i `0001/SGC` te `OBSERVED_ODD_DATA_READ False`. Ako se program zaustavi, **ne ponavljati** nakon izmjene inputa ili brisanja failure JSON-a; predaj izvorni JSON/log za neovisni audit.

## Što izračunava i što provjerava

[Runner](../scripts/run_eboss_dr16_full_eligible_mock0001_sparse.py) prvo čita samo arhivirane JSON/protokole i provjerava izvornu implementaciju. Zatim u prethodnoj A-02 funkciji *prije prvog novog FITS headera/retka* ponovno potpuno SHA256 verificira svih **72** izvorišnih gzip datoteka (zbroj **1 860 198 719** bajtova). Odbija novi download ili alternativni izvor. Za sva četiri izvora svake kape ponovno reproducira originalni A-02 `600D/1200R` sample SHA i ulaznu eligible populaciju. Za random izvore dodatno reproducira originalni E4 `4800R` index+catalogue SHA prije determinističkog izvlačenja dodatnih 43 200 redaka iz isključivo izvornog eligible komplementa.

Prije velikih parova novi sparse algoritam mora reproducirati **svih osam** originalnih fwd/reverse weighted `600D/1200R` histogramskih SHA iz A-02 na originalnom dense algoritmu, a sparse verzija mora dati identičan prihvaćeni broj parova i numerički podudarne histograme. Puni D×D, D×R, R×D i R×R zatim se računaju na istim šest s-binova × 24 signed-μ binova sa *zasebnom normalizacijom sva četiri LS člana*, za oba smjera LRG→ELG i ELG→LRG, za 4800 i 48 000 randoma. Svi ponderirani histogrami `6×24`, svi `ξ` i `ell=0,1,2,3` šest-bin projekcije zapisuju se u izlaz zajedno sa SHA i full-D izvorima. Originalni 4800R RR prihvaćeni parovi i ponderirani ukupni zbroj moraju ponoviti E4. Bez umjetnog zero-filla: multipoli se daju samo uz punu RR potporu.

[Sintetički/source-only GitHub Actions QA `36301575145` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36301575145) već potvrđuje source guard i sintetičko slaganje stare dense i nove sparse geometrije. Taj CI **ne čita lokalni WSL FITS** i ne može potvrditi stvarni NGC/SGC rezultat prije izvršavanja gornje naredbe.

**Ograničenje:** čak i ako full-D `4800→48000R` daje manji numerički pomak, rezultat je i dalje jedan mock ID s obje korelirane kape, ne fizički validiran cijeli survey prozor ili statistička značajnost. Buduća eBOSS kovarijanca zahtijeva dovoljno neovisnih matched realizacija i unaprijed definiran konačni observable; apsolutna Einstein–Vlasov amplituda zahtijeva zasebnu fizičku kalibraciju. Opaženi odd ostaje SEALED.
