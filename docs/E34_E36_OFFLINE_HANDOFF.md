# EinsteinVlasovNP — E34–E36 offline handoff (29. 9. 2026.)

## Zaključena nova matematika, bez novih fizičkih podataka

**E34, retarded source-history witness (51/51 QA PASS):** dvije pozitivne ljuske imaju istu masu, središte, radijus, profil i nultu radijalnu brzinu u oba granična snimka; njihov kontinuitetno konzistentan međuvremenski radijus se razlikuje. Uz istu *toy*, pozitivnu kauzalnu težinu dobivamo različite retardirane integrale. Dva potpuna halo snimka ne identificiraju općenitu međuvremensku povijest. Nije test stvarnoga E28 Kernela i nije samogravitirajući halo. Za ogradu interpolacijske pogreške nužan je izvana opravdan `L=sup|du/dt|`.

**E35, low-k spherical moment theorem (85/85 PASS):** `|1-u_src u_test|<=k²(<r²>src+<r²>test)/6` za pozitivne sferno usrednjene normirane profile. Unutar istog eksternog linearnog kernela slijedi integrirana apsolutna point-source replacement ograda pomoću vremenski promjenjivih drugih radijalnih momenata. Same endpoint momenti ne daju tu ogradu; toy prolazno proširenje s identičnim krajevima služi kao egzaktna negativna kontrola. Ovo nije UV ili F± fizički certifikat i ne određuje `kmax`.

**E36, independent initial kinetic state witness (179/179 PASS):** eksplicitni, svuda pozitivni početni `f_±` imaju točno istu gustoću, momentum i 1D brzinski drugi moment u svakoj prostornoj točki, ali različite više momente. Pod istim slobodnim Vlasovljevim transportom imaju kasniji `rho_±=1±eps cos x G(t)` s egzaktno `G(t)=-0.5 exp[-(t²+1)/2](cosh t-1)^2`. Dakle čak tri najniža Newtonova početna momenta ne određuju homogeni kinetički član buduće gustoće. Slobodni transport nije samogravitirajući/Einstein–Vlasov kompletan sustav; ne tvrditi isti full GR `T_mu_nu`, početnu metriku ili Einstein constraints.

## Integracija s cijelim projektom

- Zamrznuti E27/E28R1 rezultati ostaju **uvjetni**. Za niže maseno sidro, alpha=.4, 94.22589% FD finite-band ubrzanja dolazi iz k>1/Mpc; ne postoji dokazana ukupna fizička UV konvergencija niti valjana eBOSS galaxy amplitude. E25 ukupna neutrinska gravitacija smije se brojiti samo jednom, ne kao dodatni proizvoljni Euler residual.
- Stvarnih 20 E29 objekata ostaju **lokalno PID-podržani u testiranim prozorima**, ne globalni merger-tree certifikat ni reprezentativan katalog. L1 `N` nije `M200c`. L2 radial percentiles ne definiraju puni `u(k)` (E32); promjena SO mase može uključiti pseudo-evoluciju (E33). E34 dodatno dokazuje da čak i dva potpuna krajnja profila ne daju kontinuirani source; E36 pokazuje nezavisnu neodređenost početnog neutrinskoga DF-a.
- E30 precizira uvjetnu ravnotežu vjetra i razliku između `E[T_w|sign]` te linearnoga `P_delta,v`. Neuvjetovani 24D eBOSS EV predložak i high-z `chi_F=A_L c_E-A_E c_L` ostaju fizikalno neidentificirani. Devet eBOSS mockova ne može dati punorangovnu 24D uzoračku kovarijancu.

## Sljedeći stvarni fizički input, ne još jedan metadata gate

1. Dobaviti/izvesti dinamički i kozmološki dosljednu halo+environment povijest kroz cijeli retarded interval, uz nezavisno ograničene source/test second moments ako se želi koristiti E35 low-k bound. Dva halo snimka i spline ne zadovoljavaju E34.
2. Specificirati početni neutrinski DF iz prethodne kozmološke evolucije ili dati provjerljiv bound za homogeni član Vlasovljeva rješenja. Nulti raniji wake iz E27 je benchmark, ne deducirano opažanje; E36 daje konstrukciju zašto.
3. Neovisno validirati kratkovalni f2/linear Born breakdown i common halo-environment force. Tek tada pokušati fizički `beta_a,epsilon_a` i high-z LRG/ELG `chi_F`/window; ne invertiraj 24D kovarijancu iz devet mockova.

## Nezaobilazna ograničenja i reprodukcija

GitHub repo `dvlahek/stress-energy-closure`, samo branch `audit/eboss-elg-bit8-ra-orientation-20260925` i draft PR #1. Ovaj E34–36 paket NIJE mijenjao GitHub niti `main`, izvorni F0/F+/F−, E8 4000q, E16, E7 48k, E27/E28, E4, cuts, seeding ili LOS. Opaženi eBOSS galaxy rows i 24D odd ostaju **SEALED**. `A03 PHYSICAL_UNCERTIFIED`, `A04 BLOCKED`.

Paket sadrži reproducibilne standard-library Python 3 skripte i izvorne JSON-ove, protokole, SHA256 manifest i male parent dokumente. Nakon raspakiravanja: `python replay_e34_e36_offline.py`. Kod prvo verificira sve parent/output SHA256, zatim vrti QA, pa ponovno provjerava bajtni identitet. Nisu potrebni WSL, internet, ASDF/FITS, class ili posebni paketi.
