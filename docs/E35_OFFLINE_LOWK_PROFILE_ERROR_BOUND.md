# E35 — stroga niskokvalna ograda profilnoga faktora, ali ne i UV certifikat

**29. 9. 2026. · Analitički rezultat za pozitivne sferne profile i offline toy QA; nije izračun stvarne E28 sile.** E32 je pokazao da Abacusovi L2 percentili ne određuju potpuno `u(k)`. E34 je dokazao da dva čak i idealna endpoint profila ne određuju retardirani integral. Ovdje razdvajamo dodatno pitanje: pod kojim se *neovisnim ulaznim ograničenjem* na kratke prostorne skale može dobiti konzervativna kontrola profilnoga faktora?

## 1. Teorem

Neka je `dM(r)/M` nenegativna, masom normirana sferna radijalna raspodjela s konačnim drugim momentom `<r²>`. Njezin sferni Fourierov faktor iznosi `u(k)=∫sinc(kr)dM/M`. Iz identiteta `sinc(x)=1/2 ∫_{−1}^1 cos(x μ)dμ` te nejednakosti `0≤1−cos y≤y²/2` egzaktno slijedi

\[
0\le 1-u(k)\le \frac{k^2\langle r^2\rangle}{6},
\qquad |u(k)|\le1.
\]

Za pozitivni sferni source i test profil, bez pretpostavke NFW-a, slijedi

\[
\boxed{|1-u_{\rm src}(k,t)u_{\rm test}(k,t_{\rm obs})|
\le \min\!\left(2,\frac{k^2}{6}
[\langle r^2\rangle_{\rm src}(t)
+\langle r^2\rangle_{\rm test}(t_{\rm obs})]\right).}
\]

Za fiksni vanjski linearni Born kernel `K(k,t)`, pozitivni maseni omjer `m(t)` i iste source/test normalizacije, usporedba s **točkastim profilnim faktorom** `u_src*u_test=1` daje

\[
|A_{\rm profile}(k)-A_{\rm point}(k)|
\le \frac{k^2}{6}\int dt\,|K(k,t)|m(t)
[\langle r^2\rangle_{\rm src}(t)+
 \langle r^2\rangle_{\rm test}(t_{\rm obs})].
\]

Ovo je lokalna ograda po `k`, a ne dokaz male pogreške integrala preko cijeloga E28 pojasa. Stvarni gravitacijski kernel može promijeniti predznak i uključivati dodatne geometrijske čimbenike. Apsolutni iznos kernela u gornjoj ogradi ne dopušta da se oslanjamo na slučajno poništavanje potpisanih integranada.

## 2. Kontrole i kontra-primjer endpoint ograde

Na dvjema pozitivnim konačnim sfernim smjesama toy QA daje `<r²>src=0.638`, `<r²>test=0.355`. Pri bezdimenzijskom `k=1` stvarna apsolutna pogreška profilnoga proizvoda iznosi `0.153739548347`, a stroga momentna ograda `0.1655`. Na većim bezdimenzijskim k momentna ograda postaje labava; trivijalna gornja ograda iz `|u|≤1` jest `2`. Svih 85 provjera, uključujući vremenske integrale za obje E34 povijesti, prošlo je.

Kontra-primjer pokazuje zašto se endpoint momenti ne smiju koristiti kao cjelovremenska ograda. Ljuska s identičnim `r=1` na oba kraja, ali prijelaznim `r=10` u sredini, pri bezdimenzijskom `k=.15` daje stvarnu apsolutnu pogrešku produkta `0.337494275778`. **Pogrešna ograda iz endpoint radijusa** bila bi `0.0075` i očito je prekršena. Ograda iz *stvarnoga trenutnog* radijusa ostaje valjana. Ovo je kinematički, a ne samogravitacijski model haloa.

## 3. Što rezultat omogućuje, a što ne

Ako budući fizički model neovisno ograniči source i test `<r²>` na **svakom** relevantnom vremenu, možemo unaprijed označiti k područje u kojem je zamjena profilnog faktora točkastim dokazivo kontrolirana. Izbor tolerancije mora proizaći iz fizičkog error budgeta i istoga observablea, a ne iz F± kontrasta ili opažene eBOSS vrijednosti. Ovaj matematički rezultat **ne** određuje fizički `kmax`, ne zatvara anisotropiju, ne daje cijeli L1 host ili `M200c`, ne uključuje environment i ne ograničava `f2` nelinearnoga Vlasovljeva odziva.

Zbog E28 rezultata u kojem je 94.22589% uvjetne FD sile nižega sidra i alpha=.4 iz `k>1 Mpc^−1`, niskokvalni bound sam po sebi ne certificira ukupni drag. Opisne igračke dimenzijske vrijednosti iz E35 nisu u fizičkim E28 k ili R jedinicama.

Izvorni `e35_offline_lowk_profile_bound_result.json` i `e35_offline_lowk_profile_bound_qa.py` dokumentiraju roditeljske SHA256 i 85/85 PASS. `A03 PHYSICAL_UNCERTIFIED`, `A04 BLOCKED`, observed odd **SEALED**. Bez WSL-a, ASDF/FITS, novih podataka, prilagodbi F± ili GitHub izmjena.
