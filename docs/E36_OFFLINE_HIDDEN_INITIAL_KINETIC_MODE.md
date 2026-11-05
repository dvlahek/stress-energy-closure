# E36 — iste početne gustoće, brzine i tlak ne određuju kasniji slobodni Vlasovljev odziv

**29. 9. 2026. · Konstruktivni 1D free-transport rezultat s pozitivnim distribucijama.** Ovo je izdvojeni matematički svjedok neodređenoga početnog kinetičkog stanja. Nije puni Einstein–Vlasov izvod, nije samogravitacijski Vlasov–Poisson halo, nova E27/E28 numerika niti opaženi eBOSS rezultat.

## 1. Problem: halo povijest nije jedini nedostajući izvor

E34 je pokazao da čak ni dva potpuno poznata profila haloa ne određuju izvor između snimaka. Međutim, za lineariziranu Vlasovljevu evoluciju potreban je i početni faznoprostorni neutrinski poremećaj. Općeniti Duhamelov izraz ima oblik

\[
\delta f(t)=U_0(t,t_i)\,\delta f(t_i)
+\int_{t_i}^{t}dt'\,U_0(t,t')\,L_{\Psi(t')}F.
\]

Prvi član predstavlja slobodnu evoluciju prethodno postojećega kinetičkog stanja, a drugi odziv na specificirani vanjski gravitacijski izvor. E27/E28 postavljaju prvi član na nulu pri z=1 kao uvjetni benchmark. Nije izvedeno da stvarni prethodni wake mora biti nula.

## 2. Pozitivna konstrukcija s tri jednaka početna momenta

Uzmimo `x` na periodičnoj kružnici, `v∈R`, standardnu normiranu Gaussovu pozadinu `F0(v)` i funkciju

\[
h(v)=\cos v-\frac{e^{3/2}}4\cos(2v)
-\frac34e^{-1/2}.
\]

Definiramo dva pozitivna početna polja pri `eps=0.1`,

\[
f_\pm(x,v,0)=F_0(v)[1\pm\varepsilon\cos x\,h(v)].
\]

Globalna pozitivnost slijedi iz `eps |h(v)| ≤ eps[1+e^(3/2)/4+3e^(-1/2)/4] < 1`. Iz egzaktnih Gaussovih karakterističnih funkcija slijedi `∫F0 h dv=0`, `∫v F0 h dv=0` i `∫v² F0 h dv=0`. Zato oba polja u **svakoj** prostornoj točki imaju istu početnu gustoću `rho=1`, momentum `j=0` i drugi brzinski moment `P=1`, iako je puna brzinska raspodjela različita. Ovi se jednaki Newtonovi momenti ne smiju zamijeniti tvrdnjom da su sve relativističke komponente početnog `T_mu_nu` identične.

Pod **istom vanjskom nultom silom**, jednadžba `∂_t f +v∂_x f=0` ima egzaktno rješenje `f_±(x,v,t)=f_±(x-vt,v,0)`. Njegove gustoće iznose

\[
\rho_\pm(x,t)=1\pm\varepsilon\cos x\,G(t),\qquad
G(t)=-\frac12e^{-(t^2+1)/2}[\cosh(t)-1]^2.
\]

Iako je `G(0)=0` i prva tri vremenska izvoda gustoćnog moda također nestaju u nuli, `G(t)<0` za `t≠0` te je `G(t)=-(e^{-1/2}/8)t^4+O(t^6)`. Pri **isključivo bezdimenzijskom** t=1 i x=0 račun daje `G(1)=-0.054250551363639` i `rho_+-rho_-=-0.0108501102727278`. Ne radi se o mjerljivom neutrinskom ili galaktičkom kontrastu; amplituda eps odabrana je samo radi pozitivnosti i demonstracije.

## 3. Stroga granica interpretacije

Slobodna struja stvara gustoćni poremećaj koji bi zatim u samogravitacijskom Vlasov–Poisson sustavu promijenio potencijal. Zato prethodno egzaktno rješenje *slobodnoga transporta* nije egzaktna puna samogravitacijska evolucija nakon početka. U Einstein–Vlasov sustavu treba dodatno zahtijevati kompletnu početnu metriku, ograničenja, sve komponente stress-energy tenzora i vlastitu dinamiku. Ovaj svjedok ne certficira dvije pune fizički dopuštene Einstein–Vlasov kozmologije s istim svim početnim izvorima.

**Zaključak relevantan E27–E29:** čak ni početni halo podaci i tri najniža Newtonova neutrinska momenta ne određuju homogeni kinetički član linearnog forward modela. Da bismo uvjetni E27/E28 wake pretvorili u fizički ukupan odziv, raniji `delta f(t_i)` treba specificirati na temelju stvarnih fizikalnih početnih podataka ili neovisno opravdano ograničiti njegov doprinos.

Skripta `e36_offline_hidden_kinetic_initial_data_qa.py` s izvornim `e36_offline_hidden_kinetic_initial_data_result.json` ima 179/179 PASS, uključujući test analitičkog `G`, numeričku Gaussovu integraciju i negativne kontrole. Nema ASDF, FITS, WSL, GitHub promjene, novog F± modela ili observed eBOSS odd pristupa. `A03 PHYSICAL_UNCERTIFIED`, `A04 BLOCKED`, opaženi odd **SEALED**.
