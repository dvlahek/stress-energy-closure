fd=-0.00201850894575251
ap=-0.00201220917438763
am=-0.00202480872683428
high_frac=0.9422589
delta=am-ap
req=abs(delta)/abs(fd*high_frac)
assert abs(delta/fd-0.00624200971374)<1e-13
assert 0<req<.01
print("E39_PASS", delta, req, 100*req)
