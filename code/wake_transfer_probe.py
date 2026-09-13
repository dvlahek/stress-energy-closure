#!/usr/bin/env python3
from pathlib import Path
import json, tempfile
import numpy as np
from classy import Class
import class_response_optimize as cro

q=np.linspace(0,20,4000)
f0, *_ = cro.kinetic_objects(q,0.06,1100.0)
path=Path('probe_fd.dat')
np.savetxt(path,np.column_stack([q,f0]),fmt='%.14e')
p={
 'output':'mPk,dTk,vTk', 'modes':'s', 'gauge':'newtonian',
 'H0':cro.H0,'omega_b':cro.OMEGA_B,'omega_cdm':cro.OMEGA_CDM,
 'A_s':cro.A_S,'n_s':0.9649,'tau_reio':cro.TAU_REIO,
 'N_ur':cro.N_UR,'N_ncdm':1,'use_ncdm_psd_files':1,
 'ncdm_psd_filenames':str(path.resolve()),'m_ncdm':0.06,'T_ncdm':cro.T_NCDM,'deg_ncdm':1.0,
 'P_k_max_h/Mpc':0.3,'z_max_pk':1.0,'z_pk':'0,0.1,0.5'
}
c=Class(); c.set(p); c.compute()
out={}
for z in [0.0,0.1,0.5]:
 t=c.get_transfer(z=z,output_format='class')
 out[str(z)]={'keys':list(t.keys()),'shapes':{k:list(np.shape(v)) for k,v in t.items()}}
 for k in t:
  arr=np.asarray(t[k])
  if arr.ndim==1 and arr.size>5:
   out[str(z)].setdefault('samples',{})[k]=[float(arr[1]),float(arr[min(10,arr.size-1)]),float(arr[-2])]
print(json.dumps(out,indent=2))
Path('wake_transfer_probe.json').write_text(json.dumps(out,indent=2)+'\n')
c.struct_cleanup(); c.empty()
