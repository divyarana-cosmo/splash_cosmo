import numpy as np
import matplotlib.pyplot as plt
from glob import glob

flist = glob('mandc_*/mandcoutput.*')

plt.subplot(2,2,1)
for fil in flist:
    ispec,nbin1,zobs,omegam0,omegaL0,sigma8,tilt,omegab0,omeganu0,T_cmb,N_nu,h =  np.loadtxt(fil, max_rows=1, unpack=True,dtype=str)
    print(fil,sigma8)
    ziz,Miz,c,rhohiz,R,V,Ms,rhos,Rs,Vs,Miz_200c,c_200c,crhohiz_200c,Miz_200m,c_200m,rhohiz_200m,uniage_iz =  np.loadtxt(fil, skiprows=1, unpack=True)


    plt.plot(1+ziz, Miz_200m/1e13, label='$\Omega_m = %s, \sigma_8 = %s$'%(omegam0[:4], sigma8[:4]))

plt.xlim(1.0,4.0)
plt.legend(fontsize='small')
plt.xlabel('$1+z$')
plt.ylabel(r'$M_{\rm 200m}(z)/10^{13}$')

plt.xscale('log')
#plt.yscale('log')
plt.savefig('test.png', dpi=300)
