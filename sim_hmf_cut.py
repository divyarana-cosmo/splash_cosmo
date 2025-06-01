import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from colossus.halo import mass_so

def get_r200m(omgm, sigma8):
    #setup the cosmology
    params = {'flat': True, 'H0': 100, 'Om0': omgm, 'Ob0': 0.049, 'sigma8': sigma8, 'ns': 0.95}
    cosmo = cosmology.setCosmology('myCosmo', **params)
    #get the average halo mass
    xmh      = np.linspace(1e14,1e16,500)
    mfunc_so = mass_function.massFunction(xmh, 0.0, mdef = '200m', model = 'tinker08')
    avg_mh   = sum(mfunc_so*xmh)/sum(mfunc_so)
    #convert it to r200m
    return mass_so.M_to_R(avg_mh, z=0.0, mdef='200m')/1e3

print(get_r200m(0.3,0.8))

omgmarr     = np.linspace(0.15, 0.45,20)
sigma8arr   = np.linspace(0.65,1.05,20)

r200mmat = np.zeros((len(omgmarr), len(sigma8arr)))

for ii,omgm in enumerate(omgmarr):
    for jj,sigma8 in enumerate(sigma8arr):
        r200mmat[ii,jj] = get_r200m(omgm, sigma8)
        print(omgm, sigma8)

#plt.subplot(2,2,1)
plt.figure(figsize=(8, 6))
im = plt.imshow(r200mmat.T, origin='lower', aspect='auto',
                extent=[omgmarr[0], omgmarr[-1], sigma8arr[0], sigma8arr[-1]])
plt.ylabel(r'$\sigma_8$', fontsize=12)
plt.xlabel(r'$\Omega_m$', fontsize=12)
plt.colorbar(im, label=r'$r_{\rm 200m} [h^{-1}Mpc]$')
plt.savefig('test.png', dpi=600)

