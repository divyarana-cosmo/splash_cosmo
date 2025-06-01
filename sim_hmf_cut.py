import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from colossus.halo import mass_so


def get_r200m(omgm, sigma8):
    # Setup the cosmology
    params = {'flat': True, 'H0': 100, 'Om0': omgm, 'Ob0': 0.049, 'sigma8': sigma8, 'ns': 0.95}
    cosmo = cosmology.setCosmology('myCosmo', **params)

    # Use log-spaced mass array for better resolution
    xmh = np.logspace(14, 18, 500)

    # Get dn/dlnM
    mfunc_so = mass_function.massFunction(xmh, 0.0, mdef='200m', model='tinker08', q_out='dndlnM')

    # Compute average mass using proper integration over ln(M)
    lnM = np.log(xmh)
    dlnM = np.gradient(lnM)
    avg_mh = np.sum(mfunc_so * xmh * dlnM) / np.sum(mfunc_so * dlnM)

    # Convert mass to r200m and return in h⁻¹ Mpc
    return mass_so.M_to_R(avg_mh, z=0.0, mdef='200m') / 1e3

omgmarr     = np.linspace(0.15, 0.45,20)
sigma8arr   = np.linspace(0.65,1.05,20)

r200mmat = np.zeros((len(omgmarr), len(sigma8arr)))

for ii,omgm in enumerate(omgmarr):
    for jj,sigma8 in enumerate(sigma8arr):
        r200mmat[ii,jj] = get_r200m(omgm, sigma8)
        print(omgm, sigma8)

plt.figure(figsize=(8, 6))
plt.subplot(2,2,1)
im = plt.imshow(r200mmat.T, origin='lower', aspect='auto',
                extent=[omgmarr[0], omgmarr[-1], sigma8arr[0], sigma8arr[-1]])
plt.ylabel(r'$\sigma_8$', fontsize=12)
plt.xlabel(r'$\Omega_m$', fontsize=12)
plt.colorbar(im, label=r'$r_{\rm 200m} [{\rm h^{-1}Mpc}]$')
plt.savefig('test.png', dpi=600)

