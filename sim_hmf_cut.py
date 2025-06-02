import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.lss import mass_function, peaks
from colossus.halo import mass_so, splashback


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
    print(np.log10(avg_mh))
    nu200m = peaks.peakHeight(M=avg_mh, z=0.0)
    # Convert mass to r200m and return in h⁻¹ Mpc
    r200m = mass_so.M_to_R(avg_mh, z=0.0, mdef='200m') / 1e3
    gamma = splashback.modelDiemer20Gamma(nu200m, z=0.0)
    rsp   = splashback.modelDiemer20RspR200m(gamma, nu200m, z=0.0, rspdef= 'sp-apr-p50' )*r200m
    return  r200m, gamma, rsp

omgmarr     = np.linspace(0.15, 0.45,20)
sigma8arr   = np.linspace(0.65,1.05,20)

r200mmat = np.zeros((len(omgmarr), len(sigma8arr)))
gammamat = np.zeros((len(omgmarr), len(sigma8arr)))
rspmat = np.zeros((len(omgmarr), len(sigma8arr)))

for ii,omgm in enumerate(omgmarr):
    for jj,sigma8 in enumerate(sigma8arr):
        r200mmat[ii,jj], gammamat[ii,jj], rspmat[ii,jj] = get_r200m(omgm, sigma8)
        print(omgm, sigma8)

plt.figure(figsize=(8, 6))
plt.subplot(3,3,1)
im = plt.imshow(r200mmat.T, origin='lower', aspect='auto',
                extent=[omgmarr[0], omgmarr[-1], sigma8arr[0], sigma8arr[-1]])
plt.ylabel(r'$\sigma_8$', fontsize=12)
plt.xlabel(r'$\Omega_{\rm m}$', fontsize=12)
plt.colorbar(im, label=r'$r_{\rm 200m}$')
#plt.colorbar(im, label=r'$r_{\rm 200m} [{\rm h^{-1}Mpc}]$')

plt.subplot(3,3,2)
im = plt.imshow(gammamat.T, origin='lower', aspect='auto',
                extent=[omgmarr[0], omgmarr[-1], sigma8arr[0], sigma8arr[-1]])
plt.ylabel(r'$\sigma_8$', fontsize=12)
plt.xlabel(r'$\Omega_{\rm m}$', fontsize=12)
plt.colorbar(im, label=r'$\Gamma$')

plt.subplot(3,3,3)
im = plt.imshow(rspmat.T, origin='lower', aspect='auto',
                extent=[omgmarr[0], omgmarr[-1], sigma8arr[0], sigma8arr[-1]])
plt.ylabel(r'$\sigma_8$', fontsize=12)
plt.xlabel(r'$\Omega_{\rm m}$', fontsize=12)
plt.colorbar(im, label=r'$r_{\rm sp}$')


#plt.subplot(3,3,4)
#im = plt.imshow((rspmat/r200mmat).T, origin='lower', aspect='auto',
#                extent=[omgmarr[0], omgmarr[-1], sigma8arr[0], sigma8arr[-1]])
#plt.ylabel(r'$\sigma_8$', fontsize=12)
#plt.xlabel(r'$\Omega_{\rm m}$', fontsize=12)
#plt.colorbar(im, label=r'$r_{\rm sp}/r_{\rm 200m}$')
plt.tight_layout()
plt.savefig('test.png', dpi=600)

