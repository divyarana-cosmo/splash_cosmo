import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.lss import mass_function, peaks
from colossus.halo import mass_so, splashback
import os

def gamma_zhao09(omgm, sigma8, logmh):
    fname   =  'cosmo_%s_%s'%(omgm, sigma8)
    h           =   0.70
    nspec       =   0.95
    Omega_b_h2  =   0.02298
    z_0         =   0.0
    tcmb        =   2.726
    stat        =   1

    file = open('iparam.in', 'w')
    file.write('%s\n'%fname)
    file.write('%f %f\n'%(omgm,1-omgm))
    file.write('%d\n'%(1))
    file.write('%f\n'%h)
    file.write('%f\n'%sigma8)
    file.write('%f\n'%(0.96))
    file.write('%f\n'%(Omega_b_h2/h**2))
    file.write('%f\n'%tcmb)
    file.write('%d\n'%stat)
    file.write('%f\n'%logmh)
    file.write('%f\n'%z_0)
    file.close()


    os.system('./mandc-1.03main/mandc.x < %s'%'iparam.in')

    from glob import glob
    flist = glob('mchistory*')

    data = np.loadtxt(flist[0], skiprows=1)
    z   = data[:,0]
    Mz  = data[:,1]
    gamma  = (np.log10(Mz[5]) - np.log10(Mz[0]))/(np.log10(1+z[0]) - np.log10(1+z[5]))
    print('look here')
    os.system('rm mchistory*')
    os.system('mkdir -p ./outputs/')
    np.savetxt('outputs/res_%s'%fname,np.transpose([logmh, gamma]) )
    return gamma


def get_rsp(r200m, omgm, gamma):
    'using eq. 5 from arxiv:1504.05591'
    rsp_by_r200m = 0.54 * (1 + 0.53 * omgm) * (1 + 1.36 * np.e**(-gamma/3.04))
    return rsp_by_r200m*r200m


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
    # Convert mass to r200m and return in h-1 Mpc
    r200m   = mass_so.M_to_R(avg_mh, z=0.0, mdef='200m') / 1e3
    gamma   = splashback.modelDiemer17Gamma(nu200m, z=0.0)
    rsp     = get_rsp(r200m, omgm, gamma)

    # get the zhao09 gamma values
    gammazhao09 = gamma_zhao09(omgm, sigma8, np.log10(avg_mh))

    #rsp     = get_rsp(r200m, omgm, gammazhao09)

    return  r200m, gamma, rsp, gammazhao09

omgmarr     = np.linspace(0.15, 0.45,20)
sigma8arr   = np.linspace(0.65,1.05,20)

r200mmat        = np.zeros((len(omgmarr), len(sigma8arr)))
gammamat        = np.zeros((len(omgmarr), len(sigma8arr)))
rspmat          = np.zeros((len(omgmarr), len(sigma8arr)))
gamma_zhao09mat = np.zeros((len(omgmarr), len(sigma8arr)))

for ii,omgm in enumerate(omgmarr):
    for jj,sigma8 in enumerate(sigma8arr):
        r200mmat[ii,jj], gammamat[ii,jj], rspmat[ii,jj], gamma_zhao09mat[ii, jj] = get_r200m(omgm, sigma8)
        print(omgm, sigma8)



plt.figure(figsize=(7.5, 6))
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
im = plt.imshow(gammamat.T/gamma_zhao09mat.T, origin='lower', aspect='auto',
                extent=[omgmarr[0], omgmarr[-1], sigma8arr[0], sigma8arr[-1]])
plt.ylabel(r'$\sigma_8$', fontsize=12)
plt.xlabel(r'$\Omega_{\rm m}$', fontsize=12)
plt.colorbar(im, label=r'$\Gamma/\Gamma_{\rm zhao-09}$')


plt.subplot(3,3,4)
im = plt.imshow(rspmat.T, origin='lower', aspect='auto',
                extent=[omgmarr[0], omgmarr[-1], sigma8arr[0], sigma8arr[-1]])
plt.ylabel(r'$\sigma_8$', fontsize=12)
plt.xlabel(r'$\Omega_{\rm m}$', fontsize=12)
plt.colorbar(im, label=r'$r_{\rm sp}(\Gamma_{\rm Diemer+17})$')
#plt.colorbar(im, label=r'$r_{\rm sp}(\Gamma_{\rm zhao09})$')


plt.subplot(3,3,5)
im = plt.imshow((rspmat/r200mmat).T, origin='lower', aspect='auto',
                extent=[omgmarr[0], omgmarr[-1], sigma8arr[0], sigma8arr[-1]])
plt.ylabel(r'$\sigma_8$', fontsize=12)
plt.xlabel(r'$\Omega_{\rm m}$', fontsize=12)
plt.colorbar(im, label=r'$r_{\rm sp}/r_{\rm 200m}$')
plt.tight_layout()
plt.savefig('test.png', dpi=600)

