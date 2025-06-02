import os
import numpy as np
import matplotlib.pyplot as plt

lmhalo  = np.linspace(13,15,5)
gamma = np.array([])

if 0:
    for omg in [0.2,0.3,0.4]:
        for sigma8 in [0.7, 0.8, 0.9]:
            gamma = np.array([])
            for mm in lmhalo:
                fname   =  'cosmo_%s_%s_%s'%(omg, sigma8, mm)
                Omega_0     =   0.35
                h           =   0.70
                nspec       =   0.95
                Omega_b_h2  =   0.02298
                z_0         =   0.0
                tcmb        =   2.726
                stat        =   1

                file = open('iparam.in', 'w')
                file.write('%s\n'%fname)
                file.write('%f %f\n'%(omg,1-omg))
                file.write('%d\n'%(1))
                file.write('%f\n'%h)
                file.write('%f\n'%sigma8)
                file.write('%f\n'%(0.96))
                file.write('%f\n'%(Omega_b_h2/h**2))
                file.write('%f\n'%tcmb)
                file.write('%d\n'%stat)
                file.write('%f\n'%mm)
                file.write('%f\n'%z_0)
                file.close()


                os.system('./mandc.x < %s'%'iparam.in')

                from glob import glob
                flist = glob('mchistory*')

                data = np.loadtxt(flist[0], skiprows=1)
                z   = data[:,0]
                Mz  = data[:,1]
                xx  = (np.log10(Mz[5]) - np.log10(Mz[0]))/(np.log10(1+z[0]) - np.log10(1+z[5]))

                #h0 = 100/(3.086*1e+19)/3.17098e-8
                #xx = dM_dt/M_0/h0
                gamma = np.append(gamma,xx)
                print('look here')
                os.system('rm %s'%flist[0])
                print(mm, xx)
            np.savetxt('outputs/res_%s'%fname,np.transpose([lmhalo, gamma]) )

from glob import glob
omg_flist = np.sort(glob('outputs/res_cosmo_0.3_*'))
sig8_flist = np.sort(glob('outputs/res_cosmo_*_0.8_*'))
ax = plt.subplot(2,2,1)
axsigma8 = plt.subplot(2,2,2)

for fil in omg_flist:
    lmhalo, gamma = np.loadtxt(fil, unpack=1)
    omg     = float(fil.split('_')[2])
    sigma8  = float(fil.split('_')[3])
    ax.plot(lmhalo, gamma, label=r'$\Omega_{\rm m} = %2.2f, \sigma_8 =%2.2f$'%(omg,sigma8))
    ax.set_xlabel(r'$\log M_{\rm h}$')
    ax.set_ylabel(r'$\Gamma$')


ax.set_ylim(2.5,5.5)
ax.legend()

for fil in sig8_flist:
    lmhalo, gamma = np.loadtxt(fil, unpack=1)
    omg     = float(fil.split('_')[2])
    sigma8  = float(fil.split('_')[3])
    axsigma8.plot(lmhalo, gamma, label=r'$\Omega_{\rm m} = %2.2f, \sigma_8 =%2.2f$'%(omg,sigma8))
    axsigma8.set_xlabel(r'$\log M_{\rm h}$')
    axsigma8.set_ylabel(r'$\Gamma$')

axsigma8.set_ylim(2.5,5.5)
axsigma8.legend()

plt.tight_layout()



plt.savefig('outputs/test.png', dpi=300)


#omega_ms    = [0.35, 0.25, 0.25]
#sigma8s     = [0.8, 0.8, 0.7]
#
#params  =   zip(omega_ms, sigma8s)
#
#plt.subplot(2,2,1)
#for nn, (omg, sig) in enumerate(params):
#    print(nn, omg, sig)
#    data = np.loadtxt('PWGH_average.dat_%d'%nn)
#
#    idx = data[:,1]<4.0
#    plt.plot(np.log10(1+data[idx,1]), data[idx,3],label='$\Omega_m = %s, \sigma_8 = %s$'%(omg, sig))
#
#
#plt.xlabel(r'$\log[1+z]$')
#plt.ylabel(r'$\log\langle M(z)/M_0 \rangle$')
##plt.xscale('log')
##plt.ylim(-2,0)
#plt.legend()
##plt.yscale('log')
#plt.savefig('test.png', dpi=300)


