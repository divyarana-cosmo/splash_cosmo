#!/usr/bin/env python3
import numpy as np
import pandas as pd
import emcee
import sys
import time
from schwimmbad import MPIPool
import os
import matplotlib.pyplot as plt
from scipy import integrate
sys.path.append('/home/rana/github_0/splash_cosmo/')
from model_dk14 import splash

def gauss(x,mean,sig):
    val = np.exp(-(x-mean)**2/(2*sig**2))/np.sqrt(2*np.pi*sig**2)
    return val

def model_dk14(x, ):
    log_rho_s, log_alpha, log_r_s, log_rho_0, s_e, log_r_t, log_beta, log_gamma = x
    ss = splash(Log_Rho_s = log_rho_s, Log_Alpha = log_alpha, Log_R_s = log_r_s, Log_Rho_0 = log_rho_0, S_e = s_e, Log_R_t = log_r_t, Log_Beta = log_beta, Log_Gamma = log_gamma, F_cen = 1.0, R_off = 0.3, splrmin=0.1, splrbin=60)
 
    xi_2d = 10**ss.log_xi_2d(logr)
    xi_3d = 10**ss.log_xi_3d(logr)
    diff_2d = ss.diff_xi_2d(logr)
    diff_3d = ss.diff_xi_3d(logr)
    rsp_2d = ss.rsp_2d(logr)
    rsp_3d = ss.rsp_3d(logr)
    return np.concatenate([xi_2d, xi_3d, diff_2d, diff_3d, rsp_2d, rsp_3d]).flatten()
    #return xi_2d, xi_3d, diff_2d, diff_3d, rsp_2d, rsp_3d

#Here the priors are followed from ryoma's paper
def lnprior(x):
    log_rho_s, log_alpha, log_r_s, log_rho_0, s_e, log_r_t, log_beta, log_gamma = x
    if -3<=log_rho_s<=5 and np.log10(0.1) <= log_r_s <= np.log10(5.0) and -1.5 <= log_rho_0 <= 1.5 and 0.1 <= s_e <= 4 and np.log10(0.5)<=log_r_t<=0.2 :
        val = 0.0 + np.log(gauss(log_alpha,np.log10(0.2),0.6) * gauss(log_beta,np.log10(6.0),0.2) * gauss(log_gamma,np.log10(4.0),0.2))
        return val
    return -np.inf

def lnprob(x, data, icov, log_rbin):
    lp = lnprior(x)
    if not np.isfinite(lp):
       dirt = 5*np.ones(4*len(data) + 5)
       return -np.inf,dirt
    
    mod = model_dk14(x,log_rbin)
    #if np.sum(~np.isfinite(mod))>0:
    #   dirt = 5*np.ones(len(data))
    #   return -np.inf,dirt

    Delta = mod[:int(len(log_rbin))] - data
    chisq = np.dot(Delta, np.dot(icov, Delta))
    blob = mod
    blob = np.append(blob,chisq)
    print(np.shape(blob))

    print('log_rho_s, log_alpha, log_r_s, log_rho_0, s_e, log_r_t, log_beta, log_gamma, chisq')
    print(x,chisq)
    print('blob is xi2d, xi3d, diff2d, diff3d, rsp2d, rsp3d, chisq')
    print('size of blob %d'%(len(blob)))
    res = lp-chisq*0.5	
	
    return res,blob



def runchain(Ntotal,sampler,chainf,blobf,pos, nwalkers):
    blnk=[""];
    fchain=open(chainf,"w");
    fblob=open(blobf,"w");
    iterno=1;
    # Store chainfile and prednfile in the same format as before
    for result in sampler.sample(pos, iterations=Ntotal, store=1):
        posn,probn,staten,blobsn = result;
        #posn,probn,staten = result;
        for i in range(nwalkers):
            np.savetxt(fchain,posn[i],newline=' ');
            np.savetxt(fchain,[sampler.acceptance_fraction[i],-2.*probn[i]],newline=' ');
            np.savetxt(fblob,blobsn[i],newline=' ');
            np.savetxt(fchain,blnk,fmt='%s');
            np.savetxt(fblob,blnk,fmt='%s');
        print("Iteration number: %d of %d done"%(iterno,Ntotal));
        iterno=iterno+1;
        posnew=result[0];

    fchain.close();
    fblob.close();
    return posnew;


if __name__ == "__main__":
    #read the xicg  and esd
    rbins_esd, esd, esderr      = np.loadtxt('./mock_data/gen_esd.dat', unpack=1)
    rbins_xicg, xicg, xicgerr   = np.loadtxt('./mock_data/gen_xicg2d.dat', unpack=1)

    icovesd     = np.linalg.inv(np.diag(esderr**2/100))
    icovxicg    = np.linalg.inv(np.diag(xicgerr**2/100))

    # Add in a MPI pool of workers
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        print('something is shitty')
        sys.exit(0)
    
    #log_rho_s, log_alpha, log_r_s, log_rho_0, s_e, log_r_t, log_beta, log_gamma = x
    x0  = [ 1.45638209, -0.38738921, -0.26249276, -0.21003311,  1.3623482,   0.02528298, 1.05831751,  0.83872683]
    ndim = 8
    nwalkers = 512
    np.random.seed(1996)
    
    scat = np.random.uniform(-0.05,0.05,nwalkers)

    p_log_rho_s = np.random.uniform(0.8, 3.2, nwalkers) #
    p_log_alpha = np.random.uniform(-2, 0.0, nwalkers)  #
    p_log_r_s   = np.random.uniform(-1, 0.4, nwalkers)  #
    p_log_rho_0 = np.random.uniform(0.1, 0.8, nwalkers) #
    p_s_e       = np.random.uniform(0.1, 3.5, nwalkers) #
    p_log_r_t   = np.random.uniform(-0.1, 0.2, nwalkers)#)
    p_log_beta  = np.random.uniform(0.6, 0.9, nwalkers) #
    p_log_gamma = np.random.uniform(0.3, 0.9, nwalkers) #


    p_0 = list(zip(p_log_rho_s, p_log_alpha, p_log_r_s, p_log_rho_0, p_s_e, p_log_r_t, p_log_beta, p_log_gamma))

    # Initialize the sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, args=[data, icov, log_rbin])
   
    dirname = './fits/' 
    os.system("mkdir -p %s" % (dirname))

    print("Starting iteration, ndim is ",ndim);
    Nburn = 3000
    pos = runchain(Nburn,sampler,"%s/burnfile.dat"%(dirname), "%s/burnpredfile.dat"%(dirname), p_0, nwalkers);
 
    sampler.reset();

    Ntotal=3000
    posfinal = runchain(Ntotal,sampler,"%s/chainfile.dat"%(dirname),"%s/predfile.dat"%(dirname),pos, nwalkers);

    pool.close()

