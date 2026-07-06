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

# Import your optimized SplashCosmology engine
sys.path.append('/home/rana/github_0/splash_cosmo/')
from sim_peakheight import SplashCosmology

# ===================================================================
# GLOBAL ENGINE INITIALIZATION
# Initialize once globally so MPI workers inherit the fiducial cache 
# in memory without recalculating it on every single MCMC step.
# ===================================================================
z_l = 0.4
z_s = 0.6
pi_max_fid = 40.0
sim = SplashCosmology(Om0_fid=0.27, sigma8_fid=0.81, zl=z_l, zs=z_s, pi_max_fid=pi_max_fid)


def gauss(x, mean, sig):
    val = np.exp(-(x-mean)**2/(2*sig**2)) / np.sqrt(2*np.pi*sig**2)
    return val

def model(x, R_bins_ds, R_bins_xicg):
    """
    Evaluates the theoretical model by updating the global engine in-place.
    """
    Om0_tru, sigma8_tru, log_rho_s, log_alpha, log_r_s, log_rho_0, s_e, log_r_t, log_beta, log_gamma = x
    
    # Pack parameters (Include default off-centering if not sampled)
    dk14_params = {
        'Log_Rho_s': log_rho_s, 'Log_Alpha': log_alpha, 'Log_R_s': log_r_s,
        'Log_Rho_0': log_rho_0, 'S_e': s_e, 'Log_R_t': log_r_t,
        'Log_Beta': log_beta, 'Log_Gamma': log_gamma,
        'F_cen': 1.0, 'R_off': 0.1  
    }

    # Highly optimized in-place update for the true cosmology and shape
    sim.update_mcmc_step(splash_params=dk14_params, Om0_tru=Om0_tru, sigma8_tru=sigma8_tru)
    
    # If parameters push the profile into a non-physical space, 
    # the root finder will return NaN for r200m. Catch it gracefully.
    if np.isnan(sim.r200m):
        return None

    # Retrieve mapped observables
    ds_fid, xicg_fid = sim.get_esd_xicg(R_bins_ds=R_bins_ds, R_bins_xicg=R_bins_xicg)
    
    return np.concatenate([ds_fid/1e12, xicg_fid]).flatten()


def lnprior(x):
    """Ryoma's 2020 Priors"""
    Om0, sigma8, log_rho_s, log_alpha, log_r_s, log_rho_0, s_e, log_r_t, log_beta, log_gamma = x
    
    if (0.1 < Om0 < 0.4 and 0.6 < sigma8 < 1.5 and -3 <= log_rho_s <= 5 and 
        np.log10(0.1) <= log_r_s <= np.log10(5.0) and -1.5 <= log_rho_0 <= 1.5 and 
        0.1 <= s_e <= 4 and np.log10(0.5) <= log_r_t <= 0.2):
        
        val = np.log(gauss(log_alpha, np.log10(0.2), 0.6) * gauss(log_beta, np.log10(6.0), 0.2) * gauss(log_gamma, np.log10(4.0), 0.2))
        return val
    return -np.inf


def lnprob(x, data, icov, R_bins_ds, R_bins_xicg):
    lp = lnprior(x)
    if not np.isfinite(lp):
       dirt = 5*np.ones(len(data) + 1)
       return -np.inf, dirt
       
    begin = time.time() 
    mod = model(x, R_bins_ds, R_bins_xicg)
    
    # Catch both NaN arrays and the 'None' trigger from root finder failure
    if mod is None or np.sum(~np.isfinite(mod)) > 0:
        dirt = 5*np.ones(len(data) + 1)
        return -np.inf, dirt

    Delta = mod - data
    chisq = np.dot(Delta, np.dot(icov, Delta))
    
    blob = np.append(mod, chisq)
    res = lp - chisq * 0.5	
	
    return res, blob


def runchain(Ntotal, sampler, chainf, blobf, pos, nwalkers):
    blnk = [""]
    fchain = open(chainf, "w")
    fblob = open(blobf, "w")
    iterno = 1
    
    for result in sampler.sample(pos, iterations=Ntotal, store=1):
        posn, probn, staten, blobsn = result
        for i in range(nwalkers):
            np.savetxt(fchain, posn[i], newline=' ')
            np.savetxt(fchain, [sampler.acceptance_fraction[i], -2.*probn[i]], newline=' ')
            np.savetxt(fblob, blobsn[i], newline=' ')
            np.savetxt(fchain, blnk, fmt='%s')
            np.savetxt(fblob, blnk, fmt='%s')
        
        if iterno % 100 == 0:
            print("Iteration number: %d of %d done" % (iterno, Ntotal))
        iterno += 1
        posnew = result[0]

    fchain.close()
    fblob.close()
    return posnew


if __name__ == "__main__":
    # Load mock measurements
    R_bins_ds, ds, ds_err = np.loadtxt('./mock_data/gen_esd.dat', unpack=1)
    R_bins_xicg, xicg, xicg_err = np.loadtxt('./mock_data/gen_xicg2d.dat', unpack=1)

    data = np.concatenate([ds, xicg])
    # Downweighting error for mock test stability
    cov  = np.diag(np.concatenate([ds_err**2/25, xicg_err**2/25]))
    icov = np.linalg.inv(cov)

    # Add in a MPI pool of workers
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    
    ndim = 10
    nwalkers = 64
    np.random.seed(1996)
    
    xtrue = np.array([0.27, 0.81, 2.10782687, -0.46212836, -0.20334147, 0.21067979, 0.98738642, 0.1939794, 0.84172758, 0.21323246])
    
    # Initialize Walkers with a 5% Gaussian scatter around the true values
    # np.random.randn(ndim) generates standard normal variables. We multiply by 0.05 * abs(xtrue) to scale the std deviation.
    p_0 = [xtrue + 0.05 * np.abs(xtrue) * np.random.randn(ndim) for _ in range(nwalkers)]

    # Initialize the sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, args=[data, icov, R_bins_ds, R_bins_xicg])
   
    dirname = './fits/' 
    os.makedirs(dirname, exist_ok=True)

    print("Starting iteration, ndim is ", ndim)
    Nburn = 6000
    print("Running Burn-in Phase...")
    pos = runchain(Nburn, sampler, f"{dirname}/burnfile_cosmo.dat", f"{dirname}/burnpredfile_cosmo.dat", p_0, nwalkers)
 
    sampler.reset()

    Ntotal = 6000
    print("Running Main Production Chain...")
    posfinal = runchain(Ntotal, sampler, f"{dirname}/chainfile_cosmo.dat", f"{dirname}/predfile_cosmo.dat", pos, nwalkers)

    pool.close()
