import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Import GetDist
from getdist import plots, MCSamples

# Ensure paths match your MCMC implementation
sys.path.append('/home/rana/github_0/splash_cosmo/')
from sim_peakheight import SplashCosmology
from model_dk14 import splash

# ==========================================
# 1. Initialize Engine Globally
# ==========================================
z_l, z_s, pi_max_fid = 0.4, 0.6, 40.0
sim = SplashCosmology(Om0_fid=0.27, sigma8_fid=0.81, zl=z_l, zs=z_s, pi_max_fid=pi_max_fid)

def get_model_vector(theta, R_bins_ds, R_bins_xicg):
    Om0, sigma8, log_rho_s, log_alpha, log_r_s, log_rho_0, s_e, log_r_t, log_beta, log_gamma = theta
    
    ss = splash(
        Log_Rho_s=log_rho_s, Log_Alpha=log_alpha, Log_R_s=log_r_s,
        Log_Rho_0=log_rho_0, S_e=s_e, Log_R_t=log_r_t,
        Log_Beta=log_beta, Log_Gamma=log_gamma,
        F_cen=1.0, R_off=0.1, R_max=pi_max_fid, splrmin=0.1, splrbin=60
    )

    sim.update_mcmc_step(splash_obj=ss, Om0_tru=Om0, sigma8_tru=sigma8)
    
    if np.isnan(sim.r200m):
        return np.full(len(R_bins_ds) + len(R_bins_xicg), np.nan)

    ds_fid, xicg_fid = sim.get_esd_xicg(R_bins_ds, R_bins_xicg)
    return np.concatenate([ds_fid/1e12, xicg_fid]).flatten()


def compute_fisher_matrix(theta_fiducial, step_sizes, icov, R_bins_ds, R_bins_xicg):
    N_params = len(theta_fiducial)
    mu_fid = get_model_vector(theta_fiducial, R_bins_ds, R_bins_xicg)
    N_data = len(mu_fid)
    
    Jacobian = np.zeros((N_data, N_params))
    
    print("Calculating robust numerical derivatives (5-point stencil)...")
    for i in range(N_params):
        h = step_sizes[i]
        
        theta_p2, theta_p1 = np.copy(theta_fiducial), np.copy(theta_fiducial)
        theta_m1, theta_m2 = np.copy(theta_fiducial), np.copy(theta_fiducial)
        
        theta_p2[i] += 2 * h
        theta_p1[i] += h
        theta_m1[i] -= h
        theta_m2[i] -= 2 * h
        
        mu_p2 = get_model_vector(theta_p2, R_bins_ds, R_bins_xicg)
        mu_p1 = get_model_vector(theta_p1, R_bins_ds, R_bins_xicg)
        mu_m1 = get_model_vector(theta_m1, R_bins_ds, R_bins_xicg)
        mu_m2 = get_model_vector(theta_m2, R_bins_ds, R_bins_xicg)
        
        Jacobian[:, i] = (-mu_p2 + 8*mu_p1 - 8*mu_m1 + mu_m2) / (12.0 * h)
        print(f"  Derivative for parameter {i+1}/{N_params} completed.")

    return Jacobian.T @ icov @ Jacobian


def filter_samples(raw_samples, bounds):
    """Filters Monte Carlo samples strictly against uniform bounds."""
    mask = np.ones(len(raw_samples), dtype=bool)
    for i, (lower, upper) in enumerate(bounds):
        if lower is not None:
            mask &= (raw_samples[:, i] >= lower)
        if upper is not None:
            mask &= (raw_samples[:, i] <= upper)
    return raw_samples[mask]


if __name__ == "__main__":
    R_bins_ds, ds, ds_err = np.loadtxt('./mock_data/gen_esd.dat', unpack=1)
    R_bins_xicg, xicg, xicg_err = np.loadtxt('./mock_data/gen_xicg2d.dat', unpack=1)

    cov = np.diag(np.concatenate([ds_err**2/25, xicg_err**2/25]))
    icov = np.linalg.inv(cov)

    xtrue = np.array([0.27, 0.81, 2.10782687, -0.46212836, -0.20334147, 
                      0.21067979, 0.98738642, 0.1939794, 0.84172758, 0.21323246])
    
    fractional_steps = np.array([0.02, 0.02, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005])
    step_sizes = np.maximum(1e-5, fractional_steps * np.abs(xtrue))
    
    param_names = ["Om0", "sigma8", "log_rho_s", "log_alpha", "log_r_s", "log_rho_0", "s_e", "log_r_t", "log_beta", "log_gamma"]
    param_labels = [r"$\Omega_m$", r"$\sigma_8$", r"$\log\rho_s$", r"$\log\alpha$", r"$\log r_s$", r"$\log\rho_0$", r"$s_e$", r"$\log r_t$", r"$\log\beta$", r"$\log\gamma$"]

    # Exact strict uniform bounds from mpi_mcmc_cosmo.py
    param_bounds = [
        (0.1, 0.4),                           # Om0 
        (0.6, 1.5),                           # sigma8 
        (-3.0, 5.0),                          # log_rho_s
        (None, None),                         # log_alpha
        (np.log10(0.1), np.log10(5.0)),       # log_r_s
        (-1.5, 1.5),                          # log_rho_0
        (0.1, 4.0),                           # s_e
        (np.log10(0.5), 0.4),                 # log_r_t
        (None, None),                         # log_beta
        (None, None)                          # log_gamma
    ]

    # Compute Base Matrix
    F_data = compute_fisher_matrix(xtrue, step_sizes, icov, R_bins_ds, R_bins_xicg)

    # 1. Apply DK14 Gaussian Priors
    F_dk14_priors = np.zeros_like(F_data)
    F_dk14_priors[3, 3] = 1.0 / (0.6**2)  
    F_dk14_priors[8, 8] = 1.0 / (0.2**2)  
    F_dk14_priors[9, 9] = 1.0 / (0.2**2)  
    
    F_dk14_only = F_data + F_dk14_priors

    # 2. Add Correlated HSC 3x2pt Prior on Om0 and sigma8
    # Using representative constraints for the HSC S8 degeneracy
    sigma_Om0_hsc = 0.040
    sigma_sigma8_hsc = 0.045
    rho_hsc = -0.75 

    cov_prior_hsc = np.array([
        [sigma_Om0_hsc**2, rho_hsc * sigma_Om0_hsc * sigma_sigma8_hsc],
        [rho_hsc * sigma_Om0_hsc * sigma_sigma8_hsc, sigma_sigma8_hsc**2]
    ])
    
    F_hsc = np.zeros_like(F_data)
    F_hsc[0:2, 0:2] = np.linalg.inv(cov_prior_hsc) # Inject the inverted 2x2 block

    F_dk14_and_hsc = F_dk14_only + F_hsc

    try:
        cov_dk14 = np.linalg.inv(F_dk14_only)
        cov_hsc = np.linalg.inv(F_dk14_and_hsc)
        
        print("\nGenerating Monte Carlo samples to apply strict truncation...")
        np.random.seed(42)
        n_samples = 3000000
        
        # Sample & filter DK14-only
        raw_dk14 = np.random.multivariate_normal(xtrue, cov_dk14, size=n_samples)
        valid_dk14 = filter_samples(raw_dk14, param_bounds)
        
        # Sample & filter DK14 + HSC 3x2pt
        raw_hsc = np.random.multivariate_normal(xtrue, cov_hsc, size=n_samples)
        valid_hsc = filter_samples(raw_hsc, param_bounds)
        
        print(f"DK14 Only: Retained {len(valid_dk14)} samples out of {n_samples}.")
        print(f"DK14 + HSC: Retained {len(valid_hsc)} samples out of {n_samples}.")

        # Create GetDist objects
        getdist_ranges = {name: bnd for name, bnd in zip(param_names, param_bounds) if bnd != (None, None)}
        
        samples_dk14 = MCSamples(samples=valid_dk14, names=param_names, labels=param_labels, label='Data + DK14', ranges=getdist_ranges)
        samples_hsc = MCSamples(samples=valid_hsc, names=param_names, labels=param_labels, label='Data + DK14 + HSC 3x2pt', ranges=getdist_ranges)
        
        err_dk14 = np.sqrt(samples_dk14.getVars())
        err_hsc = np.sqrt(samples_hsc.getVars())

        print("\n--- True Synced Constraints (1-sigma, within Prior Volume) ---")
        print(f"{'Parameter':>12} | {'True Val':>8} | {'DK14 Only':>15} | {'+ HSC 3x2pt':>15}")
        print("-" * 61)
        for name, val, e1, e2 in zip(param_names, xtrue, err_dk14, err_hsc):
            print(f"{name:>12} | {val:8.4f} | +/- {e1:11.4f} | +/- {e2:11.4f}")
            
    except np.linalg.LinAlgError:
        print("CRITICAL: Fisher Matrix is singular. Check step sizes or model gradients.")
        sys.exit()

    os.makedirs('./fits/', exist_ok=True)

    # Visualization
    print("\nRendering GetDist Triangle Plot...")
    
    g = plots.get_subplot_plotter(width_inch=14)
    g.settings.axes_fontsize = 12
    g.settings.axes_labelsize = 16
    g.settings.lab_fontsize = 20
    g.settings.legend_fontsize = 20
    g.settings.x_label_rotation = 45 
    g.settings.figure_legend_frame = False
    g.settings.alpha_filled_add = 0.5

    # Plot both sets of contours to see the constraint difference
    g.triangle_plot(
        [samples_dk14], 
        filled=True, 
        contour_colors=['blue'],
        markers={name: val for name, val in zip(param_names, xtrue)}
    )

    out_file = './fits/fisher_corner_with_hsc.pdf'
    g.export(out_file)
    print(f"Saved comparative corner plot to {out_file}")
