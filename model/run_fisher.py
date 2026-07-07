import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Import GetDist
from getdist import plots
from getdist.gaussian_mixtures import GaussianND

# Import your optimized SplashCosmology engine
sys.path.append('/home/rana/github_0/splash_cosmo/model/')
from sim_peakheight import SplashCosmology

# ==========================================
# 1. Initialize Engine Globally
# ==========================================
z_l, z_s, pi_max_fid = 0.4, 0.6, 40.0
sim = SplashCosmology(Om0_fid=0.27, sigma8_fid=0.81, zl=z_l, zs=z_s, pi_max_fid=pi_max_fid)

def get_model_vector(theta, R_bins_ds, R_bins_xicg):
    """Evaluates the model and returns the concatenated 1D data vector."""
    Om0, sigma8, log_rho_s, log_alpha, log_r_s, log_rho_0, s_e, log_r_t, log_beta, log_gamma = theta

    dk14_params = {
        'Log_Rho_s': log_rho_s, 'Log_Alpha': log_alpha, 'Log_R_s': log_r_s,
        'Log_Rho_0': log_rho_0, 'S_e': s_e, 'Log_R_t': log_r_t,
        'Log_Beta': log_beta, 'Log_Gamma': log_gamma,
        'F_cen': 1.0, 'R_off': 0.1
    }

    sim.update_mcmc_step(splash_params=dk14_params, Om0_tru=Om0, sigma8_tru=sigma8)

    if np.isnan(sim.r200m):
        return np.full(len(R_bins_ds) + len(R_bins_xicg), np.nan)

    ds_fid, xicg_fid = sim.get_esd_xicg(R_bins_ds, R_bins_xicg)
    return np.concatenate([ds_fid/1e12, xicg_fid]).flatten()


def compute_fisher_matrix(theta_fiducial, step_sizes, icov, R_bins_ds, R_bins_xicg):
    """
    Computes the Fisher Matrix using a highly stable 5-point stencil
    finite difference method to prevent truncation and floating point errors.
    """
    N_params = len(theta_fiducial)
    mu_fid = get_model_vector(theta_fiducial, R_bins_ds, R_bins_xicg)
    N_data = len(mu_fid)

    Jacobian = np.zeros((N_data, N_params))

    print("Calculating robust numerical derivatives (5-point stencil)...")
    for i in range(N_params):
        h = step_sizes[i]

        # Create parameter copies for the 4 evaluation points
        theta_p2 = np.copy(theta_fiducial)
        theta_p1 = np.copy(theta_fiducial)
        theta_m1 = np.copy(theta_fiducial)
        theta_m2 = np.copy(theta_fiducial)

        # Apply steps (+2h, +1h, -1h, -2h)
        theta_p2[i] += 2 * h
        theta_p1[i] += h
        theta_m1[i] -= h
        theta_m2[i] -= 2 * h

        # Evaluate model at all 4 points
        mu_p2 = get_model_vector(theta_p2, R_bins_ds, R_bins_xicg)
        mu_p1 = get_model_vector(theta_p1, R_bins_ds, R_bins_xicg)
        mu_m1 = get_model_vector(theta_m1, R_bins_ds, R_bins_xicg)
        mu_m2 = get_model_vector(theta_m2, R_bins_ds, R_bins_xicg)

        # 5-point stencil derivative formula
        # f'(x) = (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / 12h
        Jacobian[:, i] = (-mu_p2 + 8*mu_p1 - 8*mu_m1 + mu_m2) / (12.0 * h)

        print(f"  Derivative for parameter {i+1}/{N_params} completed.")

    # F = J^T * icov * J
    Fisher = Jacobian.T @ icov @ Jacobian
    return Fisher


if __name__ == "__main__":
    # ---------------------------------------------------------
    # A. Data Loading & Setup
    # ---------------------------------------------------------
    R_bins_ds, ds, ds_err = np.loadtxt('./mock_data/gen_esd.dat', unpack=1)
    R_bins_xicg, xicg, xicg_err = np.loadtxt('./mock_data/gen_xicg2d.dat', unpack=1)

    # Note: Using your updated error downweighting factor (25 instead of 100)
    cov = np.diag(np.concatenate([ds_err**2/25, xicg_err**2/25]))
    icov = np.linalg.inv(cov)

    xtrue = np.array([0.27, 0.81, 2.107, -0.462, -0.203, 0.210, 0.987, 0.193, 0.841, 0.213])

    # Dynamic fractional step sizes (2% for cosmo, 0.5% for shape)
    fractional_steps = np.array([0.02, 0.02, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005])
    step_sizes = np.maximum(1e-5, fractional_steps * np.abs(xtrue))

    param_names = ["Om0", "sigma8", "log_rho_s", "log_alpha", "log_r_s", "log_rho_0", "s_e", "log_r_t", "log_beta", "log_gamma"]
    param_labels = [r"$\Omega_m$", r"$\sigma_8$", r"$\log\rho_s$", r"$\log\alpha$", r"$\log r_s$", r"$\log\rho_0$", r"$s_e$", r"$\log r_t$", r"$\log\beta$", r"$\log\gamma$"]

    # ---------------------------------------------------------
    # B. Compute Base Data Fisher Matrix
    # ---------------------------------------------------------
    F_data = compute_fisher_matrix(xtrue, step_sizes, icov, R_bins_ds, R_bins_xicg)

    # ---------------------------------------------------------
    # C. Apply DK14 Gaussian Priors (From MCMC)
    # ---------------------------------------------------------
    F_dk14_priors = np.zeros_like(F_data)

    # From lnprior: gauss(log_alpha, np.log10(0.2), 0.6) --> Index 3
    #               gauss(log_beta, np.log10(6.0), 0.2)  --> Index 8
    #               gauss(log_gamma, np.log10(4.0), 0.2) --> Index 9
    F_dk14_priors[3, 3] = 1.0 / (0.6**2)
    F_dk14_priors[8, 8] = 1.0 / (0.2**2)
    F_dk14_priors[9, 9] = 1.0 / (0.2**2)

    F_data_with_shape = F_data + F_dk14_priors

    # ---------------------------------------------------------
    # D. Apply Correlated Planck 2018 Priors
    # ---------------------------------------------------------
    sigma_Om0_planck = 0.0073
    sigma_sigma8_planck = 0.0060
    rho_planck = -0.75

    cov_prior_planck = np.array([
        [sigma_Om0_planck**2, rho_planck * sigma_Om0_planck * sigma_sigma8_planck],
        [rho_planck * sigma_Om0_planck * sigma_sigma8_planck, sigma_sigma8_planck**2]
    ])

    F_planck = np.zeros_like(F_data)
    F_planck[0:2, 0:2] = np.linalg.inv(cov_prior_planck)

    F_total = F_data_with_shape + F_planck

    # ---------------------------------------------------------
    # E. Evaluate Constraints
    # ---------------------------------------------------------
    try:
        # 1. Base Data Only
        cov_data = np.linalg.inv(F_data)
        sigma_data = np.sqrt(np.diag(cov_data))

        # 2. Data + DK14 Priors
        cov_shape = np.linalg.inv(F_data_with_shape)
        sigma_shape = np.sqrt(np.diag(cov_shape))

        # 3. Data + DK14 Priors + Planck
        cov_tot = np.linalg.inv(F_total)
        sigma_tot = np.sqrt(np.diag(cov_tot))

        # Note: These printed 1-sigma errors are derived directly from the un-truncated
        # matrices. If a Gaussian error extends beyond a prior bound, the true marginalized
        # error within the prior volume would technically be smaller.
        print("\n--- Fisher Forecast Constraints (1-sigma, Pre-Truncation) ---")
        print(f"{'Parameter':>12} | {'True Val':>8} | {'Data Only':>12} | {'+ DK14 Priors':>15} | {'+ Planck':>12}")
        print("-" * 70)
        for name, val, err_d, err_s, err_t in zip(param_names, xtrue, sigma_data, sigma_shape, sigma_tot):
            print(f"{name:>12} | {val:8.4f} | +/- {err_d:8.4f} | +/- {err_s:11.4f} | +/- {err_t:8.4f}")

    except np.linalg.LinAlgError:
        print("CRITICAL: Fisher Matrix is singular. Check step sizes or model gradients.")
        sys.exit()

    os.makedirs('./fits/', exist_ok=True)

    # ---------------------------------------------------------
    # F. GetDist Visualization with Hard Bounds Applied
    # ---------------------------------------------------------
    print("\nRendering GetDist Triangle Plot within prior volume...")
    from getdist import plots
    from getdist.gaussian_mixtures import GaussianND

    # Extract hard boundary uniform priors from mpi_mcmc_cosmo.py
    # This prevents the GetDist contour plotting from extending into non-physical space
    prior_bounds = {
        'Om0': [0.1, 0.4],
        'sigma8': [0.6, 1.5],
        'log_rho_s': [-3.0, 5.0],
        'log_r_s': [np.log10(0.1), np.log10(5.0)],
        'log_rho_0': [-1.5, 1.5],
        's_e': [0.1, 4.0],
        'log_r_t': [np.log10(0.5), 0.4]
    }

    # Initialize the Gaussians, explicitly passing the ranges argument
    fisher_base = GaussianND(xtrue, cov_data, names=param_names, labels=param_labels, label='Data Only', ranges=prior_bounds)
    fisher_shape = GaussianND(xtrue, cov_shape, names=param_names, labels=param_labels, label='+ DK14 Priors', ranges=prior_bounds)
    fisher_all = GaussianND(xtrue, cov_tot, names=param_names, labels=param_labels, label='+ Planck', ranges=prior_bounds)

    g = plots.get_subplot_plotter()
    g.settings.figure_legend_frame = False
    g.settings.alpha_filled_add = 0.5

    # Plot to see how the contours shrink at each step
    g.triangle_plot(
        [fisher_base, fisher_shape, fisher_all],
        filled=True,
        contour_colors=['gray', 'blue', 'darkred'],
        markers={name: val for name, val in zip(param_names, xtrue)}
    )

    out_file = './fits/fisher_corner_with_all_priors.png'
    g.export(out_file)
    print(f"Saved prior-overlay corner plot to {out_file}")
