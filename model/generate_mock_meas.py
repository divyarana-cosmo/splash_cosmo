import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Import your finalized engine
from sim_peakheight import SplashCosmology

# ---------------------------------------------------------
# 1. Load the Data and Covariance Matrices
# ---------------------------------------------------------
print("Loading data and covariance matrices...")
rbins_xicg, best_xicg = np.loadtxt('/net/dobbe/data2/github/splash_cosmo/DataStore/joshi_2026/xi_cg_best.dat', unpack=1)
cov_xicg = np.loadtxt('/net/dobbe/data2/github/splash_cosmo/DataStore/joshi_2026/RASS_paper_data/splashback/catalog/xi_2d_cov_lz.dat')

rbins_esd, meas_esd = np.loadtxt('/net/dobbe/data2/github/splash_cosmo/DataStore/joshi_2026/RASS_paper_data/WL_data/dsig_mean.dat', unpack=1)
idx = rbins_esd > 0.2
rbins_esd = rbins_esd[idx]
cov_esd = np.loadtxt('/net/dobbe/data2/github/splash_cosmo/DataStore/joshi_2026/RASS_paper_data/WL_data/cov.dat')
cov_esd = np.delete(cov_esd, ~idx, axis=0)
cov_esd = np.delete(cov_esd, ~idx, axis=1)

# Assign to target variables
R_bins_ds = rbins_esd
R_bins_xicg = rbins_xicg


# ---------------------------------------------------------
# 2. Initialize the SplashCosmology Engine
# ---------------------------------------------------------
print("Initializing SplashCosmology engine...")
z_l = 0.4
z_s = 0.6
pi_max_fid = 40.0
Om0_fid = 0.27
sigma8_fid = 0.81

# Initialize purely with the fiducial cosmology state
sim = SplashCosmology(
    Om0_fid=Om0_fid, sigma8_fid=sigma8_fid, 
    zl=z_l, zs=z_s, pi_max_fid=pi_max_fid
)

# ---------------------------------------------------------
# 3. Define the True Mock Parameters
# ---------------------------------------------------------
# Mock True Cosmology
Om0_tru = 0.27
sigma8_tru = 0.81

# Mock DK14 Shape Parameters
log_rho_s, log_alpha, log_r_s, log_rho_0, s_e, log_r_t, log_beta, log_gamma = [
    2.10782687, -0.46212836, -0.20334147, 0.21067979, 
    0.98738642, 0.1939794, 0.84172758, 0.21323246
]

dk14_params = {
    'Log_Rho_s': log_rho_s, 'Log_Alpha': log_alpha, 'Log_R_s': log_r_s,
    'Log_Rho_0': log_rho_0, 'S_e': s_e, 'Log_R_t': log_r_t,
    'Log_Beta': log_beta, 'Log_Gamma': log_gamma,
    'F_cen': 1.0, 'R_off': 0.1  # Ensuring off-centering params are provided
}

# ---------------------------------------------------------
# 4. Generate the Mock Observables
# ---------------------------------------------------------
print("\n--- Processing ---")
# Apply the parameters directly to the engine
sim.update_mcmc_step(splash_params=dk14_params, Om0_tru=Om0_tru, sigma8_tru=sigma8_tru)

print(f"Dynamically inferred 3D target_rsp: {sim.target_rsp:.4f} Mpc/h")

# Verify the physical parameters mapped correctly
if np.isnan(sim.r200m):
    print("CRITICAL ERROR: Root finder failed to map the target Rsp to a physical mass.")
    sys.exit()
else:
    print(f"Mapped Halo Mass (M200m): {np.log10(sim.m200m):.2f} log(Msun/h)")

# Generate the theoretical data vectors mapped to the Fiducial frame
ds_fid, xicg_fid = sim.get_esd_xicg(R_bins_ds=R_bins_ds, R_bins_xicg=R_bins_xicg)

# ---------------------------------------------------------
# 5. Save the Mock Data
# ---------------------------------------------------------
os.makedirs('./mock_data', exist_ok=True)

np.savetxt(
    './mock_data/gen_esd.dat', 
    np.transpose([R_bins_ds, ds_fid/1e12, np.diag(cov_esd)**0.5]), 
    header='rbins esd esderr'
)

np.savetxt(
    './mock_data/gen_xicg2d.dat', 
    np.transpose([R_bins_xicg, xicg_fid, np.diag(cov_xicg)**0.5]), 
    header='rbins xicg2d xicg2derr'
)
print("Saved mock data to ./mock_data/")

# ---------------------------------------------------------
# 6. Plotting
# ---------------------------------------------------------
fig = plt.figure(figsize=(12, 5))

# Delta Sigma Plot
ax1 = fig.add_subplot(1, 2, 1)
ax1.errorbar(R_bins_ds, ds_fid/1e12, yerr=np.diag(cov_esd)**0.5, fmt='o', capsize=3, color='blue')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel(r'$R_{fid} \, [h^{-1}\mathrm{Mpc}]$')
ax1.set_ylabel(r'$\Delta\Sigma \, [h M_\odot / \mathrm{pc}^2]$')
ax1.grid(True, ls="--", alpha=0.5)

# Xi_cg Plot
ax2 = fig.add_subplot(1, 2, 2)
ax2.errorbar(R_bins_xicg, xicg_fid, yerr=np.diag(cov_xicg)**0.5, fmt='s', capsize=3, color='purple')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel(r'$R_{fid} \, [h^{-1}\mathrm{Mpc}]$')
ax2.set_ylabel(r'$\xi_{cg}(R)$')
ax2.grid(True, ls="--", alpha=0.5)

plt.tight_layout()
plt.savefig("mock_meas.png", dpi=300)
print("Saved plot as 'mock_meas.png'")
