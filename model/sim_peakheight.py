import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.lss import peaks
from colossus.halo import mass_so, splashback
from scipy.optimize import root_scalar
from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.interpolate import InterpolatedUnivariateSpline as ius

class SplashCosmology:
    """
    Decoupled Cosmological Mapping Engine.
    Maps a separately initialized DK14 Splashback model from true cosmology to fiducial.
    """
    def __init__(self, Om0_fid=0.27, sigma8_fid=0.81, zl=0.4, zs=0.6, pi_max_fid=40,
                 h=1.0, nspec=0.95, Ob0=0.049, tcmb=2.726):

        self.zl = zl
        self.zs = zs
        self.pi_max_fid = pi_max_fid
        self.h = h
        self.nspec = nspec
        self.Ob0 = Ob0
        self.tcmb = tcmb

        # Setup and Permanently Cache the Fiducial Cosmology
        self.Om0_fid = Om0_fid
        self.sigma8_fid = sigma8_fid
        params_fid = {
            'flat': True, 'H0': 100.0 * self.h, 'Om0': self.Om0_fid,
            'Ob0': self.Ob0, 'sigma8': self.sigma8_fid, 'ns': self.nspec, 'Tcmb0': self.tcmb
        }
        self.cosmo_fid = cosmology.setCosmology('fid_cosmo', **params_fid, persistence='')

        # Cache invariant fiducial distances
        self.Dc_l_fid = self.cosmo_fid.comovingDistance(z_max=self.zl)
        self.Dc_s_fid = self.cosmo_fid.comovingDistance(z_max=self.zs)
        self.H_fid = self.cosmo_fid.Hz(self.zl)
        self.sigma_crit_fid = self.Dc_s_fid / (self.Dc_l_fid * (self.Dc_s_fid - self.Dc_l_fid))

        # Log grid used strictly for finding the target 3D Rsp derivative
        self.log_R_grid = np.log10(np.logspace(-1, 1, 100))


    def update_mcmc_step(self, splash_obj, Om0_tru, sigma8_tru):
        """
        Accepts a freshly instantiated splash object and scales it
        to the proposed true cosmology.
        """
        # Store the passed-in DK14 object
        self.splash = splash_obj

        # Update True Cosmology
        self.Om0_tru = Om0_tru
        self.sigma8_tru = sigma8_tru
        params_tru = {
            'flat': True, 'H0': 100.0 * self.h, 'Om0': self.Om0_tru,
            'Ob0': self.Ob0, 'sigma8': self.sigma8_tru, 'ns': self.nspec, 'Tcmb0': self.tcmb
        }
        self.cosmo_tru = cosmology.setCosmology('true_cosmo', **params_tru, persistence='')

        # Cache True Distances for the current step
        self.Dc_l_tru = self.cosmo_tru.comovingDistance(z_max=self.zl)
        self.Dc_s_tru = self.cosmo_tru.comovingDistance(z_max=self.zs)
        self.H_tru = self.cosmo_tru.Hz(self.zl)
        self.sigma_crit_tru = self.Dc_s_tru / (self.Dc_l_tru * (self.Dc_s_tru - self.Dc_l_tru))

        # Re-sync cosmological density constants inside the splash object
        self.splash.omg_m = self.Om0_tru
        self.splash.H0 = self.h * 100.0
        self.splash.rho_crt = 3 * self.splash.H0**2 / (8 * np.pi * self.splash.G)
        self.splash.rho_m = self.splash.rho_crt * self.splash.omg_m

        # Scale Line-of-Sight Limit internally
        pi_scale = self.H_fid / self.H_tru
        self.splash.R_max = self.pi_max_fid * pi_scale

        # Calculate target Rsp directly from the fresh splash profile
        self.target_rsp = self.splash.rsp_3d(self.log_R_grid)[0]

        # Resolve Halo Mass and Accretion properties
        self._calculate_halo_properties()
        return self


    def _calculate_halo_properties(self):
        """Matches target_rsp to Diemer17 to extract halo mass dynamically."""
        def rsp_residual(log10_m200m):
            m200m = 10**log10_m200m
            r200m = mass_so.M_to_R(m200m, z=self.zl, mdef='200m') / 1e3
            nu200m = peaks.peakHeight(M=m200m, z=self.zl)
            gamma_diemer = splashback.modelDiemer17Gamma(nu200m, z=self.zl)

            rsp_by_r200m = 0.54 * (1 + 0.53 * self.Om0_tru) * (1 + 1.36 * np.e**(-gamma_diemer/3.04))
            return (rsp_by_r200m * r200m) - self.target_rsp

        try:
            res = root_scalar(rsp_residual, bracket=[7.0, 17.0], method='brentq')
            if res.converged:
                self.log10_m200m = res.root
                self.m200m = 10**self.log10_m200m
                self.r200m = mass_so.M_to_R(self.m200m, z=self.zl, mdef='200m') / 1e3
                self.nu = peaks.peakHeight(M=self.m200m, z=self.zl)
                self.gamma_acc = splashback.modelDiemer17Gamma(self.nu, z=self.zl)
            else:
                self.r200m, self.m200m, self.nu, self.gamma_acc = np.nan, np.nan, np.nan, np.nan
        except ValueError:
            self.r200m, self.m200m, self.nu, self.gamma_acc = np.nan, np.nan, np.nan, np.nan


    def _get_b_hm(self):
        if np.isnan(self.r200m): return np.nan

        r_arr = np.logspace(np.log10(1e-5), np.log10(self.r200m), 1000)
        xi_3d = 10**(self.splash.log_xi_3d(np.log10(r_arr)))
        integrand = xi_3d * 4 * np.pi * (r_arr**2)

        integral_xi_3d = np.trapz(integrand, r_arr)
        numerator = (4.0 / 3.0) * np.pi * (self.r200m**3) * 199.0
        return numerator / integral_xi_3d


    def get_esd_xicg(self, R_bins_ds, R_bins_xicg):
        b_hm = self._get_b_hm()
        if np.isnan(b_hm):
            return np.full_like(R_bins_ds, np.nan), np.full_like(R_bins_xicg, np.nan)

        R_scale = self.Dc_l_tru / self.Dc_l_fid
        R_tru_ds = R_bins_ds * R_scale
        R_tru_xicg = R_bins_xicg * R_scale
        f_crit = self.sigma_crit_fid / self.sigma_crit_tru

        max_R_tru = max(R_tru_ds.max(), R_tru_xicg.max())
        min_R_tru = min(R_tru_ds.min(), R_tru_xicg.min())

        R_min_grid = min(min_R_tru / 2.0, self.splash.spl_rmin)
        R_max_grid = max(max_R_tru * 1.5, self.splash.spl_rmin * 10)

        R_fine = np.logspace(np.log10(R_min_grid), np.log10(R_max_grid), 500)
        log_R_fine = np.log10(R_fine)

        # Draw directly from the natively initialized splash object splines
        xi_2d_fine_tru = 10**self.splash.log_xi_2d(log_R_fine)
        scaled_xi_2d_tru = b_hm * xi_2d_fine_tru

        sigma_R_fine_tru = 2 * self.splash.R_max * self.splash.rho_m * scaled_xi_2d_tru
        integral_fine = cumtrapz(R_fine * sigma_R_fine_tru, R_fine, initial=0)
        integral_fine += (R_fine[0]**2 * sigma_R_fine_tru[0] / 2.0)

        mean_sigma_fine_tru = 2.0 * integral_fine / (R_fine**2)
        delta_sigma_fine_tru = mean_sigma_fine_tru - sigma_R_fine_tru

        ds_interp_tru = ius(np.log10(R_fine), delta_sigma_fine_tru, ext=3)
        xicg_interp_tru = ius(np.log10(R_fine), xi_2d_fine_tru, ext=3)

        delta_sigma_tru_eval = ds_interp_tru(np.log10(R_tru_ds))
        xicg_tru_eval = xicg_interp_tru(np.log10(R_tru_xicg))

        delta_sigma_fid = f_crit * delta_sigma_tru_eval
        xicg_fid = xicg_tru_eval

        return delta_sigma_fid, xicg_fid
# ==========================================
# TEST IMPLEMENTATION
# ==========================================
if __name__ == "__main__":
    # 1. Initialize once
    sim = SplashCosmology(Om0_fid=0.27, sigma8_fid=0.81, zl=0.4, zs=0.6, pi_max_fid=40)

    # 2. Define user bins
    R_bins_ds = np.logspace(-1, 1.5, 30)
    R_bins_xicg = np.logspace(-1, 1.5, 20)

    # 3. Emulate one MCMC step
    proposed_dk14 = {
        'Log_Rho_s': 2.1078, 'Log_Alpha': -0.4621, 'Log_R_s': -0.2033,
        'Log_Rho_0': 0.2106, 'S_e': 0.9873, 'Log_R_t': 0.1939,
        'Log_Beta': 0.8417, 'Log_Gamma': 0.2132, 'F_cen': 1.0, 'R_off': 0.1
    }

    sim.update_mcmc_step(splash_params=proposed_dk14, Om0_tru=0.300, sigma8_tru=0.81)

    # 4. Generate Signal
    if np.isnan(sim.r200m):
        print("Root finder failed, parameters yielded non-physical limits.")
    else:
        ds_model, xicg_model = sim.get_esd_xicg(R_bins_ds, R_bins_xicg)

        # Plot to verify
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        axes[0].plot(R_bins_ds, ds_model/1e12, marker='o', label='$\Delta\Sigma_{fid}$', color='blue')
        axes[0].set(xscale='log', yscale='log', xlabel=r'$R_{fid} \, [h^{-1}\mathrm{Mpc}]$', ylabel=r'$\Delta\Sigma \, [h M_\odot / \mathrm{pc}^2]$')
        axes[0].grid(True, ls="--", alpha=0.5); axes[0].legend()

        axes[1].plot(R_bins_xicg, xicg_model, marker='s', label=r'$\xi_{cg}^{fid}$', color='purple')
        axes[1].set(xscale='log', yscale='log', xlabel=r'$R_{fid} \, [h^{-1}\mathrm{Mpc}]$', ylabel=r'$\xi_{cg}(R)$')
        axes[1].grid(True, ls="--", alpha=0.5); axes[1].legend()

        plt.tight_layout()
        plt.savefig('test.png',dpi=300)
