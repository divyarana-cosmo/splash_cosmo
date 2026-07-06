import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.lss import peaks
from colossus.halo import mass_so, splashback
from scipy.optimize import root_scalar
from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline as ius

# Import the base splash class from your DK14 model file
import sys
sys.path.append('/home/rana/github_0/splash_cosmo/')
from model_dk14 import splash

class SplashCosmology(splash):
    """
    MCMC-optimized wrapper for the DK14 Splashback model.
    Handles fiducial vs. true cosmological mappings, accretion rate
    computations, and dynamic shape integrations in-place.
    """
    def __init__(self, Om0_fid=0.27, sigma8_fid=0.81, zl=0.4, zs=0.6, pi_max_fid=40,
                 h=1.0, nspec=0.95, Ob0=0.049, tcmb=2.726):

        # 1. Dummy Initialization for Parent Class
        # We pass the requested pi_max_fid to ensure the baseline spline
        # grid accommodates the requested projection length.
        dummy_dk14 = {
            'Log_Rho_s': 1.0, 'Log_Alpha': -0.5, 'Log_R_s': -0.2,
            'Log_Rho_0': 0.2, 'S_e': 1.0, 'Log_R_t': 0.2,
            'Log_Beta': 0.8, 'Log_Gamma': 0.2, 'F_cen': 1.0, 'R_off': 0.1,
            'R_max': pi_max_fid
        }
        super().__init__(**dummy_dk14)

        # 2. Store Base Hyperparameters
        self.zl = zl
        self.zs = zs
        self.pi_max_fid = pi_max_fid
        self.h = h
        self.nspec = nspec
        self.Ob0 = Ob0
        self.tcmb = tcmb

        # 3. Setup and Permanently Cache the Fiducial Cosmology
        self.Om0_fid = Om0_fid
        self.sigma8_fid = sigma8_fid
        params_fid = {
            'flat': True, 'H0': 100.0 * self.h, 'Om0': self.Om0_fid,
            'Ob0': self.Ob0, 'sigma8': self.sigma8_fid, 'ns': self.nspec, 'Tcmb0': self.tcmb
        }
        self.cosmo_fid = cosmology.setCosmology('fid_cosmo', **params_fid, persistence='')

        # Cache invariant fiducial distances to avoid repeating in MCMC
        self.Dc_l_fid = self.cosmo_fid.comovingDistance(z_max=self.zl)
        self.Dc_s_fid = self.cosmo_fid.comovingDistance(z_max=self.zs)
        self.H_fid = self.cosmo_fid.Hz(self.zl)
        self.sigma_crit_fid = self.Dc_s_fid / (self.Dc_l_fid * (self.Dc_s_fid - self.Dc_l_fid))

        # Log grid used strictly for finding the target 3D Rsp derivative
        self.log_R_grid = np.log10(np.logspace(-1, 1, 100))


    def update_mcmc_step(self, splash_params, Om0_tru, sigma8_tru):
        """
        Updates DK14 shape parameters and true cosmology in-place.
        Highly optimized to prevent memory reallocation during sampling.
        """
        # --- 1. Update Inherited DK14 Shape Parameters ---
        self.rho_s = 10**splash_params['Log_Rho_s']
        self.alpha = 10**splash_params['Log_Alpha']
        self.r_s   = 10**splash_params['Log_R_s']
        self.rho_0 = 10**splash_params['Log_Rho_0']
        self.s_e   = splash_params['S_e']
        self.r_t   = 10**splash_params['Log_R_t']
        self.beta  = 10**splash_params['Log_Beta']
        self.gamma = 10**splash_params['Log_Gamma']
        if 'F_cen' in splash_params: self.fcen = splash_params['F_cen']
        if 'R_off' in splash_params: self.roff = splash_params['R_off']

        # --- 2. Recalculate Target 3D Rsp based on new shape ---
        # The steepest slope of the new 3D profile
        self.target_rsp = self.rsp_3d(self.log_R_grid)[0]

        # --- 3. Update True Cosmology ---
        self.Om0_tru = Om0_tru
        self.sigma8_tru = sigma8_tru
        params_tru = {
            'flat': True, 'H0': 100.0 * self.h, 'Om0': self.Om0_tru,
            'Ob0': self.Ob0, 'sigma8': self.sigma8_tru, 'ns': self.nspec, 'Tcmb0': self.tcmb
        }
        # persistence='' keeps memory clean during thousands of MCMC steps
        self.cosmo_tru = cosmology.setCosmology('true_cosmo', **params_tru, persistence='')

        # Cache True Distances for the current step
        self.Dc_l_tru = self.cosmo_tru.comovingDistance(z_max=self.zl)
        self.Dc_s_tru = self.cosmo_tru.comovingDistance(z_max=self.zs)
        self.H_tru = self.cosmo_tru.Hz(self.zl)
        self.sigma_crit_tru = self.Dc_s_tru / (self.Dc_l_tru * (self.Dc_s_tru - self.Dc_l_tru))

        # Update base densities for the inherited DK14 model
        self.H0 = self.h * 100.0
        self.omg_m = self.Om0_tru
        self.rho_crt = 3 * self.H0**2 / (8 * np.pi * self.G)
        self.rho_m = self.rho_crt * self.omg_m

        # --- 4. Scale Line-of-Sight Limit and Clear Cache ---
        # pi_tru = pi_fid * (H_fid / H_tru)
        pi_scale = self.H_fid / self.H_tru
        self.R_max = self.pi_max_fid * pi_scale

        # Force the inherited splash class to rebuild integration splines
        # up to the newly scaled true pi_max limit
        self.init_spl = False
        self.init_spl_cen = False
        self.init_spl_2d_diff = False
        self.init_spl_3d_diff = False

        # --- 5. Resolve Halo Mass and Accretion properties ---
        self._calculate_halo_properties()
        return self


    def _calculate_halo_properties(self):
        """
        Dynamically infers m200m, r200m, peak height (nu), and mass accretion (Gamma)
        by matching the target_rsp to the Diemer17 empirical scaling relation.
        """
        def rsp_residual(log10_m200m):
            m200m = 10**log10_m200m
            # Colossus returns r200m in kpc/h, convert to Mpc/h
            r200m = mass_so.M_to_R(m200m, z=self.zl, mdef='200m') / 1e3
            nu200m = peaks.peakHeight(M=m200m, z=self.zl)
            gamma_diemer = splashback.modelDiemer17Gamma(nu200m, z=self.zl)

            # Using eq. 5 from arxiv:1504.05591
            rsp_by_r200m = 0.54 * (1 + 0.53 * self.Om0_tru) * (1 + 1.36 * np.e**(-gamma_diemer/3.04))
            return (rsp_by_r200m * r200m) - self.target_rsp

        try:
            # Bracket [7.0, 17.0] safely covers 10^7 to 10^17 solar masses
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
        """Replaced adaptive 'quad' with a fixed, ultra-dense trapz grid to eliminate derivative noise."""
        if np.isnan(self.r200m): return np.nan
        
        # 1000 fixed points ensures perfect smoothness for numerical derivatives
        r_arr = np.logspace(np.log10(1e-5), np.log10(self.r200m), 1000)
        xi_3d = 10**(self.log_xi_3d(np.log10(r_arr)))
        integrand = xi_3d * 4 * np.pi * (r_arr**2)

        integral_xi_3d = np.trapz(integrand, r_arr)
        numerator = (4.0 / 3.0) * np.pi * (self.r200m**3) * 199.0
        return numerator / integral_xi_3d

    def log_xi_2d_cen(self, log_R):
        """
        OVERRIDE PARENT CLASS: 
        The parent class uses a noisy Riemann sum (np.sum) with 1000 points. 
        When R_max scales with cosmology, the grid shifts and destroys the 
        Fisher finite-differences. We override it here with smooth np.trapz.
        """
        from scipy.interpolate import UnivariateSpline
        if self.fcen == 1.0:
            rrad = np.logspace(log_R.min(), log_R.max(), self.spl_rbin)
            _val = np.zeros_like(rrad)
        else:
            rrad = np.logspace(np.log10(self.spl_rmin), np.log10(10**log_R.max() + 10*self.roff), self.spl_rbin)
            _val = np.zeros_like(rrad)
            
        # Use a higher resolution trapezoidal integration
        xx = np.linspace(0, self.R_max, 2500) 
        for ii, rr in enumerate(rrad):
            integrand = 10**(self.log_xi_3d(np.log10(np.sqrt(rr**2 + xx**2))))
            _val[ii] = np.trapz(integrand, xx)
            
        val = _val / self.R_max
        
        self.spl_log_xi_2d_cen = UnivariateSpline(np.log10(rrad), np.log10(val), s=0.0, k=4)
        self.spl_log_xi_2d_cen_diff = self.spl_log_xi_2d_cen.derivative()

        self.log_xi_2d_dict = {}
        self.log_xi_2d_dict['logr_min'] = np.log10(np.min(rrad))
        self.log_xi_2d_dict['logr_max'] = np.log10(np.max(rrad))
        self.log_xi_2d_dict['log_ximax'] = np.log10(np.max(val))
        self.init_spl_cen = True
        
        return self.spl_log_xi_2d_cen(log_R)

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
        
        R_min_grid = min(min_R_tru / 2.0, self.spl_rmin)
        R_max_grid = max(max_R_tru * 1.5, self.spl_rmin * 10)
        
        # UPGRADED: 500 points instead of 100 to stop spline jitter
        R_fine = np.logspace(np.log10(R_min_grid), np.log10(R_max_grid), 500)
        log_R_fine = np.log10(R_fine)

        xi_2d_fine_tru = 10**self.log_xi_2d(log_R_fine)
        scaled_xi_2d_tru = b_hm * xi_2d_fine_tru
        
        sigma_R_fine_tru = 2 * self.R_max * self.rho_m * scaled_xi_2d_tru
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
