import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.lss import peaks
from colossus.halo import mass_so, splashback
from scipy.optimize import root_scalar
from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline as ius


# Import the splash class from the provided DK14 model file
from model_dk14 import splash

class splashsim():
    def __init__(self):
        self.h              =   1.0
        self.nspec          =   0.95
        self.Ob0            =   0.049
        self.z_0            =   0.0
        self.tcmb           =   2.726
        self.w              =   -1

    def get_rsp(self, r200m, omgm, gamma):
        '''using eq. 5 from arxiv:1504.05591'''
        rsp_by_r200m = 0.54 * (1 + 0.53 * omgm) * (1 + 1.36 * np.e**(-gamma/3.04))
        return rsp_by_r200m * r200m

    def get_r200m(self, rsp_target, omgm, sigma8, cosmo_name='myCosmo'):
        params = {'flat': True, 'H0': self.h * 100, 'Om0': omgm, 'Ob0': self.Ob0, 'sigma8': sigma8, 'ns': self.nspec, 'Tcmb0': self.tcmb}
        cosmo = cosmology.setCosmology(cosmo_name, **params)

        def rsp_residual(log10_m200m):
            m200m = 10**log10_m200m
            r200m = mass_so.M_to_R(m200m, z=self.z_0, mdef='200m') / 1e3
            nu200m = peaks.peakHeight(M=m200m, z=self.z_0)
            gamma  = splashback.modelDiemer17Gamma(nu200m, z=self.z_0)
            rsp_calc = self.get_rsp(r200m, omgm, gamma)
            return rsp_calc - rsp_target

        try:
            res = root_scalar(rsp_residual, bracket=[7.0, 17.0], method='brentq')
            if not res.converged:
                return np.nan, np.nan
        except ValueError as e:
            return np.nan, np.nan

        log10_m200m_sol = res.root
        m200m_sol = 10**log10_m200m_sol
        r200m_sol = mass_so.M_to_R(m200m_sol, z=self.z_0, mdef='200m') / 1e3

        return r200m_sol, m200m_sol

    def get_b_hm_from_profile(self, r200m, splash_params, cosmo_tru):
        """
        Calculates b_hm by enforcing M_excess(<r200m) = 199 * M_mean(<r200m).
        """
        sp = splash(**splash_params)

        def integrand(r):
            log_r = np.log10(max(r, 1e-6))
            xi_3d = 10**(sp.log_xi_3d(log_r))
            return xi_3d * 4 * np.pi * (r**2)

        integral_xi_3d, _ = quad(integrand, 1e-5, r200m, limit=100)
        numerator = (4.0 / 3.0) * np.pi * (r200m**3) * 199.0
        b_hm = numerator / integral_xi_3d

        return b_hm

    def get_delta_sigma_scaled(self, R_fid_arr, splash_params, b_hm, zl, zs, zeff, pi_max_fid, cosmo_tru, cosmo_fid):
        """
        Computes Delta Sigma scaled from the true theory cosmology to a fiducial observer cosmology.
        Uses scaling relations from More et al. 2013 (arxiv:1309.2943).
        """
        # --- 1. Compute Scaling Factors ---
        # 1a. Radial scaling: R_tru = R_fid * (Dc_tru / Dc_fid)
        Dc_l_tru = cosmo_tru.comovingDistance(z_max=zl)
        Dc_l_fid = cosmo_fid.comovingDistance(z_max=zl)
        Dc_s_tru = cosmo_tru.comovingDistance(z_max=zs)
        Dc_s_fid = cosmo_fid.comovingDistance(z_max=zs)

        R_scale = Dc_l_tru / Dc_l_fid
        R_tru_arr = R_fid_arr * R_scale

        # 1b. Pi_max scaling: pi_tru = pi_fid * (H_fid / H_tru)
        H_tru = cosmo_tru.Hz(zeff)
        H_fid = cosmo_fid.Hz(zeff)
        pi_scale = H_fid / H_tru
        pi_max_tru = pi_max_fid * pi_scale

        # 1c. Amplitude scaling: f_crit = Sigma_crit_tru / Sigma_crit_fid
        # we are working in the comoving coordinates
        # Sigma_crit is proportional to Dc_s / (Dc_l * (Dc_s - Dc_l))
        sigma_crit_tru = Dc_s_tru / (Dc_l_tru * (Dc_s_tru - Dc_l_tru))
        sigma_crit_fid = Dc_s_fid / (Dc_l_fid * (Dc_s_fid - Dc_l_fid))

        f_crit = sigma_crit_tru / sigma_crit_fid

        # --- 2. Evaluate Theoretical Model in True Cosmology ---
        # Update splash params with the True pi_max (which maps to R_max in the dk14 code)
        splash_params['R_max'] = pi_max_tru
        sp = splash(**splash_params)

        # Override the constants inside splash to perfectly match the true cosmology
        # (Otherwise it uses the hardcoded H0=100, omg_m=0.315 from model_dk14.constants)
        sp.H0 = cosmo_tru.H0
        sp.omg_m = cosmo_tru.Om0
        sp.rho_crt = 3 * sp.H0**2 / (8 * np.pi * sp.G)
        sp.rho_m = sp.rho_crt * sp.omg_m

        R_max_grid = max(R_tru_arr.max(), sp.spl_rmin * 10)
        R_fine = np.logspace(np.log10(sp.spl_rmin), np.log10(R_max_grid), 500)
        log_R_fine = np.log10(R_fine)

        xi_2d_fine = 10**sp.log_xi_2d(log_R_fine)
        scaled_xi_2d = b_hm * xi_2d_fine

        # Sigma in true physical/comoving units
        sigma_R_fine_tru = 2 * sp.R_max * sp.rho_m * scaled_xi_2d

        integral_fine = cumtrapz(R_fine * sigma_R_fine_tru, R_fine, initial=0)
        integral_fine += (R_fine[0]**2 * sigma_R_fine_tru[0] / 2.0)

        mean_sigma_fine_tru = 2.0 * integral_fine / (R_fine**2)
        delta_sigma_fine_tru = mean_sigma_fine_tru - sigma_R_fine_tru

        # Interpolate the TRUE profile
        sigma_interp_tru    = ius(np.log10(R_fine), sigma_R_fine_tru)
        ds_interp_tru       = ius(np.log10(R_fine), delta_sigma_fine_tru)

        # Evaluate at the true mapped radii
        sigma_tru_eval = sigma_interp_tru(np.log10(R_tru_arr))
        delta_sigma_tru_eval = ds_interp_tru(np.log10(R_tru_arr))

        # --- 3. Scale Back to Fiducial Observable ---
        sigma_fid = f_crit * sigma_tru_eval
        delta_sigma_fid = f_crit * delta_sigma_tru_eval

        return sigma_fid, delta_sigma_fid


if __name__ == "__main__":
    sim = splashsim()

    # Define True vs Fiducial Cosmologies
    # e.g., Simulation cosmology vs WMAP9/Planck fiducial
    params_tru = {'flat': True, 'H0': 100.0, 'Om0': 0.300, 'Ob0': 0.049, 'sigma8': 0.81, 'ns': 0.95, 'Tcmb0': 2.726}
    params_fid = {'flat': True, 'H0': 100.0, 'Om0': 0.315, 'Ob0': 0.049, 'sigma8': 0.81, 'ns': 0.95, 'Tcmb0': 2.726}

    cosmo_tru = cosmology.setCosmology('true_cosmo', **params_tru)
    cosmo_fid = cosmology.setCosmology('fid_cosmo', **params_fid)

    # Observational Redshifts
    z_l   = 0.3   # Lens plane redshift
    z_s   = 0.7   # Source plane redshift
    z_eff = 0.3   # Effective redshift for clustering/correlation function

    # Setup
    target_rsp = 1.5
    pi_max_fid = 40.0 # Standard integration limit in Mpc/h

    # 1. Get r200m in the TRUE cosmology
    r200m_sol, m200m_sol = sim.get_r200m(target_rsp, params_tru['Om0'], params_tru['sigma8'], cosmo_name='true_cosmo')

    print(f"Solved R_200m (True): {r200m_sol:.3f} Mpc/h")
    print(f"Solved M_200m (True): {m200m_sol:.2e} M_sun/h")

    # Sample DK14 profile parameters
    dk14_params = {
        'Log_Rho_s': 1.47611335, 'Log_Alpha': 0.68742786, 'Log_R_s': 0.67286521,
        'Log_Rho_0': 0.31423027, 'S_e': 2.07541219, 'Log_R_t': 0.04709024,
        'Log_Beta': 0.89337944, 'Log_Gamma': 0.13532044, 'F_cen': 0.24619855,
        'R_off': 0.33717233, 'splrmin': 0.01, 'splrbin': 150
    }

    # 2. Calculate bias in True cosmology
    computed_b_hm = sim.get_b_hm_from_profile(r200m_sol, dk14_params, cosmo_tru)
    print(f"Enforced Halo Bias (b_hm): {computed_b_hm:.3f}")

    # 3. Compute scaled Weak Lensing Delta Sigma
    R_fid_projected = np.logspace(-1, 1.5, 50)  # User requested radii in fiducial frame

    sigma_fid, ds_fid = sim.get_delta_sigma_scaled(
        R_fid_arr=R_fid_projected,
        splash_params=dk14_params,
        b_hm=computed_b_hm,
        zl=z_l, zs=z_s, zeff=z_eff,
        pi_max_fid=pi_max_fid,
        cosmo_tru=cosmo_tru,
        cosmo_fid=cosmo_fid
    )

    # Plotting
    plt.subplot(2,2,1)
    plt.plot(R_fid_projected, ds_fid/1e12, label=f'$\Delta\Sigma_{{fid}}(R_{{fid}})$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$R_{fid} \, [h^{-1}\mathrm{Mpc}]$', fontsize=14)
    plt.ylabel(r'$\Delta\Sigma \, [h M_\odot / \mathrm{pc}^2]$ (Fiducial Frame)', fontsize=14)
    plt.title('Cosmology-Scaled Weak Lensing Signal', fontsize=14)
    plt.legend(fontsize=12)
    plt.savefig("delta_sigma_scaled.png", dpi=300)

