from dark_emulator import darkemu
from scipy.interpolate import InterpolatedUnivariateSpline as ius
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.optimize import minimize
from scipy.signal import savgol_filter

emu = darkemu.base_class()

class splashsim():
    def __init__(self, H0=67.5, Om0=0.3, Ob0=0.022, Tcmb0=2.7255, Neff=3.046, sigma8=0.82, ns=0.96):
        """
        Initialize the splashback simulation with cosmological parameters.

        Parameters:
        - H0: Hubble constant in km/s/Mpc (e.g., 67.5)
        - Om0: Matter density parameter (e.g., 0.3)
        - Ob0: BARYON DENSITY PARAMETER in ωb form (Ωb·h²), not Ωb directly (e.g., 0.022)
        - Other standard cosmological parameters
        """
        # Store the full Hubble constant and density parameters
        self.H0 = 100
        self.h = H0/100.0  # Dimensionless Hubble parameter
        self.Om0 = Om0     # Matter density parameter

        # For dark emulator, we need: (ωb, ωc, Ωde, ln(10^10As), ns, w)
        # Where ωb ≡ Ωb·h² and ωc ≡ Ωc·h² are physical densities

        # Ob0 is already provided as ωb = Ωb·h²
        wb = Ob0

        # Calculate ωc = Ωc·h² = (Ωm - Ωb)·h²
        # Since Om0 is Ωm and Ωb = wb/h²
        Ob_true = wb/self.h**2  # True Ωb value
        wc = (Om0 - Ob_true)*self.h**2

        # Neutrino density parameter (approximate)
        wv = 0.00064

        # Dark energy density (assuming flat universe)
        omg_de = 1.0 - Om0 - wv/self.h**2

        # Other cosmological parameters
        ln1e10As = 3.094  # Amplitude of primordial fluctuations
        self.ns = ns      # Spectral index
        w = -1            # Dark energy equation of state

        # Set cosmology in the emulator
        cparam = np.array([wb, wc, omg_de, ln1e10As, ns, w])
        emu.set_cosmology(cparam)

        self.init_spl_M200m2nu = False
        self.G = 4.301e-9  # km^2 Mpc M_sun^-1 s^-2 gravitational constant

    def M200m2r200m(self, logmh):
        """
        Convert M200m to r200m (radius containing 200x mean density)

        Parameters:
        - logmh: log10 of the halo mass M200m in h^-1 solar masses

        Returns:
        - r200m: radius in h^-1 Mpc
        """
        # Input is already in h^-1 M_sun, which is what we want
        m_tot = 10**logmh  # total mass in h^-1 M_sun units

        # Calculate critical density in h^2 M_sun/(h^-1 Mpc)^3
        # ρ_crit = 3H²/(8πG) expressed in h^2 M_sun/(h^-1 Mpc)^3
        rho_crit_h = 3*(self.H0)**2/(8*np.pi*self.G)  # in h^2 M_sun/(h^-1 Mpc)^3

        # Calculate mean density of the universe in h^2 M_sun/(h^-1 Mpc)^3
        # ρ_mean = Ωm * ρ_crit
        rho_mean_h = self.Om0 * rho_crit_h  # in h^2 M_sun/(h^-1 Mpc)^3

        # For h^-1 Mpc units, we need to adjust the density calculation
        # Since m_tot is in h^-1 M_sun, and we want r_200m in h^-1 Mpc:
        # M200m [h^-1 M_sun] = (4π/3) * 200 * ρ_mean [h^2 M_sun/(h^-1 Mpc)^3] * (r200m [h^-1 Mpc])^3 / h^2
        # Solving for r200m:
        r_200m = (3*m_tot/(4*np.pi*200*rho_mean_h))**(1./3.)  # in h^-1 Mpc

        return r_200m

    def _spl_M200m2nu(self, z):
        """Initialize spline for peak height calculation"""
        # dark_emulator uses h^-1 M_sun units for mh
        mh, sigmh, mul_func = emu.get_f_HMF(z)
        # Create spline - no need to convert units as emulator already works in h^-1 M_sun
        self.spl_logMh2sigM = ius(np.log10(mh), sigmh)
        self.init_spl_M200m2nu = True
        self.init_redshift = z

    def logM200m2nu(self, logmh, z):
        """
        Convert M200m to peak height nu

        Parameters:
        - logmh: log10 of halo mass in h^-1 solar masses
        - z: redshift

        Returns:
        - nu: peak height
        """
        if not self.init_spl_M200m2nu or z != self.init_redshift:
            self._spl_M200m2nu(z)

        # No need to convert units - dark_emulator already works with h^-1 M_sun
        self.delta_crit = 1.686
        #return self.delta_crit/self.spl_logMh2sigM(logmh)
        return self.spl_logMh2sigM(logmh)

    def logM200m2rsp(self, logmh, z):
        """
        Find splashback radius by locating minimum of d(log xi)/d(log r)

        Parameters:
        - logmh: log10 of halo mass in h^-1 solar masses
        - z: redshift

        Returns:
        - Splashback radius in h^-1 Mpc
        """
        # No need to convert mass units - dark_emulator works with h^-1 M_sun
        M = 10**logmh  # Mass in h^-1 M_sun

        # Create a more finely sampled radius array in h^-1 Mpc
        rs_fine = np.logspace(-1, 1, 500)  # h^-1 Mpc

        # Get halo-matter cross-correlation function
        # dark_emulator expects mass in h^-1 M_sun and radius in h^-1 Mpc
        xihm = emu.get_xicross_mass(rs_fine, M, z)

        # Convert to log space
        log_rs = np.log10(rs_fine)
        log_xihm = np.log10(xihm)

        # Calculate derivative d(log xi)/d(log r)
        yy = np.diff(log_xihm)/np.diff(log_rs)
        xx = 0.5 * (log_rs[1:] + log_rs[:-1])

        # Apply stronger smoothing to get stable derivatives
        window_size = min(15, len(yy) - 2)  # Make sure window size is odd and less than array length
        if window_size % 2 == 0:
            window_size -= 1

        y_smooth = savgol_filter(yy, window_size, 3)

        # Look only at the range where we expect the splashback radius
        # Typically the splashback radius is around 1-3 times r200m
        r200m = self.M200m2r200m(logmh)  # in h^-1 Mpc

        # Find indices corresponding to reasonable range (0.5-5 times r200m)
        min_r_idx = np.argmin(np.abs(10**xx - 0.1*r200m))
        max_r_idx = np.argmin(np.abs(10**xx - 2.0*r200m))

        # Restrict search to this range
        search_range = slice(min_r_idx, max_r_idx)

        # Find the minimum in the smoothed derivative within this range
        min_idx = min_r_idx + np.argmin(y_smooth[search_range])

        # Additional robustness: check if the minimum looks reliable by ensuring it's actually negative
        if y_smooth[min_idx] > -0.5:  # If not sufficiently negative, may not be a real feature
            # Try to find any reasonable minimum in the expected range
            fallback_min_idx = min_r_idx + np.argmin(y_smooth[search_range])
            if y_smooth[fallback_min_idx] < y_smooth[min_idx]:
                min_idx = fallback_min_idx

        # Return the radius corresponding to this minimum (in h^-1 Mpc)
        rsp = 10**xx[min_idx]

        # Debug information
        print(f"M={logmh:.2f} [h^-1 M_sun], r200m={r200m:.3f} [h^-1 Mpc], rsp={rsp:.3f} [h^-1 Mpc], rsp/r200m={rsp/r200m:.3f}, min_slope={y_smooth[min_idx]:.3f}")

        return rsp

def plot_rsp_vs_peak_height_varying_omega():
    """
    Plot the relationship between splashback radius and peak height
    for different values of Omega_m.

    This function compares the splashback radius calculation from
    our implementation with the Colossus package.
    """
    # COSMOLOGICAL PARAMETERS - CLEARLY DEFINED
    h = 1.0        # Dimensionless Hubble parameter (H0/100)
    H0 = 100.0 * h   # Full Hubble constant in km/s/Mpc

    # Baryon density in ωb form (Ωb·h²)
    omega_b = 0.022  # This is ωb = Ωb·h², NOT Ωb directly

    # List of matter density parameters to test
    omega_m_values = [0.27]  # Can add more like [0.28, 0.29, 0.3, 0.31]

    # Import colossus for comparison calculations
    from colossus.cosmology import cosmology
    from colossus.halo import mass_so
    from colossus.lss import peaks
    from colossus.halo import splashback

    redshift = 0.0  # Fixed redshift for all calculations

    # Create a figure
    plt.figure(figsize=(15, 12))

    # Set up plotting areas
    ax = plt.subplot(3, 3, 1)  # Cross-correlation function
    ax1 = plt.subplot(3, 3, 2)  # Derivative of cross-correlation
    ax2 = plt.subplot(3, 3, 3)  # Radius comparison
    ax3 = plt.subplot(3, 3, 4)  # Rsp/R200m ratio
    ax4 = plt.subplot(3, 3, 5)  # Peak height vs mass
    ax5 = plt.subplot(3, 3, 6)  # Splashback radius vs mass

    for om in omega_m_values:
        print(f"--- Processing Omega_m = {om} ---")

        # Initialize our splashback model with this omega_m
        # Note: we pass omega_b (ωb = Ωb·h²) directly as Ob0
        sim = splashsim(H0=H0, Om0=om, Ob0=omega_b)

        # Set up colossus cosmology for comparison
        # Note: Colossus expects Ωb (not ωb), so we need to convert
        Ob_true = omega_b/h**2  # Convert ωb to Ωb

        params = {
            'flat': True,         # Assume flat universe
            'H0': 100,             # Full Hubble constant
            'Om0': om,            # Matter density parameter
            'Ob0': Ob_true,       # Baryon density parameter (NOT in ωb form)
            'sigma8': 0.82,       # Power spectrum normalization
            'ns': 0.96            # Spectral index
        }

        # Add and set the cosmology in Colossus
        cosmology_name = f'cosmo_Om{om}'
        cosmology.addCosmology(cosmology_name, **params)
        cosmo = cosmology.setCosmology(cosmology_name)

        print(f"  Cosmology set: H0={H0}, Om0={om}, Ob0={Ob_true}")
        print(f"  Critical density: {cosmo.rho_c(0):.3e} M_sun/Mpc^3")

        # Range of halo masses (log10 scale) in h^-1 M_sun
        logm_values = np.linspace(12.5, 15.0, 5)
        nu_values = []
        rsp_values = []
        r200m_values = []

        for logm in logm_values:
            M = 10**logm  # Mass in h^-1 M_sun

            # Calculate halo-matter cross correlation
            # The emulator expects distances in h^-1 Mpc
            rs = np.logspace(-1, 1, 100)  # h^-1 Mpc
            xihm = emu.get_xicross_mass(rs, M, redshift)

            # Plot cross-correlation
            ax.plot(rs, xihm, label=f'$\\log M_{{h}} = {logm:.2f}$ $h^{{-1}}M_\\odot$')

            # Calculate and plot derivative
            log_rs = np.log10(rs)
            log_xihm = np.log10(xihm)
            yy = np.diff(log_xihm)/np.diff(log_rs)
            xx = 0.5 * (log_rs[1:] + log_rs[:-1])

            # Apply stronger smoothing to stable derivatives
            window_size = min(11, len(yy) - 2)  # Make sure window size is odd and less than array length
            if window_size % 2 == 0:
                window_size -= 1

            y_smooth = savgol_filter(yy, window_size, 3)

            # Plot the smoothed derivative
            ax1.plot(10**xx, y_smooth, '-')

            # Calculate peak height nu
            nu = sim.logM200m2nu(logm, redshift)
            nu_values.append(nu)

            # Calculate splashback radius
            rsp = sim.logM200m2rsp(logm, redshift)  # in h^-1 Mpc
            ax1.axvline(rsp, linestyle='--', alpha=0.5)
            rsp_values.append(rsp)

            # Calculate r200m
            r200m = sim.M200m2r200m(logm)  # in h^-1 Mpc
            r200m_values.append(r200m)

        # Convert lists to arrays for easier handling
        nu_values       = np.array(nu_values)
        rsp_values      = np.array(rsp_values)
        r200m_values    = np.array(r200m_values)

        # Get r200m from colossus for comparison
        # Colossus expects masses in M_sun and returns radii in kpc - need to adjust
        colossus_masses = 10**np.array(logm_values) * h  # Convert from h^-1 M_sun to M_sun
        # Colossus returns radii in kpc, so we convert to h^-1 Mpc:
        # First to Mpc (divide by 1000), then to h^-1 Mpc (multiply by h)
        R200m_colossus = mass_so.M_to_R(colossus_masses, redshift, '200m')/1000.0 * h  # Convert to h^-1 Mpc

        # Plot radii comparison
        ax2.plot(nu_values, rsp_values, '-o', label='$r_{sp}$')
        ax2.plot(nu_values, r200m_values, '-s', label='$r_{200m}$')
        ax2.plot(nu_values, R200m_colossus, '-^', label='colossus $r_{200m}$')

        # Plot ratio of splashback to r200m
        ax3.plot(nu_values, rsp_values/r200m_values, '-o', label=f'$\\Omega_m = {om}$')
        # Compare with More et al. 2015 model
        ax3.plot(nu_values, splashback.modelMore15RspR200m(nu200m=nu_values, z=redshift, statistic='median'),
                 '--', label='More+15')

        # Plot peak height calculations
        # Colossus expects masses in M_sun
        #nu_colossus = peaks.peakHeight(colossus_masses, z=redshift)
        nu_colossus = cosmo.sigma(r200m_values, z=redshift)
        ax4.plot(logm_values, nu_values, '-o', label='Our calc')
        ax4.plot(logm_values, nu_colossus, '--s', label='colossus')

        # Plot splashback radius vs mass
        ax5.plot(logm_values, rsp_values, '-o', label='$r_{sp}$')
        ax5.plot(logm_values, r200m_values, '--s', label='$r_{200m}$')

    # Set axis properties
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$r$ [$h^{-1}$Mpc]')
    ax.set_ylabel(r'$\xi_{hm}(r)$')
    ax.set_title('Halo-matter cross-correlation')

    ax1.set_xscale('log')
    ax1.set_xlabel('$r$ [$h^{-1}$Mpc]')
    ax1.set_ylabel(r'$d\log \xi_{hm} / d \log r$')
    ax1.set_title('Derivative of cross-correlation')

    ax2.set_xlabel('$\\nu$')
    ax2.set_ylabel(r'$r$ [$h^{-1}$Mpc]')
    ax2.set_title('Radii vs. peak height')

    ax3.set_xlabel('$\\nu$')
    ax3.set_ylabel(r'$r_{sp}/r_{200m}$')
    ax3.set_title('Splashback to $r_{200m}$ ratio')

    ax4.set_ylabel('$\\nu$')
    ax4.set_xlabel(r'$\log M_{200m}$ [$h^{-1}M_\odot$]')
    ax4.set_title('Peak height vs. mass')

    ax5.set_ylabel('$r$ [$h^{-1}$Mpc]')
    ax5.set_xlabel(r'$\log M_{200m}$ [$h^{-1}M_\odot$]')
    ax5.set_title('Radii vs. mass')

    # Add legends
    ax.legend(fontsize='xx-small')
    ax2.legend(fontsize='small')
    ax3.legend(fontsize='small')
    ax4.legend(fontsize='small')
    ax5.legend(fontsize='small')

    plt.tight_layout()
    plt.savefig('rsp_vs_nu_varying_omega.png', dpi=300)

if __name__ == "__main__":
    plot_rsp_vs_peak_height_varying_omega()
