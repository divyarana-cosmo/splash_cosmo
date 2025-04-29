from dark_emulator import darkemu
from scipy.interpolate import InterpolatedUnivariateSpline as ius
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.optimize import minimize

emu = darkemu.base_class()

class splashsim():
    def __init__(self, H0=67.5, Om0=0.3, Ob0=0.022, Tcmb0=2.7255, Neff=3.046, sigma8=0.82, ns=0.96):
        # fixing the cosmology
        self.H0 = H0
        self.Om0 = Om0
        self.Ob0 = Ob0
        # (ωb, ωc, Ωde, ln(10^10As), ns, w)
        # ωb ≡ Ωbh2 and ωc ≡ Ωch2
        h = (H0/100)
        wb = Ob0  # Ob0 is already in the form of ωb (Ωbh²)
        wc = Om0*h**2 - wb  # Corrected: ωc = Ωmh² - ωb
        wv = 0.00064
        omg_de = 1 - (Om0 + wv/h**2)  # Simplified calculation
        ln1e10As = 3.094
        self.ns = ns
        w = -1
        cparam = np.array([wb, wc, omg_de, ln1e10As, ns, w])
        emu.set_cosmology(cparam)

        self.init_spl_M200m2nu = False
        self.G = 4.301e-9  # km^2 Mpc M_sun^-1 s^-2 gravitational constant

    def M200m2r200m(self, logmh):
        m_tot = 10**logmh  # total mass of the halo
        rho_crt = 3*self.H0**2/(8*np.pi*self.G)  # rho critical
        r_200m = (3*m_tot/(4*np.pi*200*rho_crt*self.Om0))**(1./3.)
        return r_200m

    def _spl_M200m2nu(self, z):
        mh, sigmh, mul_func = emu.get_f_HMF(z)
        self.spl_logMh2sigM = ius(np.log10(mh), sigmh/emu.Dgrowth_from_z(z))
        self.init_spl_M200m2nu = True
        self.init_redshift = z

    def logM200m2nu(self, logmh, z):
        if not self.init_spl_M200m2nu or z != self.init_redshift:
            self._spl_M200m2nu(z)
        self.delta_crit = 1.686
        return self.delta_crit/self.spl_logMh2sigM(logmh)

    def logM200m2rsp(self, logmh, z):
        M = 10**logmh
        rs = np.logspace(-1, 1, 100)
        xihm = emu.get_xicross_mass(rs, M, z)

        log_rs = np.log10(rs)
        log_xihm = np.log10(xihm)

        yy = (log_xihm[1:] - log_xihm[:-1])/(log_rs[1] - log_rs[0])
        xx = log_rs[1:]*0.5 + log_rs[:-1]*0.5
        from scipy.signal import savgol_filter
        y_smooth = savgol_filter(yy, 5, 3)


        # Create a spline of log(xihm) vs log(rs)
        spl = ius(xx, y_smooth,k=3)
        # Evaluate the derivative at these points
        deriv_values = spl(xx)
        # Find the index of the minimum derivative value
        min_idx = np.argmin(deriv_values)
        print(10**xx[min_idx], deriv_values[min_idx])

        res = minimize(spl, x0=xx[min_idx], method='Nelder-Mead', tol=1e-6)
        return 10**res.x

def plot_rsp_vs_peak_height_varying_omega():
    # Fixed baryon density within supported range
    Ob0 = 0.022  # Already in omega_b form (Ωbh²)
    #omega_m_values = [0.28, 0.29, 0.3, 0.31]
    omega_m_values = [0.3]
    h = 0.675  # Use h=0.675 as a baseline

    redshift = 0.0  # Fixed redshift
    plt.subplot(3,3,1)
    for om in omega_m_values:
        # Initialize model with this omega_m and compatible parameters
        sim = splashsim(H0=100*h, Om0=om, Ob0=Ob0)

        # Range of halo masses
        logm_values = np.linspace(12.5, 15.0, 5)
        nu_values = np.array([])
        rsp_r200m_values =np.array([])
        r200m_values =np.array([])

        for logm in logm_values:
            rs = np.logspace(-1, 1, 100)
            xihm = emu.get_xicross_mass(rs, 10**logm, redshift)
            ax = plt.subplot(3,3,1)
            ax.plot(rs, xihm,label=r'$\log M_{h} = %2.2f$'%logm)

            ax1  = plt.subplot(3,3,2)
            log_rs = np.log10(rs)
            log_xihm = np.log10(xihm)

            yy = (log_xihm[1:] - log_xihm[:-1])/(log_rs[1] - log_rs[0])
            xx = log_rs[1:]*0.5 + log_rs[:-1]*0.5
            from scipy.signal import savgol_filter
            y_smooth = savgol_filter(yy, 5, 3)


            # Create a spline of log(xihm) vs log(rs)
            spl = ius(xx, y_smooth,k=3)
            ax1.plot(10**xx,spl(xx))
            # Calculate peak height nu
            nu = sim.logM200m2nu(logm, redshift)
            nu_values=np.append(nu_values,nu)

            ## Calculate splashback radius
            rsp = sim.logM200m2rsp(logm, redshift)
            ax1.axvline(rsp)

            # Calculate r200m
            r200m = sim.M200m2r200m(logm)

            # Store ratio
            #rsp_r200m_values=np.append(rsp_r200m_values,rsp/r200m)
            rsp_r200m_values=np.append(rsp_r200m_values,rsp)
            r200m_values=np.append(r200m_values,r200m)
        print(logm_values)
        print(nu_values)
        # Plot for this omega_m
        #print(nu_values)
        #print(rsp_r200m_values)
        ax2 = plt.subplot(3,3,3)
        ax2.plot(nu_values, rsp_r200m_values, '-', label='$x=sp$')
        ax2.plot(nu_values, r200m_values, '-', label='$x=200m$')

        ax3 = plt.subplot(3,3,4)
        ax3.plot(nu_values, rsp_r200m_values/r200m_values, '-', label=f'$\Omega_m = {om}$')

        ax4 = plt.subplot(3,3,5)
        ax4.plot(logm_values, nu_values, '-', label=f'$\Omega_m = {om}$')

        ax5 = plt.subplot(3,3,6)
        ax5.plot(logm_values, rsp_r200m_values, '-', label='x=200m')


    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$r$')
    ax.set_ylabel(r'$\xi_{hm}$')


    ax1.set_xscale('log')
    ax1.set_xlabel('$r$')
    ax1.set_ylabel(r'$d\log \xi_{hm} / d \log r$')

    ax2.set_xlabel('$\\nu$')
    ax2.set_ylabel(r'$r_{x}$')

    ax3.set_xlabel('$\\nu$')
    ax3.set_ylabel(r'$r_{sp}/r_{200m}$')

    ax4.set_ylabel('$\\nu$')
    ax4.set_xlabel(r'$\log M_{200m}$')

    ax5.set_ylabel('$r_{sp}$')
    ax5.set_xlabel(r'$\log M_{200m}$')

    ax.legend(fontsize='xx-small')
    ax2.legend()   # for subplot (3,3,3)
    ax3.legend()   # for subplot (3,3,4)

    plt.tight_layout()
    plt.savefig('rsp_vs_nu_varying_omega.png', dpi=300)


if __name__ == "__main__":
    plot_rsp_vs_peak_height_varying_omega()
