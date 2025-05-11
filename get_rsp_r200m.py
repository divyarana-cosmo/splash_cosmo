from dark_emulator import darkemu
from scipy.interpolate import InterpolatedUnivariateSpline as ius
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.optimize import minimize
from scipy.signal import savgol_filter

emu = darkemu.base_class()

class constants():
      H0 = 100
      wv = 0.00064 # dark emulator has fixed neutrinos
      G = 4.301e-9
      delta_crit = 1.686 # critical density spherical collapse

class splashsim(constants):
    def __init__(self, wb=0.02225, wc=0.1198, omg_de=0.6844, ln1e10As=3.094, ns=0.9645, w=-1):
        self.wb         = wb
        self.wc         = wc
        self.omg_de     = omg_de
        self.ln1e10As   = ln1e10As
        self.ns         = ns
        self.w          = w

        cparam = np.array([wb, wc, omg_de, ln1e10As, ns, w])
        emu.set_cosmology(cparam)
        self.Om0    =   1-self.omg_de
        self.h      =   ((wb + wc + self.wv)/self.Om0)**0.5
        self.Ob0    =   wb/self.h**2
        self.sigma8 = emu.get_sigma8()
        self.init_spl_M200m2nu = False

    def M200m2r200m(self, logmh):
        m_tot = 10**logmh
        rho_crit_h = 3*(self.H0)**2/(8*np.pi*self.G)
        rho_mean_h = self.Om0 * rho_crit_h
        r_200m = (3*m_tot/(4*np.pi*200*rho_mean_h))**(1./3.)
        return r_200m

    def _spl_M200m2nu(self, z):
        mh, sigmh, mul_func = emu.get_f_HMF(z)
        nu = self.delta_crit/sigmh
        self.spl_logMh2nu = ius(np.log10(mh), nu)
        self.init_spl_M200m2nu = True
        self.init_redshift = z

    def logM200m2nu(self, logmh, z):
        if not self.init_spl_M200m2nu or z != self.init_redshift:
            self._spl_M200m2nu(z)
        return self.spl_logMh2nu(logmh)

    def logM200m2rsp(self, logmh, z):
        M = 10**logmh
        r200m = self.M200m2r200m(logmh)
        rs_fine = np.logspace(np.log10(0.02*r200m), np.log10(20*r200m), 80)
        xihm = emu.get_xicross_mass(rs_fine, M, z)

        from scipy.signal import savgol_filter
        slope = savgol_filter(np.log10(xihm),window_length=15,polyorder=3,deriv=1,delta=np.log10(rs_fine[1]/rs_fine[0]))
        from scipy.interpolate import InterpolatedUnivariateSpline as ius
        self.diff_func    =   ius(np.log10(rs_fine), slope)
        xx  = np.log10(rs_fine[1:-1])
        yy  = self.diff_func(xx)
        from scipy.optimize import minimize, rosen, rosen_der
        res = minimize(self.diff_func, x0=xx[np.argmin(yy)], bounds=((xx[0],xx[-1]),))
        rsp = 10**res.x[0]
        print(f"M={logmh:.2f} [h^-1 M_sun], r200m={r200m:.3f} [h^-1 Mpc], rsp={rsp:.3f} [h^-1 Mpc], rsp/r200m={rsp/r200m:.3f}, min_slope={res.fun:.3f}")
        return rsp

def plot_rsp_vs_peak_height_varying_omega():
    from colossus.cosmology import cosmology
    from colossus.halo import mass_so
    from colossus.lss import peaks
    from colossus.halo import splashback

    redshift = 0.0

    plt.figure(figsize=(9, 9))
    ax = plt.subplot(3, 3, 1)
    ax1 = plt.subplot(3, 3, 2)
    ax2 = plt.subplot(3, 3, 3)
    ax3 = plt.subplot(3, 3, 4)
    ax4 = plt.subplot(3, 3, 5)
    ax5 = plt.subplot(3, 3, 6)

    #for om in omega_m_values:
    #print(f"--- Processing Omega_m = {om} ---")
    #sim = splashsim(H0=H0, Om0=om, Ob0=omega_b)
    sim = splashsim()
    Ob_true = sim.Ob0

    params = {
        'flat': True,
        'H0': sim.H0,
        'Om0': sim.Om0,
        'Ob0': sim.Ob0,
        'sigma8': sim.sigma8,
        'ns': sim.ns
    }
    om = sim.Om0
    cosmology_name = f'cosmo_Om{om}'
    cosmology.addCosmology(cosmology_name, **params)
    cosmo = cosmology.setCosmology(cosmology_name)

    logm_values = np.linspace(13.5, 15.0, 5)
    nu_values = []
    rsp_values = []
    r200m_values = []

    for logm in logm_values:
        M = 10**logm
        rs = np.logspace(-1, 1, 200)
        xihm = emu.get_xicross_mass(rs, M, redshift)

        ax.plot(rs, xihm, label=f'$\\log M_{{h}} = {logm:.2f}$ $h^{{-1}}M_\\odot$')

        log_rs = np.log10(rs)
        log_xihm = np.log10(xihm)
        rsp = sim.logM200m2rsp(logm, redshift)
        yy = sim.diff_func(log_rs)
        xx =log_rs
        ax1.plot(10**xx, yy, '-')

        nu = sim.logM200m2nu(logm, redshift)
        nu_values.append(nu)
        rsp = sim.logM200m2rsp(logm, redshift)
        ax1.axvline(rsp, linestyle='--', alpha=0.5)
        rsp_values.append(rsp)

        r200m = sim.M200m2r200m(logm)
        r200m_values.append(r200m)

    nu_values = np.array(nu_values)
    rsp_values = np.array(rsp_values)
    r200m_values = np.array(r200m_values)

    colossus_masses = 10**np.array(logm_values)
    R200m_colossus = mass_so.M_to_R(colossus_masses, redshift, '200m')/1000.0

    ax2.plot(nu_values, rsp_values, '-o', label='$r_{sp}$')
    ax2.plot(nu_values, r200m_values, '-s', label='$r_{200m}$')
    ax2.plot(nu_values, R200m_colossus, '-^', label='colossus $r_{200m}$')

    ax3.plot(nu_values, rsp_values/r200m_values, '-o', label=f'$\\Omega_m = {om}$')
    ax3.plot(nu_values, splashback.modelMore15RspR200m(nu200m=nu_values, z=redshift, statistic='median'),
             '--', label='More+15')

    #nu_colossus = cosmo.sigma(r200m_values, z=redshift)
    nu_colossus = peaks.peakHeight(M=colossus_masses, z=redshift)
    ax4.plot(logm_values, nu_values, '-o', label='Our calc')
    ax4.plot(logm_values, nu_colossus, '--s', label='colossus')

    ax5.plot(logm_values, rsp_values, '-o', label='$r_{sp}$')
    ax5.plot(logm_values, r200m_values, '--s', label='$r_{200m}$')

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

    ax.legend(fontsize='x-small')
    ax2.legend(fontsize='small')
    ax3.legend(fontsize='small')
    ax4.legend(fontsize='small')
    ax5.legend(fontsize='small')

    plt.tight_layout()
    plt.savefig('rsp_vs_nu_varying_omega.png', dpi=300)

if __name__ == "__main__":
    plot_rsp_vs_peak_height_varying_omega()

