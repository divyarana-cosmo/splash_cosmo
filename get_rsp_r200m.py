from dark_emulator import darkemu
from scipy.interpolate import InterpolatedUnivariateSpline as ius
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.optimize import minimize
from scipy.signal import savgol_filter

emu = darkemu.base_class()

from model_fast import splash


class constants():
      H0 = 100
      wv = 0.00064 # dark emulator has fixed neutrinos
      G = 4.301e-9
      delta_crit = 1.686 # critical density spherical collapse

class splashsim(constants):
    def __init__(self, omg_de=0.7, w=-1, wb=0.02225, wc=0.1198, ln1e10As=3.094, ns=0.9645):
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
        M       = 10**logmh
        r200m   = self.M200m2r200m(logmh)
        #self.rs_fine = np.logspace(np.log10(0.1*r200m), np.log10(10*r200m), 300)
        self.rs_fine = np.logspace(-1, 1, 300)
        #rs_fine = np.logspace(np.log10(0.02), np.log10(20), 100)
        xihm    =   emu.get_xicross_mass(self.rs_fine, M, z)
        esd     =   emu.get_DeltaSigma_mass(self.rs_fine, M, z)
        nu      =   self.logM200m2nu(logmh, z)
        r200m   =   self.M200m2r200m(logmh)
        #golay filter thingy
        from scipy.signal import savgol_filter
        from scipy.optimize import minimize
        from scipy.interpolate import interp1d

        self.log_esd    =   np.log10(esd)
        self.log_xihm   =   np.log10(xihm)
        self.diff_xi3d  =   savgol_filter(np.log10(xihm), window_length=30, polyorder=3, deriv=1,delta=np.log10(self.rs_fine[1]/self.rs_fine[0])) #diff_xi3d
        self.logr       =   np.log10(self.rs_fine)
        func = interp1d(self.logr, self.diff_xi3d, kind='cubic')
        res  = minimize(func, x0=0.0, bounds=[(self.logr[0],self.logr[-1])])
        #xihm    =   10**(log_xihm)
        rsp     =   10**res.x[0]#rsp
        slope   =   res.fun#slope
        print(f"M={logmh:.2f} [h^-1 M_sun], r200m={r200m:.3f} [h^-1 Mpc], rsp={rsp:.3f} [h^-1 Mpc], rsp/r200m={rsp/r200m:.3f}, min_slope={slope:.3f}")
        return rsp

def plot_rsp_vs_peak_height():
    from colossus.cosmology import cosmology
    from colossus.halo import mass_so
    from colossus.lss import peaks
    from colossus.halo import splashback

    redshift = 0.0
    plt.figure(figsize=(8,8))
    ax  = plt.subplot(3, 3, 1)
    ax1 = plt.subplot(3, 3, 2)
    ax2 = plt.subplot(3, 3, 3)
    #ax3 = plt.subplot(2, 2, 4)

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

    logm_values = np.linspace(14.15, 15.0, 4)
    nu_values = []
    rsp_values = []
    r200m_values = []

    for nn,logm in enumerate(logm_values):
        M = 10**logm
        # peak height
        nu          = sim.logM200m2nu(logm, redshift)
        nu_values.append(nu)
        # splashback radius
        rsp         = sim.logM200m2rsp(logm, redshift)
        rsp_values.append(rsp)
        # overdensity size assignment#
        r200m = sim.M200m2r200m(logm)
        r200m_values.append(r200m)
        # plotting the delta sigma for each halo mass
        ax.plot(sim.rs_fine, 10**sim.log_esd, '-',label=f'$\log M_{{h}} = {logm:.2f}$', c='C%d'%nn, alpha=1.0)
        # plotting the xihm and the log-log slopes
        ax1.plot(sim.rs_fine, 10**sim.log_xihm, '-', c='C%d'%nn, alpha=1.0)
        # plotting the log-slope
        ax2.plot(10**sim.logr,sim.diff_xi3d, '-', c='C%d'%nn)
        ax2.axvline(rsp, linestyle='--', alpha=0.5, c='C%d'%nn)

    nu_values       = np.array(nu_values)
    rsp_values      = np.array(rsp_values)
    r200m_values    = np.array(r200m_values)

    #ax3.plot(nu_values, rsp_values/r200m_values, '-o', label=f'$\\Omega_m = {om:2.2f}$')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("$R\,[h^{-1}\mathrm{Mpc}]$")
    ax.set_ylabel("$\Delta\Sigma(R)\,[h M_\odot \mathrm{pc}^{-2}]$")
    ax.legend(fontsize='x-small')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('$r$ [$h^{-1}$Mpc]')
    ax1.set_ylabel(r'$\xi_{\rm hm}(r)$')

    ax2.set_xscale('log')
    ax2.set_xlabel('$r$ [$h^{-1}$Mpc]')
    ax2.set_ylabel(r'$d\log \xi_{\rm hm} / d \log r$')
    #ax3.set_xlabel(r'$\nu_{\rm 200m}$')
    #ax3.set_ylabel(r'$r_{\rm sp}/r_{\rm 200m}$')
    #ax3.legend(fontsize='small')
    plt.tight_layout()
    plt.savefig('rsp_vs_nu.png', dpi=300)
    plt.clf()
    return 0


def plot_rsp_vs_peak_height_varying_omega():
    redshift = 0.0
    #ax  = plt.subplot(2, 2, 1)
    #ax = plt.subplot(3, 3, 1)
    ax1 = plt.subplot(3, 3, 4)

    omega_m_values = np.array([0.20,0.25,0.30,0.35,0.40])
    for mm,om in enumerate(omega_m_values):
        sim     = splashsim(omg_de=1-om)
        rsp     = sim.logM200m2rsp(14.0, redshift)
        # plotting the delta sigma for each halo mass
        #ax.plot(sim.rs_fine, 10**sim.log_esd, '-', c='C%d'%mm, alpha=1.0)

        logm_values = np.linspace(14.15, 15.0, 20)
        nu_values       = []
        rsp_values      = []
        r200m_values    = []

        for nn,logm in enumerate(logm_values):
            # peak height
            nu          = sim.logM200m2nu(logm, redshift)
            nu_values.append(nu)
            # splashback radius
            rsp         = sim.logM200m2rsp(logm, redshift)
            rsp_values.append(rsp)
            # overdensity size assignment#
            r200m       = sim.M200m2r200m(logm)
            r200m_values.append(r200m)

        nu_values       = np.array(nu_values)
        rsp_values      = np.array(rsp_values)
        r200m_values    = np.array(r200m_values)
        idx = (nu_values>2.0) & (nu_values<3.0)
        #ax.plot(nu_values, rsp_values, '-',color='C%d'%mm)
        ax1.plot(nu_values[idx], rsp_values[idx]/r200m_values[idx], '-', label=f'$\\Omega_m = {om:2.2f}$',color='C%d'%mm)
        #ax1.plot(nu_values, rsp_values/r200m_values, '-', label=f'$\\Omega_m = {om:2.2f},\sigma_8={sim.sigma8:2.2f}$')

    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #ax.set_xlabel(r'$\nu_{\rm 200m}$')
    #ax.set_ylabel(r'$r_{\rm sp}$')
    ax1.set_ylim(1.05, 1.25)
    ax1.set_xlabel(r'$\nu_{\rm 200m}$')
    ax1.set_ylabel(r'$r_{\rm sp}/r_{\rm 200m}$')
    ax1.legend(fontsize='xx-small')
    #plt.tight_layout()
    plt.savefig('rsp_vs_nu_varying_omega.png', dpi=300)

#def plot_rsp_vs_peak_height_varying_w():
#    redshift = 0.0
#    #ax  = plt.subplot(2, 2, 1)
#    #ax = plt.subplot(3, 3, 1)
#    ax1 = plt.subplot(3, 3, 4)
#
#    w_values = np.array([-1.2, -1.1, -1.0, -0.9, -0.8])
#    for mm,ww in enumerate(w_values):
#        om      = 0.3
#        sim     = splashsim(omg_de=1-om, w=ww)
#        rsp     = sim.logM200m2rsp(14.0, redshift)
#        # plotting the delta sigma for each halo mass
#        #ax.plot(sim.rs_fine, 10**sim.log_esd, '-', c='C%d'%mm, alpha=1.0)
#
#        logm_values = np.linspace(14.15, 15.0, 20)
#        nu_values       = []
#        rsp_values      = []
#        r200m_values    = []
#
#        for nn,logm in enumerate(logm_values):
#            # peak height
#            nu          = sim.logM200m2nu(logm, redshift)
#            nu_values.append(nu)
#            # splashback radius
#            rsp         = sim.logM200m2rsp(logm, redshift)
#            rsp_values.append(rsp)
#            # overdensity size assignment#
#            r200m       = sim.M200m2r200m(logm)
#            r200m_values.append(r200m)
#
#        nu_values       = np.array(nu_values)
#        rsp_values      = np.array(rsp_values)
#        r200m_values    = np.array(r200m_values)
#
#        idx = (nu_values>2.0) & (nu_values<3.0)
#        #ax.plot(nu_values, rsp_values, '-',color='C%d'%mm)
#        ax1.plot(nu_values[idx], rsp_values[idx]/r200m_values[idx], '-', label=f'$w = {ww:2.2f}$',color='C%d'%mm)
#
#        #ax1.plot(nu_values, rsp_values/r200m_values, '-', label=f'$\\Omega_m = {om:2.2f},\sigma_8={sim.sigma8:2.2f}$')
#
#    #ax.set_xscale('log')
#    #ax.set_yscale('log')
#    #ax.set_xlabel(r'$\nu_{\rm 200m}$')
#    #ax.set_ylabel(r'$r_{\rm sp}$')
#
#    ax1.set_xlabel(r'$\nu_{\rm 200m}$')
#    ax1.set_ylabel(r'$r_{\rm sp}/r_{\rm 200m}$')
#    ax1.legend(fontsize='xx-small')
#    #plt.tight_layout()
#    plt.savefig('rsp_vs_nu_varying_w.png', dpi=300)




if __name__ == "__main__":
    #plot_rsp_vs_peak_height()
    plot_rsp_vs_peak_height_varying_omega()
    #plot_rsp_vs_peak_height_varying_w()
        #def model(x,Log_Rho_s, Log_R_s, Log_Rho_0, S_e, Log_R_t):
        #    alpha       = 0.155 + 0.0095*nu**2
        #    R_out       =  5*r200m
        #    Log_Alpha   =   np.log10(alpha)
        #    Log_Beta    =   np.log10(6)
        #    Log_Gamma   =   np.log10(4)
        #    F_cen       =   1.0
        #    R_off       =   0.1
        #    sp = splash(Log_Rho_s, Log_Alpha, Log_R_s, Log_Rho_0, S_e, Log_R_t, Log_Beta, Log_Gamma, F_cen, R_off, R_out)

        #    #val = 0.0*x
        #    #idx = (x>np.log10(0.5*r200m))

        #    #logxihm = np.log10(np.concatenate([
        #    #                    sp.rho_in(10**x[~idx]) * sp.f_trans(10**x[~idx]),
        #    #                    10**sp.log_xi_3d(x[idx])
        #    #                    ]))

        #    logxihm     =   sp.log_xi_3d(x)
        #    diff_xi_3d  =   sp.diff_xi_3d(x)
        #    ans         =   sp.rsp_3d(x)
        #    rsp         =   ans[0]
        #    slope       =   ans[1]
        #    return logxihm, diff_xi_3d, rsp, slope

        #def model_logxihm(x, Log_Rho_s, Log_R_s, Log_Rho_0, S_e, Log_R_t):
        #      logxihm, _, _, _ = model(np.log10(x), Log_Rho_s, Log_R_s, Log_Rho_0, S_e, Log_R_t)
        #      return 10**logxihm * x**2

        #from scipy.optimize import curve_fit
        #popt, pcov = curve_fit(model_logxihm, rs_fine, rs_fine**2 * xihm, p0=[4,0,0,2.0,0.0])
        #log_xihm, diff_xi3d, rsp, slope = model(np.log10(rs_fine),*popt)
