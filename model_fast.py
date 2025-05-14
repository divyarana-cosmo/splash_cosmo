# Given projected radial range will predict the 2pcf_2D and other functions
# Formalism is followed from ryoma's 2020 paper - https://arxiv.org/abs/2001.01160
import numpy as np
from scipy.integrate import quad,simps,cumtrapz,fixed_quad
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from scipy.misc import derivative
from scipy.optimize import fmin,minimize
import sys

class constants():
    """Useful constants"""
    G = 4.301e-9 #km^2 Mpc M_sun^-1 s^-2 gravitational constant
    H0 = 100. #h kms-1 Mpc-1 hubble constant at present
    omg_m = 0.315 #omega_matter
    not_so_tiny = 1e-24


class splash(constants):
    """Useful functions for predicting Rsp for both 2d and 3D case"""
    def __init__(self, Log_Rho_s, Log_Alpha, Log_R_s, Log_Rho_0, S_e, Log_R_t, Log_Beta, Log_Gamma, F_cen, R_off, R_out=1.5, R_max=40, splrmin = 0.001, splrbin=500):

        #parameters for density in inner regions
        self.rho_s = 10**Log_Rho_s
        self.alpha = 10**Log_Alpha
        self.r_s = 10**Log_R_s

        #parameters for density in outer region
        self.r_out = R_out
        self.rho_0 = 10**Log_Rho_0
        self.s_e = S_e

        #parameters for transition
        self.r_t = 10**Log_R_t
        self.beta = 10**Log_Beta
        self.gamma = 10**Log_Gamma

        #off-centering parameters
        self.fcen = F_cen
        self.roff = R_off

        #sigma computation variables
        self.rho_crt = 3*self.H0**2/(8*np.pi*self.G) # rho critical
        self.rho_m   = self.rho_crt*self.omg_m       # mean matter density in units h-1 Msun/(h-3Mpc^3)

        #R_max for project 3D-2pcf to 2D-2pcf
        self.R_max = R_max

        #putting spline on
        self.spl_rmin = splrmin
        self.spl_rbin = splrbin

        self.init_spl = False
        self.init_spl_cen = False
        self.init_spl_2d_diff = False
        self.init_spl_3d_diff = False

        # Add these precomputed values
        self.angular_weights = np.linspace(0, 2*np.pi, 100)
        self.integration_r = np.linspace(0, 10*self.roff, 400)
        self.p_off_values = self.p_off(self.integration_r)


    def rho_in(self,r):
        num = (-2.0/self.alpha)*((r*1.0/self.r_s)**self.alpha - 1)
        val = self.rho_s * np.exp(num)
        return val

    def rho_out(self,r):
        val = self.rho_0*(r*1.0/self.r_out)**(-self.s_e)
        return val

    def f_trans(self,r):
        val = (1 + (r*1.0/self.r_t)**self.beta)**(-self.gamma*1.0/self.beta)
        return val

    def log_xi_3d(self,log_r):
        r = 10**log_r
        val = self.rho_in(r)*self.f_trans(r) + self.rho_out(r)
        val = np.log10(val)
        return val

    def p_off(self,R):
        val = R*1.0/self.roff**2 * np.exp(-R**2/(2.0*(self.roff)**2))
        return val

    def log_xi_2d_off(self, log_R):
        """Optimized 2D correlation function for off-centered profile"""
        if not self.init_spl_cen:
            self.log_xi_2d_cen(log_R)

        # Convert to 1D array and get dimensions
        R = np.atleast_1d(10**log_R)
        r_values = np.linspace(0, 10*self.roff, 561)
        theta_values = np.linspace(0, 2*np.pi, 101)

        # Create reshaped views for broadcasting
        R_3d = R[:, np.newaxis, np.newaxis]          # Shape (N, 1, 1)
        r_3d = r_values[np.newaxis, :, np.newaxis]   # Shape (1, M, 1)
        theta_3d = theta_values[np.newaxis, np.newaxis, :]  # Shape (1, 1, K)

        # Calculate R_eff using vectorized operations
        R_eff = np.sqrt(R_3d**2 + r_3d**2 + 2*R_3d*r_3d*np.cos(theta_3d))

        # Vectorized spline evaluation
        log_R_eff = np.log10(R_eff)
        valid_mask = (log_R_eff >= self.log_xi_2d_dict['logr_min']) & \
                     (log_R_eff <= self.log_xi_2d_dict['logr_max'])

        xi_2d_values = np.full_like(R_eff, 10**self.log_xi_2d_dict['log_ximax'])
        xi_2d_values[valid_mask] = 10**self.spl_log_xi_2d_cen(log_R_eff[valid_mask])

        # Average over theta and integrate over r
        theta_avg = simps(xi_2d_values, theta_values,axis=2)/(2*np.pi)
        p_off = self.p_off(r_values)
        integral = simps(p_off * theta_avg, r_values, axis=1)
        return np.log10(integral.squeeze())



    def log_xi_2d_cen(self,log_R):
        if self.fcen == 1.0:
            rrad = np.logspace(log_R.min(), log_R.max(), self.spl_rbin)
            _val = 0.0*rrad
        else:
            rrad = np.logspace(np.log10(self.spl_rmin), np.log10(10**log_R.max() + 10*self.roff), self.spl_rbin)
            _val = 0.0*rrad
        for ii,rr in enumerate(rrad):
            _val[ii] = quad(lambda x: 10**(self.log_xi_3d(np.log10(np.sqrt(rr**2+x**2)))), 0, self.R_max)[0]
        val = _val*1.0/self.R_max

        self.spl_log_xi_2d_cen = InterpolatedUnivariateSpline(np.log10(rrad), np.log10(val),  k=4)
        self.spl_log_xi_2d_cen._data = {  # Cache evaluation points
        'x': np.log10(rrad),
        'y': np.log10(val)}
        self.spl_log_xi_2d_cen_diff = self.spl_log_xi_2d_cen.derivative()

        self.log_xi_2d_dict = {}
        self.log_xi_2d_dict['logr_min'] = np.log10(np.min(rrad))
        self.log_xi_2d_dict['logr_max'] = np.log10(np.max(rrad))
        self.log_xi_2d_dict['log_ximax'] = np.log10(np.max(val))
        self.init_spl_cen = True
        val = self.spl_log_xi_2d_cen(log_R)
        return val

    def log_xi_2d(self,log_R):
        if self.fcen==1.0:
            val = self.log_xi_2d_cen(log_R)
            self.spl_log_xi_2d      = self.spl_log_xi_2d_cen
            self.spl_log_xi_2d_diff = self.spl_log_xi_2d_cen_diff
        else:
            val = self.fcen * 10**self.log_xi_2d_cen(log_R) + (1 - self.fcen) * 10**self.log_xi_2d_off(log_R)

            self.spl_log_xi_2d      = InterpolatedUnivariateSpline(log_R, np.log10(val),  k=4)
            self.spl_log_xi_2d_diff = self.spl_log_xi_2d.derivative()
            self.init_spl = True
            val = np.log10(val)
        return val

    def log_sigma(self, log_R):
        "current sigma implementation is for the scalar log_R using shin et al 2021 units"
        rhoeff = 2*self.R_max *self.rho_m/1e12 +  2*self.R_max*10**self.log_xi_2d(log_R)
        if 10**log_R<self.spl_rmin:
            val = 2*self.R_max *self.rho_m/1e12 +  2*self.R_max*10**self.log_xi_2d(np.log10(self.spl_rmin))
            return np.log10(val)
        return np.log10(rhoeff)



    def diff_xi_2d(self,log_R):
        if not self.init_spl:
            dirt = self.log_xi_2d(log_R)
        val = self.spl_log_xi_2d_diff(log_R)
        return val

    def diff_xi_3d(self,log_R):
        R = 10**log_R
        xi_3d = 10**self.log_xi_3d(log_R)
        d_rhoin_dr = self.rho_in(R)*(-2)*((R/self.r_s)**(self.alpha -1))/self.r_s
        d_ftran_dr = -(self.gamma/self.beta)*((1+(R/self.r_t)**self.beta)**(-self.gamma/self.beta -1))*(self.beta*(R/self.r_t)**(self.beta -1))/self.r_t
        d_rhoout_dr = -(self.rho_0 * self.s_e) * ((R/self.r_out)**(-self.s_e -1))/self.r_out
        val = self.rho_in(R) * d_ftran_dr + self.f_trans(R)*d_rhoin_dr + d_rhoout_dr
        val = val * R/xi_3d
        return val

    def rsp_2d(self,log_R):
        if not self.init_spl:
            dirt = self.log_xi_2d(log_R)
        val = self.spl_log_xi_2d_diff(log_R)
        idx = np.argmin(val)
        gval = log_R[idx]
        bnds = ((log_R.min(), log_R.max()),)
        val  = minimize(self.spl_log_xi_2d_diff, gval, bounds=bnds)
        return np.append(10**val.x,val.fun)

    def rsp_3d(self,log_R):
        val = self.diff_xi_3d(log_R)
        idx = np.argmin(val)
        gval = log_R[idx]
        bnds = ((log_R.min(), log_R.max()),)
        val  = minimize(self.diff_xi_3d, gval, bounds=bnds)
        return np.append(10**val.x,val.fun)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    plt.figure(figsize=[10.0,10.0])
    plt.subplot(2,2,1)

#log_rho_s, log_alpha, log_r_s, log_rho_0, s_e, log_r_t, log_beta, log_gamma
    #2.05628643, -0.74687378, -0.32809567, -3.32641305,  5.71911196, 0.13970109,  0.79832924,  0.70256978

    #parameters for density in inner regions
    log_rho_s = 1.47611335 #1.10#np.log10(27.08)#np.log10(27.08)#1.97
    log_alpha = 0.68742786 #-0.95#np.log10(0.16) #-0.68#np.log10(0.16)#-0.68
    log_r_s = 0.67286521 #-0.32#np.log10(0.32) #-0.58#np.log10(0.32)#-0.58
    #parameters for density in outer region
    log_rho_0 = 0.31423027 #0.349#np.log10(0.43) #-0.48#np.log10(0.43)#-0.48
    s_e = 2.07541219 #1.601 #1.3#1.61#1.3

    #parameters for transition
    log_r_t   = 0.04709024 #-0.082#np.log10(1.31) #0.16#np.log10(1.31)#0.16
    log_beta  = 0.89337944 #0.762#np.log10(3.71) #0.77#np.log10(3.71)#0.77
    log_gamma = 0.13532044 #0.66#np.log10(6.42) #0.64#np.log10(6.42)#0.64

    import sys
    #off-centering parameters
    f_cen = 1.0
    r_off = 0.33717233

    ss = splash(Log_Rho_s = log_rho_s, Log_Alpha = log_alpha, Log_R_s = log_r_s, Log_Rho_0 = log_rho_0, S_e = s_e, Log_R_t = log_r_t, Log_Beta = log_beta, Log_Gamma = log_gamma, F_cen = f_cen, R_off = r_off, splrmin=0.001, splrbin=100)

    #rr = 0.001018207976387967
    #rtest = np.linspace(0,50,10)
    #vtest = 0.0*rtest
    #for i in range(len(rtest)):
    #    vtest[i] =  10**(ss.log_xi_3d(np.log10(np.sqrt(rr**2+rtest[i]**2))))
    #print vtest

    #print quad(lambda x: 10**(ss.log_xi_3d(np.log10(np.sqrt(rr**2+x**2)))), 0, ss.R_max)[0]

    #print ss.log_xi_2d(np.log10(rr))

    rad = np.logspace(np.log10(0.05),np.log10(20), 9)

    logr = np.log10(rad)
    print((ss.log_xi_2d_cen(logr)))
    print((ss.log_xi_2d_off(logr)))
    print((ss.log_xi_2d(logr)))
    print((ss.diff_xi_2d(logr)))
    print((ss.diff_xi_3d(logr)))
    #plt.plot(10**logr, ss.diff_xi_2d(logr))
    #plt.xscale('log')
    #plt.savefig('test.png', dpi=300)
    print((rad.min(), rad.max()))
    #print ss.log_xi_2d(logr)
    #plt.plot(10**logr, 10**ss.log_xi_2d(logr))
    ####plt.plot(10**logr, np.log10(ss.rho_in(10**logr)))
    ####plt.plot(10**logr, np.log10(ss.f_trans(10**logr)))
    ####plt.plot(10**logr, np.log10(ss.rho_out(10**logr)))
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.savefig('test.png', dpi=300)
    ##print ss._int_xi_2d_off(10**logr[0])
    print('Rsp3d = %s, val = %s'%(*ss.rsp_3d(logr),))
    print('Rsp2d = %s, val = %s'%(*ss.rsp_2d(logr),))




