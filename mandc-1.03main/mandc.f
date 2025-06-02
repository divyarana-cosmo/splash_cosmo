c=============================================================================
c  Compiling:
c
c     You can make the executable file by typing
c
c        > make            or by typing       > f77 -O -o mandc.x mandc.f
c
c     Other tested compilers are g77, f90, f95 and ifort, so you can choose
c     a compiler on your computer with proper flags either in above command
c     or in the file 'makefile'. 
c
c
c=============================================================================
c  Running:
c
c     You may run the executable file by typing
c
c        > ./mandc.x
c     then you are prompted to input parameters interactively,
c
c     or by typing    > ./mandc.x < inputfilename
c     when you have saved input parameters in a file inputfilename.
c
c     Input parameters and output results are expalined below in detail.
c
c
c=============================================================================
c  Plotting results:
c
c     You may find some sample Super Mongo macros in file 'mandc.sm', which
c     can be used to plot both 1) mass & redshift dependence of halo median
c     concentration c(z,M) and 2) evolution of halo properties along its
c     mass accretion history, such as M(z|zobs,Mobs) and c(z|zobs,Mobs).
c
c     Some example postscript figures are provided in example/ directory for
c     you to check results you obtain against them.
c


c=============================================================================
c  This code was written by
c  Donghai Zhao
c  Copyright (C) 2007
c
c     Permission is granted to anyone to make or distribute verbatim
c     copies of this document provided that the copyright notice and
c     this permission notice are preserved, and that the distributor
c     grants the recipient permission for further redistribution as
c     permitted by this notice.
c
c  Version 1.00 (2008.11.02):
c     Only a main program version was released at that time.
c  Version 1.03 (2010.09.14):
c     1) We add a new output, the age of the universe, to the end of each row
c     in the output file;
c     2) For clarity, in the case of power spectrum type 1, we make use of a
c     new fitting formula to calculate the transfer function, which was provided
c     also by Eisentein & Hu (1998) but does NOT take into account effects of
c     massive neutrinos on the transfer function.
c     Actually, this modification will induce only TINY changes in the final
c     results for universe without massive neutrino, like our concordant LCDM
c     universe.
c     A main program version and two subroutine versions are released at the
c     same time. This is the MAIN PROGRAM version and on the same web site you
c     can also find the subroutine versions which perform the same as this main
c     program version but are much easier for you to call.
c
c
c=============================================================================
c  Reporting bug:
c
c     This code is maintained by Donghai Zhao. If you find any bug or if
c     you have any suggestion or question, please write to
c     dhzhao@shao.ac.cn
c
c
c=============================================================================
c  Getting the latest version:
c
c     The latest version (with bug fixes) can be obtained from the web site:
c     http://www.shao.ac.cn/dhzhao/mandc.html
c
c     There a calculator which allows one to interactively generate data for
c     any given cosmological model and a paper with an appendix describing
c     the step-by-step implementation of our models are also provided.
c
c
c=============================================================================
c  Credits:
c
c     Any publication which is benefited from the current code is supposed 
c     to contain reference to the paper Zhao, Jing, Mo & Boerner 2009, ApJ,
c     707, 354 / arXiv: 0811.0828 and acknowledgement of use of the code.
c     In addition, we will appreciate preprints of publications based on it,
c     and please send them to one of the following addresses:
c     dhzhao@shao.ac.cn
c     ypjing@shao.ac.cn
c     hjmo@astro.umass.edu
c


c=============================================================================
c  According to Zhao, Jing, Mo & Boerner (2009), this program predicts
c     1) the median mass accretion history (the median main branch of merger
c        trees) for halos of any given mass Mobs at any given redshift zobs,
c        M(z|zobs,Mobs), and
c     2) evolution of halo structural properties along this median main branch
c        (such as concentraion evolution c(z|zobs,Mobs), see below for details)
c
c     You may choose to output
c     a) time evolution of the mass and structural properties along the median
c        main branch of final halos with a given mass Mobs (such as above
c        M(z|zobs,Mobs) and c(z|zobs,Mobs)), OR
c     b) mass dependence of these structural properties at final redshift zobs
c        (such as that of median concentration, c(zobs,M); and you can run
c        repeatedly with many other zobs to obtain the redshift dependence).
c
c     Mass accretion history and concentration evolution history with different
c     halo definitions are also available.
c     
c     This version can be used for Einstein de Sitter, open, or flat
c     universe with cosmological costant, as some fitting formulae from
c     the literature (such as formulae for halo virial density and for
c     linear growth factor) are limited to these cosmologies.
c
c-----------------------------------------------------------------------------
c
c  1.Input format:
c
c     You may save input parameters in an input file beforehand.
c
c     Here is a sample input file that may be cut out and used
c     (for the so-called concordance LCDM model; evolution along MAH).
cx-------------------------------------------------------------------
c        LC0001                     run (any 6-character string)
c        0.3  0.7                   omegam0,omegaL0       
c        2                          ispec
c        0.666667                   h
c        0.90                       sigma8
c        1.0                        tilt
c        1                          ioutput
c        0.0                        zobs
c        12.0                       lgMobs in [M_sun/h]
cx-------------------------------------------------------------------
c
c     Here is another sample input file that may be cut out and used
c     (for the so-called concordance LCDM model; mass dependence).
cx-------------------------------------------------------------------
c        LC0001                     run (any 6-character string)
c        0.3  0.7                   omegam0,omegaL0       
c        2                          ispec
c        0.666667                   h
c        0.90                       sigma8
c        1.0                        tilt
c        2                          ioutput
c        0.0                        zobs
cx-------------------------------------------------------------------
c
c     Here is another sample input file that may be cut out and used
c     (for WMAP5 cosmology).
cx-------------------------------------------------------------------
c        WMAP05                     run (any 6-character string)
c        0.258  0.742               omegam0,omegaL0
c        1                          ispec
c        0.72                       h
c        0.796                      sigma8
c        0.963                      tilt
c        0.0438       2.726         omegab0,T_cmb (from version 1.03)     
c        2                          ioutput
c        0.0                        zobs
cx-------------------------------------------------------------------
c
c     Here is another sample input file that may be cut out and used
c     (for WMAP3 cosmology).
cx-------------------------------------------------------------------
c        WMAP03                     run (any 6-character string)
c        0.238  0.762               omegam0,omegaL0
c        1                          ispec
c        0.73                       h
c        0.75                       sigma8
c        0.95                       tilt
c        0.042        2.726         omegab0,T_cmb (from version 1.03)     
c        2                          ioutput
c        0.0                        zobs
cx-------------------------------------------------------------------
c
c     Here is another sample input file that may be cut out and used
c     (for WMAP1 cosmology).
cx-------------------------------------------------------------------
c        WMAP01                     run (any 6-character string)
c        0.268  0.732               omegam0,omegaL0
c        1                          ispec
c        0.71                       h
c        0.90                       sigma8
c        1.0                        tilt
c        0.044        2.726         omegab0,T_cmb (from version 1.03)     
c        2                          ioutput
c        0.0                        zobs
cx-------------------------------------------------------------------
c
c     Here is another sample input file that may be cut out and used
c     (for you to supply your own transfer function).
cx-------------------------------------------------------------------
c        Mine06                     run (any 6-character string)
c        0.3  0.7                   omegam0,omegaL0       
c        -1                         ispec
c        mytransferfilename         myfile
c        0.90                       sigma8
c        1.0                        tilt
c        2                          ioutput
c        0.0                        zobs
cx-------------------------------------------------------------------
c
c     Here is another sample input file that may be cut out and used
c     (for you to supply your own power spectrum).
cx-------------------------------------------------------------------
c        Mine08                     run (any 6-character string)
c        0.3  0.7                   omegam0,omegaL0       
c        -2                         ispec
c        myspectrumfilename         myfile
c        0.90                       sigma8     
c        2                          ioutput
c        0.0                        zobs
cx-------------------------------------------------------------------
c
c     Here is another sample input file that may be cut out and used
c     (for power-law power spectrum).
cx-------------------------------------------------------------------
c        PL0005                     run (any 6-character string)
c        1.000  0.000               omegam0,omegaL0
c        0                          ispec
c        -0.5                       tilt
c        2                          ioutput
c        0.0                        zobs
cx-------------------------------------------------------------------
c
c     Here is another sample input file that may be cut out and used
c     (for power-law power spectrum).
cx-------------------------------------------------------------------
c        PL0005                     run (any 6-character string)
c        1.000  0.000               omegam0,omegaL0
c        0                          ispec
c        -0.5                       tilt
c        1                          ioutput
c        0.0                        zobs
c        0.5                        lgMobs in [M_star0]
cx-------------------------------------------------------------------
c
c-----------------------------------------------------------------------------
c
c  2.Input:
c
c        run      --  any character string of length 6 to name the current run
c                     and to be embedded in the output file name
c
c     cosmological parameters:
c        omegam0  --  matter density at redshift 0 in critical units; including
c                     dark matter and baryon; 0<omegam0<=1 
c        omegaL0  --  cosmological constant at redshift 0 in critical units;
c                     omegaL0=1-omegam0 or omegaL0=0
c        ispec    --  power spectrum type indicator
c                     (-2) for you to supply your own power spectrum
c                     (-1) for you to supply your own transfer function
c                     (0)  power-law power spectrum
c                     (1)  Eisenstein & Hu 1998 power spectrum
c                     (2)  BBKS 1986 power spectrum
c                     (3)  Bond & Efstathiou 1987 power spectrum
c                     (4)  DEFW 1985 power spectrum
c        myfile   --  name of the file in which your own power spectrum or
c                     transfer function is saved; needed when ispec=-2/-1
c        h        --  Hubble constant at redshift 0 in [100 km/s/Mpc];
c                     needed when ispec>0
c        sigma8   --  rms of the linear density field smoothed in a sphere of
c                     radius 8 Mpc/h at redshift 0; not needed when ispec=0
c        tilt     --  (when ispec=-2) needn't be input
c                     (when ispec=0) power index of linear power spectrum;
c                     -3 < tilt < 1
c                     (for other cases) power index of primodial power spectrum;
c                     -3 < tilt
c
c     following 4 are NOT needed unless ispec=1 (Eisentein & Hu spectrum), and
c     from the release of version 1.03, omeganu0 and N_nu of the 4 parameters are
c     NOT needed any longer even for ispec=1, as we adopt a new fitting formula
c     to calculate the transfer function, which was provided also by Eisentein &
c     Hu (1998) but does not consider effects of massive neutrinos on the transfer
c     function: 
c        omegab0  --  baryon density at redshift 0 in critical units
c        omeganu0 --  massive neutrino density at redshift 0 in critical units;
c                     0<=omegab0+omeganu0<=omegam0
c        T_cmb    --  CMB temperature at redshift 0
c        N_nu     --  number of degenerate massive neutrinos
c
c     following 3 are used to specify halo population to be studied
c        ioutput  --  output type indicator
c                     (1) properties along the median mass accretion history of
c                     final offspring halos of mass lgMobs WILL be output
c                     (2) properties of final offspring halos of many different
c                     masses at redshift zobs WILL be output
c        zobs     --  redshift at which you choose final offspring halos
c        lgMobs   --  logarithmic VIRIAL mass of halos you choose; need when
c                     ioutput=1
c                     (when ispec=0) in [M_star0], i.e., in characteristic 
c                     nonlinear mass at redshift 0
c                     (for other cases) in [M_sun/h]
c
c     If you choose to supply your own power spectrum or transfer funtion by
c     setting ispec=-2/-1, column 1 should be log10(k/[h/Mpc]), and column 2
c     should also be logarithmic value. Its normalization is arbitary.
c
c-----------------------------------------------------------------------------
c
c  3.Output format:
c
c     output file name: 
c
c        when input inoutput=1, properties along the median mass accretion
c        history of final offspring halos of mass lgMobs are output to file
c        'mchistory_run.intzstring.intlgMstring', with
c        intzstring=int(zobs*100) and intlgMstring=int(lgMobs*100)
c
c        when input ioutput=2, properties of final offspring halos of many
c        different masses at redshift zobs are output to file
c        'mc_run.intzstring'
c
c     both the files have the same format:
c
c        row 1:  ispec,nbin1,zobs,omegam0,omegaL0,sigma8,tilt,omegab0,
c                omeganu0,T_cmb,N_nu,myfile/h/note ('---' will take place
c                of those parameters that needn't be input)
c
c        row >1: ziz,Miz,c,rhohiz,R,V,Ms,rhos,Rs,Vs,Miz_200c,c_200c,
c                rhohiz_200c,Miz_200m,c_200m,rhohiz_200m,uniage_iz
c
c-----------------------------------------------------------------------------
c
c  4.Output:
c
c        nbin1    --  final offspring number
c                     (when input ioutput=1) nbin1=1, output a main branch
c                     (when input ioutput=2) nbin1>1, output offsprings 
c        note     --  a note to remind reader that mass, radius and velocity
c                     etc. have been scaled when ispec=0
c        ziz      --  redshift for output
c                     (when input ioutput=1) ziz>= zobs
c                     (when input ioutput=2) ziz = zobs
c        uniage_iz--  universe age at ziz, here in [year/h] 
c
c     following are for the spherical virial halo definition 
c        Miz      --  halo virial mass, in [M_sun/h] (except for the special
c                     case with ispec=0, see below)
c                     (when input ioutput=1) for the final offspring halo of
c                     logarithmic mass lgMobs and its main progenitors;
c                     Miz=M(ziz|zobs,Mobs)
c                     (when input ioutput=2) for final offspring halos only;
c                     lgMobs is renewedly uniformly choosed in (7,16) or so;
c                     Miz=M(zobs|zobs,Mobs)=Mobs
c        c        --  halo virial concentration
c                     (when input ioutput=1) c=c(ziz|zobs,Mobs)
c                     (when input ioutput=2) c=c(zobs|zobs,Mobs)=c(zobs,Mobs)
c        rhohiz   --  corresponding halo virial density at redshift ziz, in 
c                     [h^2 M_sun/Mpc^3]
c        R        --  halo virial radius, in [Mpc/h]
c        V        --  halo circular velocity at virial radius, in [km/s]
c
c     following are some characteristic inner quntities
c        Rs       --  halo characteristic inner radius Rs, in [Mpc/h]
c        Ms       --  mass enclosed in Rs,                 in [M_sun/h]
c        rhos     --  density at Rs,                 in [h^2 M_sun/Mpc^3]
c        Vs       --  circular velocity at Rs,             in [km/s]
c
c     following are for some other halo definitions 
c      rhohiz_200c--  200 times critical density at redshift ziz, in
c                     [h^2 M_sun/Mpc^3]
c        Miz_200c --  corresponding halo mass with halo density equal to 
c                     rhohiz_200c, in [M_sun/h]
c        c_200c   --  corresponding halo concentration with halo density
c                     equal to rhohiz_200c
c      rhohiz_200m--  200 times universe average matter density at redshift
c                     ziz, in [h^2 M_sun/Mpc^3]
c        Miz_200m --  corresponding halo mass with halo density equal to 
c                     rhohiz_200m, in [M_sun/h]
c        c_200m   --  corresponding halo concentration with halo density
c                     equal to rhohiz_200m
c
c     For the case of power-law power spectrum (ispec=0), the age of the
c     universe and all halo properties, except concentration, have been
c     scaled respectively with universe age at redshift 0 and corresponding
c     quantities of halos of mass M_star0 at redshift 0. For example, in
c     this special case, universe age, halo mass and radius are in [uniage0],
c     [M_star0] and [R_star0] respectively. 
c
c     For main progenitors that are the earliest (also the smallest)
c     in the output file, if Miz is less than 25 times of the smallest
c     main progenitor, negative c, Ms, rhos, Rs, Vs, Miz_200c, c_200c,
c     Miz_200m and c_200m will be output (because the characteristic time
c     t_0.04 has to be computed, see Zhao et al. 2009).
c
c     For users interested in different halo definitions: For the second
c     output type, mass dependences of halo properties in DIFFERENT HALO
c     DEFINITIONs are all output as shown above. For the first output type,
c     the halo mass input, lgMobs, should always be VIRIAL mass, even though
c     masses in other halo definitions are output again as shown above. If
c     you want to output property evolution for final offspring halos of
c     a given mass in ANOTHER HALO DEFINITION, i.e., if you want to start
c     from a mass in another halo definition, you may do it in two steps:
c     1. by setting ioutput=2, you can obtain one-to-one correspondence 
c     among halo masses in different definitions at the redshift in
c     consideration for the current cosmology; 2. by setting ioutput=1
c     and specifying the corresponding virial halo mass, you will get the
c     history you want.
c


      Program mandc

      implicit none
      integer*4 NMAX,NMAXT,nm1,nm2
      integer*4 NMSTP,NMSTPT,nz,iz,izmax,nbin,nbin1,ibin
      integer*4 ispec,ioutput,iargc,intzobs,intlgMobs,j,lblnk
      parameter(NMAX=5001,NMSTP=400)
      real*8 omegam0,omegaL0,h,tilt,sigma8
      real*8 omegab0,omeganu0,T_cmb,N_nu
      real*8 omhh,f_baryon,uniage0,uniage1d
      real*8 lgamass2
      real*8 lggrmmin,lggrmmax,n_lggrm,dlggrm
      real*8 zobs,lgMobs,dlga,dlgm,dlgat,dlgmt
      real*8 lgmmin1,lgmmax1,dm1,lgsig1,dydx1
      real*8 lgmmin2,lgmmax2,dm2,lgsig2,dydx2
      real*8 lgz1min,lgz1max,dlgz1,lgd,zt,dczt,lgz1t
      real*8 z,lgdcz,lgm,lgsigma1,lgs
c      parameter(n_lggrm=10.,dlga=0.01,dlgm=dlog10(0.5d0)/8.)
      parameter(n_lggrm=10.,dlga=0.01,dlgm=-0.037628749)
      real*8 lgz1(0:NMSTP),c,Ms,Rs,Vs,rhos
      real*8 lgMiz,Miz,R,V,ziz,rhohiz,uniage_iz,dlgmtp,ztp,uniage_tp
      real*8 ez2iz,ez2,omegamiz,omegaz
      real*8 rhohiz_200c,c_200c,rMiz_200c,Miz_200c
      real*8 rhohiz_200m,c_200m,rMiz_200m,Miz_200m
      real*8 M8,dcz0,rhoh0,M_star0,R_star0,V_star0
      character inf1*50,myfile*50
      character run*8,intzstring*5,intlgMstring*5
      character*17 sigma8char,tiltchar
      character*17 omegab0char,omeganu0char,T_cmbchar,N_nuchar
      common /lgd/ lgz1min,lgz1max,dlgz1,lgd(NMAX),NMAXT
      common /lgsig1/ lgmmin1,lgmmax1,dm1,lgsig1(NMAX),nm1
      common /lgsig2/ lgmmin2,lgmmax2,dm2,lgsig2(NMAX),nm2
      common /mhist/ dlgat,dlgmt,z(0:NMSTP),lgdcz(0:NMSTP),
     &   lgm(0:NMSTP),lgsigma1(0:NMSTP),lgs(0:NMSTP),nz,NMSTPT


      write(*,*)' '
      write(*,*)'Please give the current run a 6-character name.'
      read(*,*)run

      write(*,*)' '
      write(*,*)'omegam0,   omegaL0  ?'
      read(*,*)omegam0,omegaL0
      if(.NOT.(omegam0.gt.0.and.omegam0.le.1.and.(omegaL0.eq.0
     &   .or.omegaL0+omegam0.eq.1)))then
         write(*,*)'Sorry.'
         write(*,*)
     &   'This version can only be used for Einstein de Sitter, '
         write(*,*)'   open, or flat universe with cosmology costant.'
         stop 'Pls check omegam0 and omegaL0 you input.'
      endif

      write(*,*)' '
      write(*,*)'power spectrum type indicator  ?'
      write(*,*)'   (-2) for you to supply your own power spectrum'
      write(*,*)'   (-1) for you to supply your own transfer function'
      write(*,*)'   (0)  power-law power spectrum'
      write(*,*)'   (1)  Eisenstein & Hu 1998 power spectrum'
      write(*,*)'   (2)  BBKS 1986 power spectrum'
      write(*,*)'   (3)  Bond & Efstathiou 1987 power spectrum'
      write(*,*)'   (4)  DEFW 1985 power spectrum'
      read(*,*)ispec
      if(ispec.lt.-2.or.ispec.gt.4)stop 'ispec out of range'

      write(*,*)' '
      if(ispec.eq.-2)then
         write(*,*)'Name of your own powerspectrum file?'
         read(*,*)myfile
      elseif(ispec.eq.-1)then
         write(*,*)'Name of your own transferfunction file?'
         read(*,*)myfile
      elseif(ispec.eq.0)then
         h=1.    ! no use
         myfile='Mass, radius and velocity etc. have been scaled'
      else
         write(*,*)'h  ?'
         read(*,*)h
         if(h.gt.10.or.h.le.0)stop 'H_0/[100km/s/Mpc] needed'
         write(myfile,'(e17.8)')h
      endif

      if(ispec.eq.0)then
         sigma8=1.0d0
         sigma8char='  ---------------'
      else
         write(*,*)' '
         write(*,*)'sigma8,  ?'
         read(*,*)sigma8
         if(sigma8.le.0)stop 'sigma8 less than/equal to 0'
         write(sigma8char,'(e17.8)')sigma8
      endif

      if(ispec.ge.-1)then
         write(*,*)' '
         if(ispec.eq.0)then
            write(*,*)'tilt (power index of linear power spectrum)?'
            read(*,*)tilt
            if(tilt.le.-3.or.tilt.ge.1)stop 'tilt out of range'
         else
            write(*,*)'tilt (power index of primodial power spectrum)?'
            read(*,*)tilt
            if(tilt.le.-3)stop 'tilt out of range'
         endif
         write(tiltchar,'(e17.8)')tilt
      else
         tilt=0.  ! no use
         tiltchar='  ---------------'
      endif

      if(ispec.eq.0)then
         M8=2.7754e+11*omegam0*3.1415926535*4./3.*8.0d0**3
         call deltacrit(0.0d0,omegam0,omegaL0,dcz0)
         M_star0=(sigma8/dcz0)**(6/(tilt+3))*M8
         call halodensv(omegam0,omegaL0,0.0d0,rhoh0)
         call rvh(M_star0,rhoh0,R_star0,V_star0)
         uniage0=uniage1d(omegam0,omegaL0,100.0d0,0.0d0) ! in [yr/h]
      else
         M_star0=1.0d0
         R_star0=1.0d0
         V_star0=1.0d0
         rhoh0=1.0d0
         uniage0=1.0d0
      endif

      if(ispec.eq.1)then
         write(*,*)' '
c         write(*,*)'omegab0,   omeganu0,   T_cmb,   N_nu  ?'  ! version 1.00
c         read(*,*)omegab0,omeganu0,T_cmb,N_nu
c         if(omegab0.lt.0)stop 'omegab0 less than 0'
c         if(omeganu0.lt.0)stop 'omeganu0 less than 0'
c         if(omegab0+omeganu0.gt.omegam0)stop
c     &      'omegab0+omeganu0 larger than omegam0'
c         call TFmdm_hu(omegam0,omegaL0,omegab0,omeganu0,h,T_cmb,N_nu)
c         write(omegab0char,'(e17.8)')omegab0
c         write(omeganu0char,'(e17.8)')omeganu0
c         write(T_cmbchar,'(e17.8)')T_cmb
c         write(N_nuchar,'(e17.8)')N_nu
         write(*,*)'omegab0,   T_cmb  ?'       ! from version 1.03
         read(*,*)omegab0,T_cmb
         if(omegab0.lt.0)stop 'omegab0 less than 0'
         if(omegab0.gt.omegam0)stop 'omegab0 larger than omegam0'
         omhh=omegam0*h*h
         f_baryon=omegab0/omegam0
         call TFset_parameters(omhh, f_baryon, T_cmb)
         write(omegab0char,'(e17.8)')omegab0
         write(T_cmbchar,'(e17.8)')T_cmb
         omeganu0char='  ---------------'
         N_nuchar=omeganu0char
      else
         omegab0char='  ---------------'
         omeganu0char=omegab0char
         T_cmbchar=omegab0char
         N_nuchar=omegab0char
      endif

      write(*,*)' '
      write(*,*)
     & 'If ioutput=1, properties along the median accretion history of '
      write(*,*)'   final offspring halos of given mass WILL be output.'
      write(*,*)
     & 'If input ioutput=2, properties of final offspring halos of many'
      write(*,*)'   different masses at redshift zobs WILL be output.'
      write(*,*)'ioutput  ?'
      read(*,*)ioutput
      if(.NOT.(ioutput.eq.1.or.ioutput.eq.2))stop 'ioutput out of range'

      write(*,*)' '
      write(*,*)'zobs  ?'
      read(*,*)zobs
      if(zobs.le.-1)stop 'zobs less than/equal to -1'

      if(ioutput.eq.1)then
         write(*,*)' '
         if(ispec.eq.0)then
            write(*,*)
     &        'lgMobs (in characteristic nonlinear mass at redshift 0)?'
         else
            write(*,*)'lgMobs (in [M_sun/h])?'
         endif
         read(*,*)lgMobs
         lgMobs=lgMobs+dlog10(M_star0)
         if(lgMobs.lt.-3.8.or.lgMobs.gt.23.8)stop
     &      'lgMobs out of range'
      endif


      NMAXT=NMAX
      NMSTPT=NMSTP
      dlgat=dlga
      dlgmt=dlgm
      lgamass2=1.0d0
      lggrmmin=7.
      lggrmmax=16.0
      if(ispec.eq.0.and.tilt.ge.-0.5)lggrmmin=9.
      print*,'preparing tables for given cosmology and power spectrum'
      print*,'...'

c preparing sigma(M) table, extrapolated to redshift 0
c With ispec=-3, you can supply your own sigma(M) table in a file named
c 'sigma_'//run in the working directory, using the same format as in
c subroutine s_m().
c You can also supply your own power spectrum table or transfer function
c table in 'myfile' with ispec=-2/-1 respectively.
      if(ispec.ne.-3)call s_m(omegam0,h,tilt,sigma8,ispec,run,myfile)
      inf1='sigma_'//run
      call init_read(inf1,lgmmin1,lgmmax1,nm1,dm1,NMAX,lgsig1)

c preparing s(M) table
      call s_m_modify(run,lgamass2)
      inf1='s_'//run
      call init_read(inf1,lgmmin2,lgmmax2,nm2,dm2,NMAX,lgsig2)

c preparing delta_c(z) table
      lgz1min=0.
      lgz1max=dlga*NMSTP+1.
      dlgz1=(lgz1max-lgz1min)/(NMAX-1)
      do j=1,NMAX
         lgz1t=lgz1min+dlgz1*(j-1)
         zt=10.**lgz1t-1.
         call deltacrit(zt,omegam0,omegaL0,dczt)
         lgd(j)=dlog10(dczt)  ! delta_c(z)
      enddo

c preparing average density table as a function of r/Rs in NFW halo,
c (ln(1+x)-x/(1+x))/x**3
      call densinr

c preparing output file
      intzobs=zobs*100
      write(intzstring,'(i5.5)')intzobs
      if(ioutput.eq.1)then
         intlgMobs=(lgMobs-dlog10(M_star0))*100
         write(intlgMstring,'(i5.5)')intlgMobs
         inf1='mchistory_'//run//'.'//intzstring//'.'//intlgMstring
         nbin=0
      else
         inf1='mc_'//run//'.'//intzstring
         nbin=(lggrmmax-lggrmmin)*n_lggrm+1.0001
         dlggrm=1./n_lggrm
         lggrmmax=lggrmmin+dlggrm*nbin
      endif
      write(*,*)' '
      write(*,*)'Results will be output to file: ',inf1
      write(*,*)' '
      write(*,*)'Sample Super Mongo mcaros for ploting the results'
      write(*,*)'   can be found in file: mandc.sm.'
      write(*,*)' '
      nbin1=nbin+1

      open(88,file=inf1,status='unknown')
c writing input parameters to the file
      write(88,801)ispec,nbin1,zobs,omegam0,omegaL0,sigma8char,tiltchar,
     &   omegab0char,omeganu0char,T_cmbchar,N_nuchar,myfile

c predicting median halo mass accretion history
      do ibin=0,nbin
         if(nbin.gt.0)lgMobs=lggrmmin+ibin*dlggrm
         call median_mh(omegam0,omegaL0,zobs,lgMobs)
         do iz=0,nz
            lgz1(iz)=dlog10(1.+z(iz))
         enddo
         if(nbin.gt.0)then
            izmax=0
         else
            izmax=nz
         endif
         do iz=0,izmax
c predicting median concentration for a given POINT on the history
            call mah2c1(nz,lgz1,lgm,omegam0,omegaL0,iz,c,
     &         lgMiz,dlgmtp,ziz,ztp,uniage_iz,uniage_tp)
c calculating halo global properties
            Miz=10.**lgMiz
            call halodensv(omegam0,omegaL0,ziz,rhohiz)
            call rvh(Miz,rhohiz,R,V)
            if(c.gt.0.2.and.c.lt.800000.)then
c calculating characteristic halo inner properties
               call rvs(Miz,rhohiz,c,Ms,Rs,rhos,Vs)
c converting to other halo definitions
               ez2iz=ez2(omegam0,omegaL0,ziz)
               omegamiz=omegaz(omegam0,omegaL0,ziz)
               rhohiz_200c= 2.7754d+11*ez2iz*200
               rhohiz_200m= 2.7754d+11*ez2iz*omegamiz*200
               call cconvert(rhohiz,rhohiz_200c,c,c_200c,rMiz_200c)
               call cconvert(rhohiz,rhohiz_200m,c,c_200m,rMiz_200m)
               Miz_200c=Miz*rMiz_200c
               Miz_200m=Miz*rMiz_200m
            else
               Ms=-2.0
               Rs=-2.0
               rhos=-2.0
               Vs=-2.0
               c_200c=-1.5
               Miz_200c=-1.
               c_200m=-1.5
               Miz_200m=-1.
            endif
c writing results to the file
            write(88,802)ziz,
     &         Miz/M_star0,c,rhohiz/rhoh0,R/R_star0,V/V_star0,
     &         Ms/M_star0,rhos/rhoh0,Rs/R_star0,Vs/V_star0,
     &         Miz_200c/M_star0,c_200c,rhohiz_200c/rhoh0,
     &         Miz_200m/M_star0,c_200m,rhohiz_200m/rhoh0,
     &         uniage_iz/uniage0
         enddo
      enddo

      close(88)

 801  format(2i10,3e17.8,6a,'   ',a)
 802  format(17e17.8)

      end


c=============================================================================
c This subrotine predicts the median mass accretion history (median main 
c branch of merger trees) for halos of given mass at given redshift in
c a given cosmology, according to Zhao, Jing, Mo & Boerner (2009).
c
c Input:
c
c     omega0    --  matter density at redshift 0 in critical units 
c     lambda0   --  cosmological constant at redshift 0 in critical units
c     zobs      --  redshift at which you choose the final offspring halo 
c     lgmobs    --  logarithmic mass of the halo you choose
c     lgd()     --  log10(delta_c(z)) table with log10(1+z) uniformly spaced 
c     NMAXT     --  number of lgd()
c     lgsig1()  --  log10(sigma(M)) table with log10(M) uniformly spaced
c                     between lgmmin1 and lgmmax1
c     nm1       --  number of lgsig1()
c     lgsig2()  --  log10(s(M)) table with log10(M) uniformly spaced
c                     between lgmmin2 and lgmmax2
c     nm2       --  number of lgsig2()
c     dlga      --  logarithmic interval of expansion factor for output
c     dlgm      --  logarithmic interval of halo mass for output (need 
c                    only in some extreme case)
c
c Output:
c
c     nz        --  number of progenitors predicted
c     z()       --  redshifts of the final offspring halo and predicted
c                     progenitors in the output history, with number nz+1
c     lgdcz()   --  logarithmic delta_c() corresponding to z()
c     lgm()     --  logarithmic mass of halos in the history
c     lgsigma1()--  rms of linear density field smoothed at scale lgm()
c     lgs()     --  corresponding logarithmic s as defined in the paper
c
c

      subroutine median_mh(omega0,lambda0,zobs,lgmobs)
      implicit none

      integer*4 NMAX,NMAXT,NMSTP,NMSTPT
      parameter(NMAX=5001,NMSTP=400) 
      integer*4 iargc,nm1,nm2,imstp,i,j
      integer*4 nz,iz
      real*8 omega0,lambda0,yacc
      real*8 lgmmin1,lgmmax1,dm1,lgsig1,dydx1
      real*8 lgmmin2,lgmmax2,dm2,lgsig2,dydx2
      real*8 lgz1min,lgz1max,dlgz1,lgd
      real*8 z,lgdcz,lgm,lgsigma1,lgs
      real*8 dlgsig_dcz,sigma_growth,dczt,s,w,p,pobs,wp,lgmt,zt,lgz1t
      real*8 zobs,lgmobs,finlgdcz     !,dnu2sqr,nusqrshft_end
      real*8 p_obs,p_z,fin_lgdcz
      real*8 dlga,dlgm         !,mstarzobs
c      real*8 ziz,ztp,c,c1,dlgmtp,uniage_iz,uniage_tp,denshiz,denshtp
c      real*8 lgmiz,ztp2,c2,dlgmtp2,uniage_tp2,denshtp2
c      real*8 dcziz,dcztp,dcztp2
c      real*8 brmass(0:NMSTP),lgz1(0:NMSTP)
c      real*8 lgsigma254,zf,v
c      real*8 omegaiz,omegaz,denshiz_200c,denshiz_200m,ez2,ez2iz
c      real*8 omegatp,denshtp_200c,denshtp_200m,ez2tp
c      real*8 omegatp2,denshtp2_200c,denshtp2_200m,ez2tp2
c      real*8 c2_200c,c2_200m,c1_200c,c1_200m,c_200c,c_200m
c      real*8 rmass2_200c,rmass2_200m,rmass1_200c,rmass1_200m
c      real*8 rmass_200c,rmass_200m
c      real*8 rtp,rhogrowth
      common /lgd/ lgz1min,lgz1max,dlgz1,lgd(NMAX),NMAXT
      common /lgsig1/ lgmmin1,lgmmax1,dm1,lgsig1(NMAX),nm1
      common /lgsig2/ lgmmin2,lgmmax2,dm2,lgsig2(NMAX),nm2
      common /mhist/ dlga,dlgm,z(0:NMSTP),lgdcz(0:NMSTP),
     &   lgm(0:NMSTP),lgsigma1(0:NMSTP),lgs(0:NMSTP),nz,NMSTPT

      if(NMAXT.ne.NMAX)stop 'NMAXT.ne.NMAX in median_mh'
      if(NMSTPT.ne.NMSTP)stop 'NMSTPT.ne.NMSTP in median_mh'


      z(0)=zobs
      lgm(0)=lgmobs

c to prepare delta_c(zobs),sigma(Mobs),s(Mobs)
      if(lgm(0).le.lgmmin1+dm1.or.lgm(0).ge.lgmmax1-dm1)goto 902
      call deltacrit(z(0),omega0,lambda0,dczt)
      lgdcz(0)=dlog10(dczt)  ! delta_c(z)
      call ftable(lgm(0),lgsigma1(0),dydx1,lgmmin1,lgmmax1,nm1,dm1,
     &   lgsig1)
      lgmt=lgm(0)
      if(lgmt.le.lgmmin2+dm2.or.lgmt.ge.lgmmax2-dm2)goto 902
      call ftable(lgmt,lgs(0),dydx2,lgmmin2,lgmmax2,nm2,dm2,
     &   lgsig2)

c to prepare w(zobs,Mobs), p(zobs,zobs,Mobs) and final nonlinear time scale
      w=10.d0**(lgdcz(0)-lgs(0))
      pobs=p_obs(w)
      finlgdcz=fin_lgdcz(w)

c to trace the median MAH backword iterately
      yacc=0.00001d0
      do imstp=1,NMSTP

c to get growth rate of sigma, dlg(sigma(M))/dlg(delta_c(z)) 
         w=10.d0**(lgdcz(imstp-1)-lgs(imstp-1))
         p=p_z(lgdcz(imstp-1),lgdcz(0),finlgdcz,pobs)
         wp=w-p
         dlgsig_dcz=sigma_growth(wp)

c to get z and delta_c(z) of this new step
         z(imstp)=10.**(dlog10(z(0)+1)+imstp*dlga)-1.
         call deltacrit(z(imstp),omega0,lambda0,dczt)
         lgdcz(imstp)=dlog10(dczt)  ! delta_c(z) 

c to get sigma(m) of this step 
         lgsigma1(imstp)=lgsigma1(imstp-1)+dlgsig_dcz*
     &      (lgdcz(imstp)-lgdcz(imstp-1))

c in another direction (one may skip this part at first reading)
         if(lgsigma1(imstp).le.lgsig1(nm1-1)*1.05.or.lgsigma1(imstp)
     &      .ge.lgsig1(2)/1.05)then
            do i=imstp,NMSTP
               lgm(i)=lgm(i-1)+dlgm
               if(lgm(i).le.lgmmin1+dm1.or.lgm(i).ge.lgmmax1-dm1)then
                  nz=i-1
                  goto 902
               endif
               call ftable(lgm(i),lgsigma1(i),dydx1,lgmmin1,lgmmax1,
     &            nm1,dm1,lgsig1)
               w=10.0d0**(lgdcz(i-1)-lgs(i-1))
               p=p_z(lgdcz(i-1),lgdcz(0),finlgdcz,pobs)
               wp=w-p
               dlgsig_dcz=sigma_growth(wp)
               lgdcz(i)=lgdcz(i-1)+(lgsigma1(i)-lgsigma1(i-1))
     &            /dlgsig_dcz
               lgz1t=dlog10(z(i-1)+1.)    !guess
               call findx(lgz1t,lgdcz(i),yacc,lgz1min,lgz1max,
     &            NMAX,dlgz1,lgd)
               z(i)=10.0**lgz1t-1
               
               lgmt=lgm(i)
               if(lgmt.le.lgmmin2+dm2.or.lgmt.ge.lgmmax2-dm2)then
                  nz=i-1
                  goto 902
               endif
               call ftable(lgmt,lgs(i),dydx2,lgmmin2,lgmmax2,
     &            nm2,dm2,lgsig2)
            enddo
            nz=i-1
            goto 902
         endif

c to get mass of this step with sigma(M)
         lgm(imstp)=lgm(imstp-1)    !guess
         call findx(lgm(imstp),lgsigma1(imstp),yacc,lgmmin1,lgmmax1,
     &      nm1,dm1,lgsig1)

c to get s(M) of this step for further iteration
         lgmt=lgm(imstp)
         if(lgmt.le.lgmmin2+dm2.or.lgmt.ge.lgmmax2-dm2)then
            nz=imstp-1
            goto 902
         endif
         call ftable(lgmt,lgs(imstp),dydx2,lgmmin2,lgmmax2,
     &      nm2,dm2,lgsig2)

      enddo
      nz=imstp-1

 902  end


c-------------------------------------------------------------
      function sigma_growth(wp)
c by Donghai Zhao
      implicit none
      real*8 wp,sigma_growth
      sigma_growth=wp/5.85 
      end
 
c-------------------------------------------------------------
      function p_obs(wobs)
c by Donghai Zhao
      implicit none
      real*8 wobs,wtp,p_obs
      wtp=4.0 
      p_obs=wobs/2/(wobs**6+wtp**6)*wtp**6
      end

c-------------------------------------------------------------
      function fin_lgdcz(wobs)
c by Donghai Zhao
      implicit none
      real*8 wobs,fin_lgdcz
      fin_lgdcz=0.272d0/wobs
      end
c-------------------------------------------------------------
      function p_z(lgdcz,lgdczobs,finlgdcz,pobs)
c by Donghai Zhao
      implicit none
      real*8 lgdcz,lgdczobs,dlgdcz,finlgdcz,pobs,p_z
      dlgdcz=lgdcz-lgdczobs
      if(dlgdcz.lt.finlgdcz)then
         p_z=(finlgdcz-dlgdcz)/finlgdcz*pobs
      else
         p_z=0.d0
      endif
      end

c-------------------------------------------------------------
c This subroutine changes sigma(M) to s(M) according to 
c Zhao, Jing, Mo & Boerner (2009)

      subroutine s_m_modify(model,lgamass)
      implicit none
      integer*4 i,nsample,nsample1,NMAX
      parameter(NMAX=10000)
      real*8 xmin,xmax,dx,lgamass
      real*8 yy(NMAX),xx(NMAX)
      character model*(*),inf1*50

      inf1='sigma_'//model
      open(29,file=inf1,status="old")
      read(29,*)xmin,xmax,nsample,dx
      nsample1=nsample+1
      if(nsample1.gt.NMAX)stop 'nsample1.gt.NMAX in s_m_modify'
      do i=1,nsample1
         read(29,*)xx(i),yy(i)
      enddo
      close(29)

      inf1='s_'//model
      open(89,file=inf1,status="unknown")
      write(89,701)xmin+dx,xmax-dx,nsample-2,dx ! 080709
      do i=2,nsample1-1
         write(89,702)xx(i),max(-20.,yy(i)+(yy(i+1)-
     &      yy(i-1))/2/dx*lgamass)
      enddo
      close(89)
 701  format(2(e17.8),i8,e17.8)
 702  format(2(e17.8))

      end



c=============================================================================
c This subroutine predict concentration for a halo at ANY POINT (i.e., either
c for the final offspring or for any progenitor) of the mass accretion history
c you provid according to Zhao, Jing, Mo & Boerner (2009) if its mass is larger
c than 25 times that of the earlist main progenitor in the history, otherwise
c return negative value. 
c
c In order to specify the POINT, you can choose to input either mass, 
c redshift or snapshot number of the POINT.
c
c As for the input mass accretion history, neither lgz1() nor lgm() is 
c required to be regularly spaced.
c
c Here virial halo defitnion is allowed only. If you need concentration
c in other definitions, please change the mass and concentration with
c subroutines 'densinr' and 'cconvert'. 
c
c     nz        --  snapshot number of progenitores in the mass accretion 
c                     history provided
c     lgz1()    --  log10(z+1) for the final offspring and progenitors in 
c                     the history, here z representing redshift
c     lgm()     --  logarithmic halo mass along the history
c     omega0    --  matter density at redshift 0 in critical units
c     lambda0   --  cosmological constant at redshift 0 in critical units
c     iz        --  indicator of the POINT on the history
c                   (<= -3) do not run this subroutine
c                   (= -2) for you to input halo mass of the POINT, lgmiz
c                   (=-1) for you to input redshift of the POINT, ziz
c                   (>= 0) snapshot number of the POINT
c     lgmiz     --  logarithmic halo mass of the POINT
c     ziz       --  redshift of the POINT
c     uniage_iz --  universe age at ziz, here in [yr/h] 
c     dlgmtp    --  log10(0.04), model parameter in Zhao et al. (2009)
c     ztp       --  redshift when progenitor mass is 4% of that at ziz
c     uniage_tp --  universe age at ztp, here in [yr/h] 
c     c         --  predicted concentration for offspring/progenitor at
c                     redshift ziz, i.e., for the POINT you specified
c

      subroutine mah2c1(nz,lgz1,lgm,omega0,lambda0,iz,c,
     &   lgmiz,dlgmtp,ziz,ztp,uniage_iz,uniage_tp)
      implicit none
      integer*4 NZMAX,nz,nz1,iz,i,imdiffmin
      parameter(NZMAX=10000)
      real*8 lgz1(0:nz),lgm(0:nz),omega0,lambda0,lgz1iz,lgmiz
      real*8 dlgmtp,lgmtp,lgz1tp,ztp,t1t,dt1,t1(0:NZMAX)
      real*8 mdiff,mdiffmin,uniage_iz,uniage_tp,ziz,c,yacc,dydx,m2c1
      real*8 uniage1d

      dlgmtp=dlog10(4.0d-2) 

      if(nz.lt.0)stop 'nz less than 0 in mah2c1'
      if(nz.gt.NZMAX)stop 'nz larger than NZMAX in mah2c1'
      dt1=1.
      do i=0,nz
         t1(i)=dt1*i
      enddo
      nz1=nz+1
      yacc=0.00001d0

c You can choose to specify either mass, redshift or snapshot number.

      if(iz.le.-3)then           ! do not run mah2c1
         print*,'iz.le.-3, you have choosen not to run mah2c1.' 
         goto 105

      elseif(iz.eq.-2)then       ! given lgmiz
         if(lgmiz.gt.lgm(0).or.lgmiz.lt.lgm(nz))
     &      stop 'Sorry, the mass you specified is out of range.'
         imdiffmin=nz/2
         mdiffmin=1.0e+30
         do i=0,nz
            mdiff=dabs(lgm(i)-lgmiz)
            if(mdiff.lt.mdiffmin)then
               mdiffmin=mdiff
               imdiffmin=i
            endif
         enddo
         t1t=t1(imdiffmin)                     ! guess
         call findx(t1t,lgmiz,yacc,t1(0),t1(nz),nz1,dt1,lgm)
         call ftable(t1t,lgz1iz,dydx,t1(0),t1(nz),nz1,dt1,lgz1)
         ziz=10.**lgz1iz-1.d0

      elseif(iz.eq.-1)then       ! given redshift ziz
         lgz1iz=dlog10(ziz+1.d0)
         if(lgz1iz.lt.lgz1(0).or.lgz1iz.gt.lgz1(nz))
     &      stop 'Sorry, the redshift you specified is out of range.'
         imdiffmin=nz/2
         mdiffmin=1.0e+30
         do i=0,nz
            mdiff=dabs(lgz1(i)-lgz1iz)
            if(mdiff.lt.mdiffmin)then
               mdiffmin=mdiff
               imdiffmin=i
            endif
         enddo
         t1t=t1(imdiffmin)                     ! guess
         call findx(t1t,lgz1iz,yacc,t1(0),t1(nz),nz1,dt1,lgz1)
         call ftable(t1t,lgmiz,dydx,t1(0),t1(nz),nz1,dt1,lgm)

      else                       ! given snapshot number iz
         if(iz.gt.nz)
     &      stop 'Sorry, iz you specified is out of range.'
         ziz=10.**lgz1(iz)-1.d0
         lgmiz=lgm(iz)
      endif

      lgmtp=lgmiz+dlgmtp
      if(lgmtp.le.lgm(0).and.lgmtp.ge.lgm(nz))then

c redshift for progenitor mass to be given fraction of that at given POINT

         imdiffmin=nz/2
         mdiffmin=1.0e+30
         do i=0,nz
            mdiff=dabs(lgm(i)-lgmtp)
            if(mdiff.lt.mdiffmin)then
               mdiffmin=mdiff
               imdiffmin=i
            endif
         enddo
         t1t=t1(imdiffmin)                     ! guess

         yacc=0.00001d0
         call findx(t1t,lgmtp,yacc,t1(0),t1(nz),nz1,dt1,lgm)
         call ftable(t1t,lgz1tp,dydx,t1(0),t1(nz),nz1,dt1,lgz1)
         ztp=10.**lgz1tp-1.

c universe age at these two redshifts

         uniage_iz=uniage1d(omega0,lambda0,100.0d0,ziz) ! in [yr/h]
         uniage_tp=uniage1d(omega0,lambda0,100.0d0,ztp)
         c=m2c1(uniage_iz,uniage_tp)

      else
         uniage_iz=uniage1d(omega0,lambda0,100.0d0,ziz)
         c=-0.5
         if(nz.eq.0)c=-0.6
         ztp=-2.
         uniage_tp=-1.
      endif

 105  return
      end


c---------------------------------------------------------
      function m2c1(t0,t1)
      implicit none
      real*8 t0,t1,m2c1
      m2c1=(4.**8+(t0/t1)**8.4)**(1./8)
      end





c################ SUBROUTINTE GROUP1: on space ##########################

c========================================================================

        subroutine s_m(omega0,h0,tilt0,sigma_8,ispec0,model,myfile)
c-----------------------------------------------------------------------
c     Generates sigma(M) look-up table - D.H. Zhao May 2003.
c     using routines of A. Jenkins.
c-----------------------------------------------------------------------

        implicit none

        real*8 omega,omega1,h,tilt,sigma_8,omega0,h0,tilt0
        real*8 rsphere,const,unnsigma
        real*8 dummy
        real*8 logmin,logmax,dm,logm,sig
        integer*4 nsample,im,ispec,ispec0
        character model*(*),myfile*(*),sigfile*30

        common /cosmo/ omega
        common /shape/ omega1,h,tilt,ispec
        common /rad/ rsphere
        common /norm/ const

        omega=omega0
        omega1=omega0
        h=h0
        tilt=tilt0
        ispec=ispec0

        logmin=-4.   ! minimum mass in lg h^{-1}M_sun
        logmax=24.   ! maximum mass in lg h^{-1}M_sun
        nsample=2800 ! bin number

c  Use your own power spectrum or transfer function.
        if(ispec.eq.-2.or.ispec.eq.-1)call read_function(myfile)   

c-------------Find normalisation constant-------------------------
        rsphere = 8.0d0
        const = sigma_8**2/unnsigma(dummy)
c-----------------------------------------------------------------

        dm=(logmax-logmin)/nsample

        sigfile='sigma_'//model
        open(25,file=sigfile,status='unknown')
        write(25,701)logmin,logmax,nsample,dm
        do im=0,nsample
           logm = logmin+im*dm
           call sm(logm,sig)
           if(sig.lt.10.**(-20.))sig=10.**(-20.) !!!!
           write(25,702)logm,dlog10(sig)
        enddo
        close(25)
 701    format(2(e17.8),i8,e17.8)
 702    format(2(e17.8))

        end

c-----------------------------------------------------------------
        subroutine sm(logm,sig)
        implicit none
        real*8 rho_crit
        real*8 logm,rm,rmass,sig
        real*8 omega

        common /cosmo/ omega
        parameter (rho_crit=2.7754e11)


        rm = 10**logm
        rmass = rm/rho_crit/omega

        call sigma_m(rmass,sig)

        end

c---------------------------------------------------------------
        subroutine  sigma_m(m,sig)
        implicit none
c
c   Use unit of mass where 1h^{-1}Mpc^3 has mass 1

        real*8 m,sig,rsphere,const,unnsigma,dummy
        common /rad/ rsphere
        common /norm/ const

        rsphere = (3.*m/4./3.1415926535)**0.33333333333

        sig = sqrt(const*unnsigma(dummy))

        return
        end
c--------------------------------------------------------------
c modified by Donghai Zhao

        function unnsigma(dummy)
        implicit none
        real*8 dummy
        real*8 dxk,sum2,xk,evar2,unnsigma
        parameter (dxk = 0.01)
        real*8 omega,h,tilt,rsphere
        integer*4 ispec
        common /shape/ omega,h,tilt,ispec
        common /rad/ rsphere

        if(ispec.eq.0)then
           unnsigma=rsphere**(-tilt-3.)
           return
        endif

        sum2 = 0.0
        do xk= -20.,20.0,dxk
           sum2 = sum2 + evar2(xk)*dxk
        enddo
        unnsigma = sum2
        end


        real*8 function evar2(x)
        implicit none
        real*8 rk,x,var2
        rk = exp(x)
        evar2 = var2(rk)*rk
        return
        end


        real*8 function var2(x)
        implicit none
        real*8 x,powspec,weight
        var2 = weight(x)*powspec(x)
        return
        end


c   Power Spectrum

        real*8 function powspec(xk)
        implicit none
        integer*4 ncontinue,ispec
        real*8 xk,xklow,xkup
        real*8 omega,h,tilt,q,func_eval
        real*8 TFtransfer
        common /shape/ omega,h,tilt,ispec


      if(ispec.eq.-2)then 
c--Select one's own power spectrum----------------------------------------
         powspec = func_eval(xk)/xk**3
c-------------------------------------------------------------------------

      elseif(ispec.eq.-1)then 
c--Select one's own transfer function-------------------------------------
         powspec = xk**tilt*func_eval(xk)**2
c-------------------------------------------------------------------------

      elseif(ispec.eq.0)then
c--Power-law power spectrum-----------------------------------------------
         powspec = xk**tilt
c-------------------------------------------------------------------------

      elseif(ispec.eq.1)then
c--Eisenstein & Hu 1998---------------------------------------------------
c         T_master = TF_master(xk*h)
c         powspec=xk**tilt*T_master**2
         powspec=xk**tilt*TFtransfer(xk*h)**2

      elseif(ispec.eq.2)then
c--BBKS 1986--------------------------------------------------------------
         q = xk/(omega*h)
         powspec=xk**tilt*((dlog(1.d0+2.34*q)/2.34/q)/
     &   (1.d0+3.89*q+(16.1*q)**2.+(5.46*q)**3.+(6.71*q)**4)**0.25)**2
c-------------------------------------------------------------------------

      elseif(ispec.eq.3)then
c--Bond and Efstathiou 1987-----------------------------------------------
         q = xk/(omega*h)
         powspec = xk**tilt/(1.+ (6.4*q + (3.0*q)**1.5 +
     &      (1.7*q)**2)**1.13)**(2./1.13)
c-------------------------------------------------------------------------

      elseif(ispec.eq.4)then
c--DEFW 1985--------------------------------------------------------------
         q = xk/(omega*h)
         powspec=xk**tilt/(1.+1.7*q+9.0*q**1.5+1.0*q**2)**2

      endif


        return
        end


        real*8 function weight(x) ! Top-hat filter function x 4pi k^2
        implicit none
        real*8 x,y,rsphere
        common /rad/ rsphere

        y = rsphere*x        
        weight = 3.1415926535*x*x*36.*(sin(y)/y/y/y-cos(y)/y/y)**2

        return
        end


      subroutine read_function(myfile)
      implicit none
      integer*4 npoints,i,np
      parameter (npoints=6000)
      real*8 rkvec(npoints), rpow(npoints)
      character myfile*(*)
      common /inter/ rkvec,rpow,np


      open (10,file=myfile,status='old')
      do i=1,npoints
       read (10,*,end=10) rkvec(i),rpow(i)
      enddo
      stop 'File too long!'
 10   continue
      close (10)
      np = i-1

      return
      end


      real*8 function func_eval(tval)
      implicit none
      integer*4 npoints,np,i
      real*8 val,tval,s,pv
      parameter (npoints=6000)
      real*8 rkvec(npoints), rpow(npoints)
      common /inter/ rkvec,rpow,np

      val = dlog10(tval)

      if ((val.le.rkvec(1)).or.(val.gt.rkvec(np))) then
       print*,tval,val,rkvec(1),rkvec(np),np
       stop 'val out of range'
      endif
      do i=1,np-1     ! Bisection method would be great improvement.

        if ((val.gt.rkvec(i)).and.(val.le.rkvec(i+1))) then
            s = (val - rkvec(i))/(rkvec(i+1)-rkvec(i))
            pv = (1.-s)*rpow(i) + s * rpow(i+1)
            func_eval = 10**pv
            return
        endif
      enddo

      stop 'Should not reach this point!'
      end


c====Eisenstein & Hu 1998 subroutines====================================

c There are two routines and a set of functions. 
c   TFset_parameters() sets all the scalar parameters, while
c   TFtransfer_function() calculates various transfer functions
c
c Global variables -- We've left many of the intermediate results as
c global variables in case you wish to access them, e.g. by declaring
c them as a common block in your main program.
c
c Note that all internal scales are in Mpc, without any Hubble constants!
c

        subroutine TFset_parameters(omhh0,f_baryon0,Tcmb0)

        implicit none
        real*8 omhh0,f_baryon0,Tcmb0
        real y,omhh,obhh,Tcmb
        real theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality,
     &       sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,
     &       f_baryon,beta_node
        common/GLOBALVARIABLES/theta_cmb,z_equality,k_equality,z_drag,
     &       R_drag,R_equality,sound_horizon,k_silk,alpha_c,beta_c,
     &       alpha_b,beta_b,beta_node
        common/GLOBALVAR2/omhh,f_baryon,Tcmb

c Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula */
c Input omhh -- The density of CDM and baryons, in units of critical dens,
c                multiplied by the square of the Hubble constant, in units
c                of 100 km/s/Mpc */
c       f_baryon -- The fraction of baryons to CDM */
c       Tcmb -- The temperature of the CMB in Kelvin, 2.728(4) is COBE and is
c               the default reached by inputing Tcmb=0 -- reset on output. */
c Output nothing, but set many global variables in common block
c       GLOBALVARIABLES. You can access them yourself, if you want:
c
c       theta_cmb,      /* Tcmb in units of 2.7 K */
c       z_equality,     /* Redshift of matter-radiation equality, really 1+z */
c       k_equality,     /* Scale of equality, in Mpc^-1 */
c       z_drag,         /* Redshift of drag epoch */
c       R_drag,         /* Photon-baryon ratio at drag epoch */
c       R_equality,     /* Photon-baryon ratio at equality epoch */
c       sound_horizon,  /* Sound horizon at drag epoch, in Mpc */
c       k_silk,         /* Silk damping scale, in Mpc^-1 */
c       alpha_c,        /* CDM suppression */
c       beta_c,         /* CDM log shift */
c       alpha_b,        /* Baryon suppression */
c       beta_b,         /* Baryon envelope shift */



        omhh=omhh0
        f_baryon=f_baryon0
        Tcmb=Tcmb0

c Are inputs reasonable?

        if (f_baryon.le.0) f_baryon=1.e-5
        if (Tcmb.le.0) Tcmb=2.726   !!!
c        if (omhh.le.0.0) then
c           write(6,*) 'TFset_parameters(): Illegal input' 
c           pause
c        end if
c
c        if (hubble.gt.10.0) then
c           write(6,*) 'TFset_parameters(): WARNING, Hubble constant in
c     &                 100km/s/Mpc desired'
c        end if

c Auxiliary variables
        obhh = omhh*f_baryon
        theta_cmb = Tcmb/2.7

c Main variables
        z_equality = 2.50e4*omhh*theta_cmb**(-4.) - 1.D0
        k_equality = 0.0746*omhh*theta_cmb**(-2.)

          z_drag = 0.313*omhh**(-0.419)*(1.+0.607*omhh**(0.674))
          z_drag = 1e0 + z_drag*obhh**(0.238*omhh**(0.223))
        z_drag = 1291e0 * omhh**(0.251)/
     &           (1e0 + 0.659*omhh**(0.828)) * z_drag

        R_drag = 31.5*obhh*theta_cmb**(-4.)*1000e0/(1e0 + z_drag)
        R_equality = 31.5*obhh*theta_cmb**(-4.)
     &               *1000e0/(1e0 + z_equality)

        sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*
     &      log(( sqrt(1.+R_drag)+sqrt(R_drag+R_equality) )
     &       /(1.+sqrt(R_equality)))

        k_silk = 1.6*obhh**(0.52)*omhh**(0.73)*
     &           (1e0 + (10.4*omhh)**(-0.95))

          alpha_c = ((46.9*omhh)**(0.670)*(1e0+(32.1*omhh)**(-0.532)))
          alpha_c = alpha_c**(-f_baryon)
        alpha_c = alpha_c*((12.0*omhh)**(0.424)*(1e0 +
     &             (45.0*omhh)**(-0.582)))**(-f_baryon**3.)


          beta_c = 0.944/(1+(458.*omhh)**(-0.708))
          beta_c = 1.+beta_c*((1.-f_baryon)**((0.395*omhh)**(-0.0266))
     &          - 1e0)
          beta_c = 1./beta_c

          y = (1e0+z_equality)/(1e0+z_drag)
          alpha_b = y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)
     &              /(sqrt(1.+y)-1.)))
        alpha_b = 2.07*k_equality*sound_horizon*
     &            (1.+R_drag)**(-0.75)*alpha_b


        beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*
     &           sqrt((17.2*omhh)**2.+1e0)

        beta_node = 8.41*omhh**(0.435)

        return

        end


        function TFtransfer(k0)
c        subroutine TFtransfer_function(k,omhh,f_baryon,tf_full,
c     &                tf_baryon,tf_cdm)

c  Calculate transfer function from the fitting parameters stored in
c  GLOBALVARIABLES.
c
c  Input:
c        k -- wavenumber in Mpc^{-1} 
c        omhh -- The density of CDM and baryons, in units of critical dens,
c                multiplied by the square of the Hubble constant, in units
c                of 100 km/s/Mpc */
c        f_baryon -- The fraction of baryons to CDM */
c      
c  Output:
c        tf_full -- The full fitting formula, eq. (16), for the matter
c                   transfer function.
c        tf_baryon -- The baryonic piece of the full fitting formula, eq. 21.
c        tf_cdm -- The CDM piece of the full fitting formula, eq. 17.
c


        implicit none
        real*8 TFtransfer,k0
        real k,tf_full,tf_baryon,tf_cdm,q,ks
        real theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality,
     &       sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,
     &       f_baryon,beta_node
        common/GLOBALVARIABLES/theta_cmb,z_equality,k_equality,z_drag,
     &       R_drag,R_equality,sound_horizon,k_silk,alpha_c,beta_c,
     &       alpha_b,beta_b,beta_node
        real omhh,Tcmb,TF_pressureless,s_tilde
        common/GLOBALVAR2/omhh,f_baryon,Tcmb


        k=k0

cc  Reasonable k?
c
c        if (k.le.0) then
c           write(6,*) 'TFtransfer_function(): Illegal k'
c           pause
c        end if


c  Auxiliary Variables

            q = k/13.41/k_equality
            ks = k*sound_horizon

c  Main Variables

              tf_cdm = 1./(1.+(ks/5.4)**4.)
            tf_cdm = tf_cdm*TF_pressureless(q,1.,beta_c) +
     &           (1.-tf_cdm)*TF_pressureless(q,alpha_c,beta_c)


              s_tilde = sound_horizon/(1.+(beta_node/ks)**3.)**(1./3.)
              tf_baryon = TF_pressureless(q,1.,1.)/(1.+(ks/5.2)**2.)
              tf_baryon = tf_baryon + alpha_b/(1.+(beta_b/ks)**3)
     &                       *exp(-(k/k_silk)**(1.4))
              tf_baryon = tf_baryon *(sin(k*s_tilde)/(k*s_tilde))
            tf_full = f_baryon*tf_baryon + (1-f_baryon)*tf_cdm

            TFtransfer=tf_full   !!!

         return

         end

c       auxiliary function: Pressureless TF

        real function TF_pressureless(q,a,b)

          implicit none
          real q,a,b

          TF_pressureless = Log(exp(1.)+1.8*b*q)
          TF_pressureless = TF_pressureless/(TF_pressureless +
     &                      (14.2/a + 386/(1.+69.9*q**1.08))*q**2)

          return

        end    


c====End of Eisenstein & Hu 1998 subroutines=============================




c################ SUBROUTINTE GROUP2: on time ###########################

c========================================================================

      subroutine deltacrit(zredshift,omega,omega_lambda,delta_cz)
c---------------------------------------------------------------------
c A. Jenkins routine.
c Modifies the parameter delta_c. For Open universes we take a factor
c from Lacey \& Cole 1993, and for Omega+Lambda universes apply
c a fit to results from Eke et al 1996.
c---------------------------------------------------------------------

      implicit none
      real*8 omega,omega_lambda,zredshift,delta_cz,fdeltac
      real*8 omega_z,rfac,ling

        if (omega.eq.1.0.and.omega_lambda.eq.0)then
           fdeltac = 1.0
        elseif (omega+omega_lambda.eq.1.0) then
            omega_z = omega*(1.+zredshift)**3/
     1                (omega*(1.+zredshift)**3 + omega_lambda)            

            fdeltac = 1. - 0.0052*(1.-omega_z) - 0.009*(1.-omega_z)**3
     1                - 0.01*(1.-omega_z)**18
        elseif (omega.lt.1.0.and.omega_lambda.eq.0) then
            omega_z = omega*(1.+ zredshift)/(1.+omega*zredshift)
            fdeltac = rfac(omega_z)
        else
           stop 'deltacrit: unknown option'
        endif

        call growth_factor(omega,omega_lambda,zredshift,ling)
 
        delta_cz = 1.6864702*fdeltac/ling

      end


      real*8 function rfac(omega)
c---------------------------------------------------
c  Code to compute delta_c for Open universe.
c---------------------------------------------------

      implicit none
      real*8 omega,delta_c,pi,ceta,seta,eta
      parameter (delta_c=1.6864702,pi=3.1415926535)
      ceta = 2./omega-1
      seta = sqrt(ceta**2-1.)
      eta = dlog(ceta+seta)

      rfac = 1.5/delta_c*(1.+(2.*pi/(seta-eta))**(2./3.))*
     1  (3.*seta*(seta-eta)/(ceta-1.)**2-2.)

      return
      end 

c=============================================================================

      subroutine growth_factor(omega,lambda,z,ling)
C------------------------Routines to calculate linear growth factor ----------
c------------------------by V.R. Eke -----------------------------------------

      implicit none
      real*8 a,omega,lambda,z,lingro,ling
      a = 1./(1.+ z)
      ling = lingro(a,omega,lambda)
      end


      function lingro(a,omega,lambda)
c----------------------------------------------------------------------------
c     To calculate the linear growth factor D(a) for different cosmological 
c     models. Normalised such that D(1) = 1. (Doesn't include closed models 
c     or lambda models where omega+lambda isn't one.)
c----------------------------------------------------------------------------

      implicit none
      real*8 a,omega,lambda,func0,func1,func2,x,int,aofx,aofxn
      real*8 xn,w,dn,lingro,sum
      external func0,func1,func2

      w = omega**(-1.0) - 1.0
      sum = omega + lambda
      if (sum.gt.1 .or. omega.le.0 .or.sum.ne.1.and.lambda.gt.0) then
       if (abs(sum-1.0d0).gt.1.e-10) then
         write(*,*) 'Cannot cope with this cosmology!'
         stop
          endif
      endif
      if (omega .eq. 1) then
         lingro = a
      else if (lambda .gt. 0) then
         xn = (2.0*w)**(1.0/3)
         call simp(func2,0.0D0,xn,int)
         aofxn = ((xn**3.0+2.0)**0.5)*(int/xn**1.5)
         if (a .eq. 1.64) write(*,*) xn,aofxn
         x = a*xn
         call simp(func2,0.0D0,x,int)
         aofx = ((x**3+2)**0.5)*(int/x**1.5)
         lingro = aofx/aofxn
         if (a .eq. 1.64) write(*,*) x,aofx,lingro
      else
         dn = func1(w)
         x = w*a
         lingro = func1(x)/dn
      endif
 
      end

      function func0(x)
      real*8 x,func0
      func0 = 3/x + (3*((1+x)**0.5)/x**1.5)*dlog((1+x)**0.5-x**0.5)
      end

      function func1(x)
      real*8 x,func1
      func1 = 1 + 3/x + (3*((1+x)**0.5)/x**1.5)*dlog((1+x)**0.5-x**0.5)
      end

      function func2(x)
      real*8 x,func2
      func2 = (x/(x**3+2))**1.5
      end

      subroutine simp(func,a,b,s)
      implicit none
      real*8 func,a,b,eps,ost,os,s,st
      integer*4 jmax,j
      external func
      ost = -1.e-30
      os = -1.e-30
      eps = 1.e-5
      jmax = 20
      do j = 1,jmax
         call trapzd(func,a,b,st,j)
         s = (4.*st - ost)/3.
         if (abs(s-os) .lt. eps*abs(os)) goto 3
         os = s
         ost = st
      enddo
      stop 'Too many steps.'
3     end


      
      function omegaz(omega0,lambda0,z)
c---------------------------------------------------
c Caculates omega(z) - Donghai Zhao, 2001.
c---------------------------------------------------
      implicit none
      real*8 omega0,lambda0,z,ztp1,omegaz,ez2

      ztp1=z+1.
      omegaz=omega0*ztp1**3/ez2(omega0,lambda0,z)

      end

      function lambdaz(omega0,lambda0,z)
c---------------------------------------------------
c Caculates lambda(z) - Donghai Zhao, 2001.
c---------------------------------------------------
      implicit none
      real*8 omega0,lambda0,z,ztp1,lambdaz,ez2

      ztp1=z+1.
      lambdaz=lambda0/ez2(omega0,lambda0,z)

      end

      function ez2(omega0,lambda0,z)
c---------------------------------------------------
c Caculates (H(z)/H0)^2 - Donghai Zhao, 2001.
c---------------------------------------------------
      implicit none
      real*8 omega0,lambda0,z,ztp1,ez2

      ztp1=z+1.
      ez2=((omega0*ztp1+1.-lambda0
     &     -omega0)*ztp1**2+lambda0)

      end

      subroutine halodensv(omega0,lambda0,z,densh)
c---------------------------------------------------
c Caculates virial halo density - Donghai Zhao, 2001.
c---------------------------------------------------
      implicit none
      real*8 omega0,lambda0,z
      real*8 PI,omegat,xomega,omegaz,ztp1,con,densh
      parameter(PI=3.141592654)

      ztp1=z+1.
      omegat=omegaz(omega0,lambda0,z)
      xomega=omegat-1.

      if((omega0+lambda0).eq.1)then
         con=(18.*PI**2.+82.*xomega-39.*xomega**2.)/omegat
      elseif(lambda0.eq.0)then
         con=(18.*PI**2.+60.*xomega-32.*xomega**2.)/omegat
      else
         stop 'neither flat nor zero lambda'
      endif
      densh=2.7754e+11*omega0*ztp1**3.*con

      end


c This routine calculates the age of the universe at
c redshift z for a general FRW model.
c    omega0--the density parameter at x=1.;
c    lambda0-the vaccum density parameter at x=1.;
c    h0--------Hubble constant in km/s/Mpc;
c    z---------the redshift
c
c Output variables:
c    uniage1d----age in years;
c
c Other routines:
c    qromb.f-integration routine
c
      function uniage1d(omega0,lambda0,h0,z)
      implicit none
      real*8 lambda0,lambda01,omega0,omega01,h0,z
      real*8 uniage1d,ul,dl,ctfunc
      common /fluc/omega01,lambda01
      external ctfunc

      omega01=omega0
      lambda01=lambda0
      ul=0.
      dl=1./(1.+z)**1.5
      call qromb(ctfunc,ul,dl,uniage1d)
      uniage1d=uniage1d*9.78e11/h0*2./3.
      end

      function ctfunc(y)
      implicit none
      real*8 omega0,lambda0
      real*8 ctfunc,y
      common /fluc/omega0,lambda0

      ctfunc=1./sqrt(omega0+lambda0*y**2+
     &     (1.-omega0-lambda0)*y**(2./3.))
      end

c========integrator=======================================

      SUBROUTINE qromb(func,a,b,ss)
      implicit none
      INTEGER*4 JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS             ! zhao
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER*4 j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)  !zhao
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      stop 'too many steps in qromb'
      END

      SUBROUTINE trapzd(func,a,b,s,n)
      implicit none
      INTEGER*4 n
      REAL*8 a,b,s,func          ! zhao
      EXTERNAL func
      INTEGER*4 it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END

      SUBROUTINE polint(xa,ya,n,x,y,dy)
      implicit none
      INTEGER*4 n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)       ! zhao
      PARAMETER (NMAX=10)
      INTEGER*4 i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)stop 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

c====================================================





c################ SUBROUTINTE GROUP3: on NFW function ###################

c========================================================================
c by Donghai Zhao

       subroutine rvh(m,rho,r,v)
c m------mass in [M_sun/h]
c r------radius in [Mpc/h]
c rho----average density in [h^2 M_sun/Mpc^3]
c v------circular velocity in [km/s]
c G------gravity constant in [Mpc/M_sun (km/s)^2]

       implicit none
       real*8 PI,G,m,r,rho,v
       parameter (PI=3.141592654,G=4.302d-9)

       r=(3/PI/4*m/rho)**(1.0d0/3)
       v=(G*m/r)**0.5d0

       end

c========================================================================
c This subroutine calculates some characteristic inner properties for
c a NFW halo of given mass and concentration.
c by Donghai Zhao

       subroutine rvs(m,rho,c,ms,rs,rhos,vs)
c m------mass in [M_sun/h]
c c------NFW concentration
c r------radius in [Mpc/h]
c rs-----NFW characteristic inner radius in [Mpc/h]
c rho----average density in [h^2 M_sun/Mpc^3]
c rhos---density at rs in [h^2 M_sun/Mpc^3]
c v------circular velocity in [km/s]
c vs-----circular velocity at rs in [km/s]
c G------gravity constant in [Mpc/M_sun (km/s)^2]

       implicit none
       real*8 PI,G,m,r,rho,v,c,ms,rs,rhos,vs
       parameter (PI=3.141592654,G=4.302d-9)

       call rvh(m,rho,r,v)
       ms=m*(dlog(2.0d0)-0.5)/(dlog(1+c)-c/(1+c))
       rs=r/c
       vs=(G*ms/rs)**0.5d0
       rhos=rho*c**3/(dlog(1+c)-c/(1+c))/12

       end

c========================================================================
       subroutine densinr
c------------------------------------------------------------------------------
c Calculates average density as a function of r/rs, in [average density in rs].
c by Donghai Zhao in Jan., 2007
c------------------------------------------------------------------------------
       implicit none
       real*8 lgxmin,lgxmax,dlgx,lgx,x,const
       integer*4 nsample,nsamplet,i
       parameter (nsample=2000)
       real*8 lgdensinr
       common /densinr1/lgxmin,lgxmax,dlgx,lgdensinr(0:nsample),nsamplet

       nsamplet=nsample
       lgxmin=-2.0d0
       lgxmax=8.0d0
       dlgx=(lgxmax-lgxmin)/nsample
       lgxmax=lgxmin+dlgx*nsample
       const=1./(dlog(2.0d0)-0.5d0)

       do i=0,nsample
          lgx=lgxmin+dlgx*i
          x=10.0d0**lgx
          lgdensinr(i)=dlog10(const*(dlog(1.+x)-x/(1.+x))/x**3)
       enddo

       end


c========================================================================
c This subroutine convert concentration between different halo definitions.
c by Donghai Zhao in Jan., 2007
c
c Input:
c    dens1   --  halo definition density 1
c    dens2   --  halo definition density 2
c    c1      --  concentration corresponding to definition 1
c Output:
c    c2      --  concentration corresponding to definition 2
c    rmass21 --  mass of definiton 2 devided by mass of definition 1
c
       subroutine cconvert(dens1,dens2,c1,c2,rmass21)
       implicit none
       integer*4 nsample,nsample1,NMAX
       parameter(NMAX=2001)
       real*8 dens1,dens2,c1,c2,rmass21
       real*8 lgc1,lgc2,lgdinr1,lgdinr2,lgxmin,lgxmax,dlgx,yacc,dydx
       real*8 lgdensinr
       common /densinr1/ lgxmin,lgxmax,dlgx,lgdensinr(NMAX),nsample

       yacc=1.0d-5
       lgc1=dlog10(c1)
       nsample1=nsample+1
       if(nsample1.ne.NMAX)stop 'nsample1.ne.NMAX in cconvert'
       call ftable(lgc1,lgdinr1,dydx,lgxmin,lgxmax,nsample1,dlgx,
     &   lgdensinr)
       lgdinr2=lgdinr1+dlog10(dens2/dens1)
       lgc2=lgc1                               !guess
       call findx(lgc2,lgdinr2,yacc,lgxmin,lgxmax,nsample1,dlgx,
     &   lgdensinr)
       c2=10.0d0**lgc2

       rmass21=dens2*c2**3/(dens1*c1**3)

       end
       

c################ SUBROUTINTE GROUP4: on table ####################

      subroutine init_read(inf1,xmin,xmax,nsample1,dx,NMAX,yy)
c by Donghai Zhao
      implicit none
      integer*4 i,nsample,nsample1,NMAX
      real*8 xmin,xmax,dx,dum
      real*8 yy(NMAX)
      character inf1*(*)

      open(27,file=inf1,status="old")
      read(27,*)xmin,xmax,nsample,dx
      nsample1=nsample+1
      if(nsample1.gt.NMAX)stop 'nsample1.gt.NMAX in init_read'
      do i=1,nsample1
         read(27,*)dum,yy(i)
      enddo
      close(27)

      end


      subroutine findx(xguess,y0,yacc,xmin,xmax,nsample1,dx,yy)
c by Donghai Zhao
      implicit none
      real*8 xguess,diff,y0,yacc,xmin
      integer*4 nsample1
      real*8 xmax,dx,yy(nsample1),y,dydx

      diff=2*yacc
      do while(abs(diff).gt.yacc)
         call ftable(xguess,y,dydx,xmin,xmax,nsample1,dx,yy)
         diff=y0-y
         xguess=diff/dydx+xguess
      enddo
      end


      subroutine ftable(x,y,dydx,xmin,xmax,nsample1,dx,yy)
c by Donghai Zhao
      implicit none
      real*8 x,xmin,xmax,dx,y,dydx,dydxt
      integer*4 nsample1,ibin1,ibin2
      real*8 yy(nsample1)
      real*8 x1

      ibin1=(x-xmin)/dx+1
      ibin2=ibin1+1
      if(ibin1.ge.1.and.ibin2.le.nsample1)then
         x1=xmin+(ibin1-1)*dx
         dydx=(yy(ibin2)-yy(ibin1))/dx
         y=yy(ibin1)+dydx*(x-x1)
         dydxt=(yy(nsample1)-yy(1))/dx/(nsample1-1)
      else
         write(*,*)ibin1,x,xmin,xmax,nsample1
         stop 'x out of range in ftable'
      endif

      end




