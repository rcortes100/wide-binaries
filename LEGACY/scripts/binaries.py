#!/usr/bin/env python3
#
# Modulo para estrellas binarias
# Ricardo Adán Cortés Martín
# 8 de julio de 2019
#
# NOTAS:
# 2019.09.05 - La función RMS_histogram_1D tiene mal implementado el montecarlo.
#
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inp
import scipy.stats as stats

#=======================================
# Fonts
#=======================================
plt.rcParams['font.family'] = 'serif'

#=======================================
# Variables
#=======================================

kappa = 1000*(np.pi/(180*3600*1000))*(3.086E13)/(365*24*3600.)

#=======================================
# Geometric functions
#=======================================

def del_ra(x,y):
    a = np.abs(x-y)
    b = np.abs(a-360.0)
    return np.minimum(a,b)

def pm2kms(pm,dist):
    """Proper motions in mas/yr to km/s, using distance in pc"""
    return (np.pi/(180*3600*1000))*(1./(365*24*3600))*(3.086E13)*pm*dist
def pm2kms_error(pm,pm_err,dist,dist_err):
    """Proper motions error in mas/yr to km/s, using distance in pc"""
    return pm2kms(pm,dist)*(pm_err/pm + dist_err/dist)

def vel_tan_diff(pmra_1,pmdec_1,dist_1,pmra_2,pmdec_2,dist_2):
    """Tangential velocity difference using pm in ra and dec directions, and distances in pc,
    in km/s
    Only valid for small angles.
    """
    return np.sqrt((pm2kms(pmra_2,dist_2)-pm2kms(pmra_1,dist_1))**2 +
                   (pm2kms(pmdec_2,dist_2)-pm2kms(pmdec_1,dist_1))**2)

def vel_tan_diff_err(pmra_1,pmra_1_err,pmdec_1,pmdec_1_err,dist_1,dist_1_err,
                     pmra_2,pmra_2_err,pmdec_2,pmdec_2_err,dist_2,dist_2_err):
    """Tangential velocity difference error using pm in ra and dec directions, and distances in pc,
    in km/s
    Only valid for small angles
    """
    diffvra  = np.abs((pm2kms(pmra_2,dist_2)-pm2kms(pmra_1,dist_1))*
                      (pm2kms_error(pmra_1,pmra_1_err,dist_1,dist_1_err)+pm2kms_error(pmra_2,pmra_2_err,dist_2,dist_2_err)))
    diffvdec = np.abs((pm2kms(pmdec_2,dist_2)-pm2kms(pmdec_1,dist_1))*
                      (pm2kms_error(pmdec_1,pmdec_1_err,dist_1,dist_1_err)+pm2kms_error(pmdec_2,pmdec_2_err,dist_2,dist_2_err)))
    return (diffvra+diffvdec)/vel_tan_diff(pmra_1,pmdec_1,dist_1,pmra_2,pmdec_2,dist_2)

def cos_distang(ra1,dec1,ra2,dec2):
    dra  = np.deg2rad(del_ra(ra1,ra2))
    dec1 = np.deg2rad(dec1)
    dec2 = np.deg2rad(dec2)
    return np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(dra)

def dcd(ra1,dec1,ra2,dec2,dra1,ddec1,dra2,ddec2):
    dra    = np.abs(np.deg2rad(del_ra(ra1,ra2)))
    cd1sd2 = np.abs(np.cos(np.deg2rad(dec1))*np.sin(np.deg2rad(dec2)))
    sd1cd2 = np.abs(np.sin(np.deg2rad(dec1))*np.cos(np.deg2rad(dec2)))
    dra1   = np.deg2rad(dra1)
    dra2   = np.deg2rad(dra2)
    ddec1  = np.deg2rad(ddec1)
    ddec2  = np.deg2rad(ddec2)
    
    return cd1sd2*ddec1 + sd1cd2*ddec2 + np.cos(dra)*(
    sd1cd2*ddec1 + cd1sd2*ddec2) + np.cos(np.deg2rad(dec1))*np.cos(np.deg2rad(dec2))*np.sin(dra)*(dra1+dra2)

def R3D(r1,r2,cost):
    """
    Separación tridimensional exacta entre dos puntos. Toma las distancias radiales de
    cada punto, y el coseno del ángulo de separación entre los puntos.
    """
    return np.sqrt((r2-r1)**2 + 2.*r1*r2*(1.-cost))

def dR3D(r1,r2,cost,dr1,dr2,dcost):
    """
    Error en la separación tridimensional
    """
    return (np.abs(r2-r1)*(dr1+dr2) + r2*r1*dcost + (1.-cost)*(r1*dr2 + r2*dr1))/R3D(r1,r2,cost)

############################################################################
# GAIA DR2 MODEL CORRECTIONS
############################################################################
def D_model_correction(ra,dec,G_band):
    """
    See Lindegren 2018
    """
    omegax = -0.059 #mas yr^-1
    omegay = -0.112 #mas yr^-1
    omegaz =  0.010 #mas yr^-1

    F = np.floor(-0.03333*(G_band-41)) \
       -np.floor(-0.033333*(G_band-11))*\
        np.floor(-0.0333333*(G_band-43))*0.5*(13-G_band)
#    if (G_band > 13):
#        F=0
#    elif (G_band <= 11):
#        F=1
#    else:
#        F = 0.5*(13 - G_band)
        
    return (F*omegax*np.sin(np.deg2rad(dec))*np.cos(np.deg2rad(ra))+
            F*omegay*np.sin(np.deg2rad(dec))*np.sin(np.deg2rad(ra))-
            F*omegaz*np.cos(np.deg2rad(dec)),
           -F*omegax*np.sin(np.deg2rad(ra)) + F*omegay*np.cos(np.deg2rad(ra)))

def D_model_correction_err(ra,dec,G_band):
    """
    See Lindegren 2018
    """
    domegax = 0.025 #mas yr^-1
    domegay = 0.029 #mas yr^-1
    domegaz = 0.030 #mas yr^-1
    
    F = np.floor(-0.03333*(G_band-41)) \
       -np.floor(-0.033333*(G_band-11))*\
        np.floor(-0.0333333*(G_band-43))*0.5*(13-G_band)
    
    return (F*domegax*np.sin(np.deg2rad(dec))*np.cos(np.deg2rad(ra))+
            F*domegay*np.sin(np.deg2rad(dec))*np.sin(np.deg2rad(ra))+
            F*domegaz*np.cos(np.deg2rad(dec)),
           +F*domegax*np.sin(np.deg2rad(ra)) + F*domegay*np.cos(np.deg2rad(ra)))



############################################################################
###SPHERICAL CORRECTION
############################################################################
def spherical_corr(ra_sec,dec_sec,dist_sec,pmra_sec,pmdec_sec,vrad_sec,
                   ra_pri,dec_pri,U=False,V=False,W=False,useUVW=False):
    """
Converts the proper motions (ra and dec in mas/yr) and radial velocity (in km/s) from the frame of reference
of the secondary to the frame of reference of the primary. pmra is defined as pmra*cos(dec) as in GAIA.
    """

    #Fill na with zeroes
    vrad_sec = vrad_sec.fillna(0)
    
    C = kappa*dist_sec*1E-3# Kappa/parallax
    
    if not useUVW:
        U = -C*pmra_sec*np.sin(np.deg2rad(ra_sec)) \
            -C*pmdec_sec*np.cos(np.deg2rad(ra_sec))*np.sin(np.deg2rad(dec_sec)) \
            +vrad_sec*np.cos(np.deg2rad(ra_sec))*np.cos(np.deg2rad(dec_sec))
        V =  C*pmra_sec*np.cos(np.deg2rad(ra_sec)) \
            -C*pmdec_sec*np.sin(np.deg2rad(ra_sec))*np.sin(np.deg2rad(dec_sec)) \
            +vrad_sec*np.sin(np.deg2rad(ra_sec))*np.cos(np.deg2rad(dec_sec))
        W = C*pmdec_sec*np.cos(np.deg2rad(dec_sec)) \
           +vrad_sec*np.sin(np.deg2rad(dec_sec))
    
    sinra  = np.sin(np.deg2rad(ra_pri))
    cosra  = np.cos(np.deg2rad(ra_pri))
    sindec = np.sin(np.deg2rad(dec_pri))
    cosdec = np.cos(np.deg2rad(dec_pri))

    #Return the corrected pmra, pmdec and vrad, as three different objects    
    return (-U*sinra + V*cosra)/C,\
           (-U*cosra*sindec - V*sinra*sindec + W*cosdec)/C,\
             U*cosra*cosdec + V*sinra*cosdec + W*sindec

def test_spherical_corr(ra_sec,dec_sec,dist_sec,pmra_sec,pmdec_sec,vrad_sec,
                        ra_pri,dec_pri,
                        U,V,W):
    """
Converts the proper motions (ra and dec in mas/yr) and radial velocity (in km/s) from the frame of reference
of the secondary to the frame of reference of the primary. pmra is defined as pmra*cos(dec) as in GAIA.
    """

    #Fill na with zeroes
    vrad_sec = vrad_sec.fillna(0)
    
    C = kappa*dist_sec*1E-3# Kappa/parallax
    
    sinra  = np.sin(np.deg2rad(ra_pri))
    cosra  = np.cos(np.deg2rad(ra_pri))
    sindec = np.sin(np.deg2rad(dec_pri))
    cosdec = np.cos(np.deg2rad(dec_pri))

    #Return the corrected pmra, pmdec and vrad, as three different objects    
    return (-U*sinra + V*cosra)/C,\
           (-U*cosra*sindec - V*sinra*sindec + W*cosdec)/C,\
             U*cosra*cosdec + V*sinra*cosdec + W*sindec

#========================================
# Selection functions
#========================================

def seleccion(good,sigma=3.0,par=0.16,rvel=3.):
    x = good.copy()
    return x[(np.abs(1.-x.pmRA_y/x.pmra_y) < sigma*np.abs(
                 (x.pmRA_y/x.pmra_y)*(x.e_pmRA_y/x.pmRA_y+x.pmra_error_y/x.pmra_y))) &
             (np.abs(1.-x.pmDE_y/x.pmdec_y) < sigma*np.abs(
                 (x.pmDE_y/x.pmdec_y)*(x.e_pmDE_y/x.pmDE_y+x.pmdec_error_y/x.pmdec_y))) &
#          (np.abs(1.-x.Plx_y/x.parallax_y) < sig*np.abs(
#              (x.Plx_y/x.parallax_y)*(x.e_Plx_y/x.Plx_y+x.parallax_error_y/x.parallax_y))) &
             (np.abs(1.-x.pmRA_x/x.pmra_x) < sigma*np.abs(
                 (x.pmRA_x/x.pmra_x)*(x.e_pmRA_x/x.pmRA_x+x.pmra_error_x/x.pmra_x))) &
             (np.abs(1.-x.pmDE_x/x.pmdec_x) < sigma*np.abs(
                 (x.pmDE_x/x.pmdec_x)*(x.e_pmDE_x/x.pmDE_x+x.pmdec_error_x/x.pmdec_x))) &
#          (np.abs(1.-x.Plx_x/x.parallax_x) < sig*np.abs(
#              (x.Plx_x/x.parallax_x)*(x.e_Plx_x/x.Plx_x+x.parallax_error_x/x.parallax_x))

             (x.parallax_y > (1.-par)*x.parallax_x) &
             (x.parallax_y < (1.+par)*x.parallax_x) &

             ((np.abs(x.radial_velocity_x-x.radial_velocity_y) < rvel) | x.radial_velocity_x.isna() | x.radial_velocity_y.isna()) &
#      ~x.Ind.isin([737,626,352,492,426,839,691,309,549,89,65,682,546,633,677,201,749,815,678]) &
             True]

#####################################################
#
# Aquí definimos las matrices para calcular el efecto de proyección
#
def matA(ra,dec):
    """
    Es la matriz correspondiente a la transformación del sistema de ecuaciones 5,6,7 de la p.16
    del libro Stellar Kinematics de Smart
    """
    return np.array([[-np.sin(ra)            , np.cos(ra)             , 0.],
                     [-np.cos(ra)*np.sin(dec), -np.sin(ra)*np.sin(dec), np.cos(dec)],
                     [ np.cos(ra)*np.cos(dec),  np.sin(ra)*np.cos(dec), np.sin(dec)]])

def matB(ra,dec):
    """
    Es la matriz correspondiente a la transformación del sistema de ecuaciones 8,9,10 de la p.16
    del libro Stellar Kinematics de Smart
    """        
    return np.array([[-np.sin(ra)*np.cos(dec), -np.cos(ra)*np.sin(dec), np.cos(ra)*np.cos(dec)],
                     [ np.cos(ra)*np.cos(dec), -np.sin(ra)*np.sin(dec), np.sin(ra)*np.cos(dec)],
                     [                     0.,             np.cos(dec), np.sin(dec)]])

#PREPROCESING
def pickleizar(input_file,output_file='main_list.pkl',skprws=47):
    pd.read_csv(input_file,
                skiprows=skprws,false_values=['*','          *'],
                na_values=['*','          *']).to_pickle(output_file)
    return 0

def rot_z(vector,angulo):
    """
    Rota el vector (x,y,z) alrededor del Eje Z, dado un ángulo en grados.
    """
    angulo = np.radians(angulo)
    return np.dot(vector,
                      [[np.cos(angulo),-np.sin(angulo),0],
                       [np.sin(angulo),np.cos(angulo),0],
                       [0,0,1]])

def rot_y(vector,angulo):
    """
    Rota el vector (x,y,z) alrededor del Eje Y, dado un ángulo en grados.
    """
    angulo = np.radians(angulo)
    return np.dot(vector,
                      [[np.cos(angulo),0,-np.sin(angulo)],
                       [0,1,0],
                       [np.sin(angulo),0,np.cos(angulo)]])
    

picklefile = 'jimenez-esteban2019.pkl'

#====================================================
# Tremaine data
#====================================================

rjacobi = 1.7 #pc (M1+M2/2Msol)^1/3
Omega_g_rj = 0.05 #km/s (M1+M2/2Msol)^1/3
#Coordenadas en papel milimétrico
tremaine=np.array([
    [ 10,35,   39,   53, 56,   67,77,  81,  85,  90,97,  100,  105,  110,  115,118,127,135,142,150,160,172,180,190,198,206,213,220],
    [160,135,131.2,117.9,115,104.5,95,91.2,87.4,82.6,76, 73.1, 68.4, 63.6, 58.8, 56,52,57,68,69,70,72,74,77,80,83,85,85]])

#Conversión de la escala del papel a valores
DIST_trem = 7./203*tremaine[0]-4.5  #7 dex(rp/rj)/203mm
RMS_trem = 3./167*tremaine[1]-1. # 3 dex(rms V en rJ*omegaG)/167mm

f1 = inp.interp1d(DIST_trem,RMS_trem)
f2 = inp.interp1d(DIST_trem,RMS_trem,kind='cubic')
xspan = np.linspace(DIST_trem.min(),DIST_trem.max(),100)
def trem_line(r,a=Omega_g_rj,b=rjacobi):
    return a*10**(f2(np.log10(r/b)))
#======================================================
# Histograms
#======================================================

def N2_model(x):
    return 2*0.5**(np.log10(x))
def N4_model(x):
    return 20*3**(-np.log10(x)/3)

def make_histogram(table,rmin=False,rmax=False,rv_max=10,bins=10,poiss=False,
                   proy_sep = 'proy_sep_GDR2'):

#    x = table[((np.abs(table.rv2 - table.rv1 + table.rv1_diff) < rv_max) & ~table.rv1.isna() & ~table.rv2.isna()) &
#              (table.N2 <= N2_max) &
#              (table.N4 <= N4_max)].copy()
#              (table.N2 <= N2_model(table.proy_sep)) &
#              (table.N4 <= N4_model(table.proy_sep))].copy()
    #pc*mas/yr to km/s
    x = table.copy()
    k = 1000*(np.pi/(180*3600*1000))*(3.086E13)/(365*24*3600.)
    
    
    if rmin==False:
        rmin = np.log10(x[proy_sep].min())
    if rmax==False:
        rmax = np.log10(x[proy_sep].max())

    step   = (rmax-rmin)/bins

    rms_ra=[]
    rms_ra_error=[]
    rms_dec=[]
    rms_dec_error=[]
    rms_tan=[]
    rms_tan_error=[]
    rms_rad=[]
    rms_rad_error=[]
    sep=[]
    samp=[]

    for i in range(bins):
        _temp = x[(x[proy_sep] > 10**(rmin+i*step)) &
                  (x[proy_sep] < 10**(rmin+(i+1)*step))]
#        print('Intervalo %f %f\n'%(10**(rmin+i*step),10**(rmin+(i+1)*step)))
    #CORRECION ESFERICA SOBRE LA PRIMARIA    

        pmrax = _temp.pmra_x - _temp.pmra_corr
        pmdex = _temp.pmdec_x - _temp.pmde_corr
        pmray = _temp.pmra_y
        pmdey = _temp.pmdec_y
        distx = 1000./_temp.parallax_x
        disty = 1000./_temp.parallax_y
        rvelx = _temp.radial_velocity_x - _temp.rvel_corr
        rvely = _temp.radial_velocity_y

#       pmrax = (_temp.v1_ra - _temp.v1_ra_diff)*1000/(k*_temp.d1)
#       pmdex = (_temp.v1_dec - _temp.v1_dec_diff)*1000/(k*_temp.d1)
#       pmray = _temp.v2_ra*1000/(k*_temp.d2)
#       pmdey = _temp.v2_dec*1000/(k*_temp.d2)
#       distx = _temp.d1
#       disty = _temp.d2

#       rvelx = _temp.rv1 - _temp.rv1_diff
#       rvely = _temp.rv2
        pmrax_err = _temp.pmra_error_x
        pmdex_err = _temp.pmdec_error_x
        pmray_err = _temp.pmra_error_y
        pmdey_err = _temp.pmdec_error_y
#        distx_err = distx*_temp.parallax_error_x/_temp.parallax_x
        distx_err = distx/_temp.parallax_over_error_x
#        disty_err = disty*_temp.parallax_error_y/_temp.parallax_y
        disty_err = disty/_temp.parallax_over_error_y
        rvelx_err = _temp.radial_velocity_error_x
        rvely_err = _temp.radial_velocity_error_y
    
#       pmrax_err = pmrax*(_temp.v1_ra_err/_temp.v1_ra   + _temp.d1_err/_temp.d1)
#       pmdex_err = pmdex*(_temp.v1_dec_err/_temp.v1_dec + _temp.d1_err/_temp.d1)
#       pmray_err = pmray*(_temp.v2_ra_err/_temp.v2_ra   + _temp.d2_err/_temp.d2)
#       pmdey_err = pmdey*(_temp.v2_dec_err/_temp.v2_dec + _temp.d2_err/_temp.d2)
#       distx_err = _temp.d1_err
#       disty_err = _temp.d2_err
#       rvelx_err = _temp.rv1_err
#       rvely_err = _temp.rv2_err
        
        rmsa = np.sqrt((vel_tan_diff(pmrax,0.,distx,
                                     pmray,0.,disty)**2).mean())
        rmsa_error = (np.abs(
            vel_tan_diff(pmrax,0.,distx,pmray,0.,disty)*
            vel_tan_diff_err(pmrax,pmrax_err,1E-20,1E-30,distx,distx_err,
                             pmray,pmray_err,1E-20,1E-30,disty,disty_err))).mean()/rmsa
#    rmsa_error = c**2*(np.abs((disty*pmray-distx*pmrax)*(
#                                 disty_err*pmray+
#                                 disty*pmray_err+
#                                 distx_err*pmrax+
#                                 distx*pmrax_err)).mean())/rmsa

        rmsd = np.sqrt((vel_tan_diff(0.,pmdex,distx,
                                     0.,pmdey,disty)**2).mean())
        rmsd_error = (np.abs(
            vel_tan_diff(0.,pmdex,distx,0.,pmdey,disty)*
            vel_tan_diff_err(1E-20,1E-30,pmdex,pmdex_err,distx,distx_err,
                             1E-20,1E-30,pmdey,pmdey_err,disty,disty_err))).mean()/rmsd

        rmst = np.sqrt(rmsa**2+rmsd**2)
        
        rmst_error = (np.abs(
            vel_tan_diff(pmrax,pmdex,distx,pmray,pmdey,disty)*
            vel_tan_diff_err(pmrax,pmrax_err,pmdex,pmdex_err,distx,distx_err,
                             pmray,pmray_err,pmdey,pmdey_err,disty,disty_err))).mean()/rmst
    
        rmsr = np.sqrt(((rvely-rvelx)**2).mean())
        rmsr_error = np.abs((rvely-rvelx)*(rvelx_err+rvely_err)).mean()/rmsr
    
        # Agregar error poissoniano:
        
        if poiss:
            rmsa_error = np.sqrt(_temp.shape[0])*rmsa_error
            rmsd_error = np.sqrt(_temp.shape[0])*rmsd_error
    
        rms_ra.append(rmsa)
        rms_ra_error.append(rmsa_error)
        rms_dec.append(rmsd)
        rms_dec_error.append(rmsd_error)
        rms_tan.append(rmst)
        rms_tan_error.append(rmst_error)
        rms_rad.append(rmsr)
        rms_rad_error.append(rmsr_error)
        sep.append(rmin+(i+0.5)*step)
        samp.append(_temp.shape[0])
    
#        print(rmsa, rmsa_error)
#        print(rmsd, rmsd_error)
#    print(rmst, rmst_error)
#        print(step)
#        print('x')
#    print(sep)
    return sep,samp,step,rms_ra,rms_ra_error,rms_dec,rms_dec_error,x

def RMS_histogram_1D(dataframe,rmin=False,rmax=False,bins=10,poiss=False,mc=True, N=999,
                     seplabel='proy_sep',vellabel='vx1',velerrlabel='vx1_err'):

    """
Makes an RMS histogram for the proyected velocities of binaries in function of a projected separation
    
Parameters
----------
  dataframe : Dataframe like with velocities, separations and velocity errors.
       rmin : Base 10 logarithm of the minimum separation. If not stated, automatically uses the minimum
              separation available in dataframe.
       rmax : Base 10 logarithm of the maximum separation. If not stated, automatically uses the maximum
              separation available in dataframe.
       bins : Number of bins of the histogram.
      poiss : Adds Poisson error, when the sample is small. Default is False.
   seplabel : Label of the projected separations in the dataframe.
   vellabel : Label of the velocity differences in the dataframe.
velerrlabel : Label of the velocity differences errors in the dataframe.
     
    
    Returns:
        rmin = log10 of minimun separation,
        rmax = log10 of maximum separation,
        step,sep,samp_x,samp_y,samp_z, rms_x,rms_x_error,rms_y,rms_y_error,rms_z,rms_z_error, l
    """
    
    l = dataframe
    k = 1000*(np.pi/(180*3600*1000))*(3.086E13)/(365*24*3600.)
    
    if rmin==False:
        rmin = np.log10(l[seplabel].min())
    if rmax==False:
        rmax = np.log10(l[seplabel].max())

    step   = (rmax-rmin)/bins

    rms       = np.array([])
    rms_error = np.array([])
    sep       = np.array([])
    samp      = np.array([])

    #Calculation of each bin
    for i in range(bins):
        _temp = l[(l[seplabel] > 10**(rmin+i*step)) &
                  (l[seplabel] < 10**(rmin+(i+1)*step))]

        #Calculate the RMS of each bin
        
        if mc:
            #Hacer N experimentos
            rms_mc = np.array([])
            
            for j in range(N):
                rms_mc = np.append(rms_mc,np.sqrt(
                    (np.random.normal(_temp[vellabel],_temp[velerrlabel])**2).mean()))
                
            rms1D     = rms_mc.mean()
            rms1D_err = rms_mc.std()
            
        else:
            rms1D     = np.sqrt( (_temp[vellabel]**2).mean() )
            rms1D_err = np.abs(_temp[vellabel]*_temp[velerrlabel]).mean()/(rms1D*np.sqrt(len(_temp)))

#            rms1D_err = np.sqrt(((np.abs(_temp[vellabel]) + _temp[velerrlabel])**2).mean()) -\
#                        np.sqrt(((np.abs(_temp[vellabel]) - _temp[velerrlabel])**2).mean())

        # Agregar error poissoniano:        
        if poiss:
            rms1D_err = np.sqrt(_temp.shape[0])*rms1D_err

        #Agregar puntos a los histogramas
        rms       = np.append(rms,rms1D)
        rms_error = np.append(rms_error,rms1D_err)
        
        sep       = np.append(sep,rmin+(i+0.5)*step)
#        sep.append(rmin+(i+0.5)*step)
        samp      = np.append(samp,_temp.shape[0])
#        samp.append(_temp.shape[0])
    
    return rmin,rmax,step,sep,samp, rms,rms_error

def make_histogram_XYZ(table,rmin=False,rmax=False,rv_max=10,bins=10,poiss=False,montecarlo=False,
                   proy_sep = 'proy_sep_GDR2'):

    """
        Table format:    Positions in (pc):     x1,y1,z1 for pri
                                                x2,y2,z2 for sec
                         Velocities in (km/s):  vx1,vy1,vz1 for pri
                                                vx2,vy2,vz2 for sec
                         And errors with suffix _err
    Returns:
        rmin = log10 of minimun separation,
        rmax = log10 of maximum separation,
        step,sep,samp_x,samp_y,samp_z, rms_x,rms_x_error,rms_y,rms_y_error,rms_z,rms_z_error, l
    """
    
    l = table.copy()
    k = 1000*(np.pi/(180*3600*1000))*(3.086E13)/(365*24*3600.)
    
    l['xy_sep'] = np.sqrt((l.x2-l.x1)**2+(l.y2-l.y1)**2)
    l['yz_sep'] = np.sqrt((l.y2-l.y1)**2+(l.z2-l.z1)**2)
    l['zx_sep'] = np.sqrt((l.z2-l.z1)**2+(l.x2-l.x1)**2)
    
    l['dvx'] = l.vx2-l.vx1
    l['dvy'] = l.vy2-l.vy1
    l['dvz'] = l.vz2-l.vz1
    
    l['dvx_err'] = l.vx2_err+l.vx1_err
    l['dvy_err'] = l.vy2_err+l.vy1_err
    l['dvz_err'] = l.vz2_err+l.vz1_err

    if rmin==False:
        rmin = np.log10(min(l.xy_sep.min(),
                            l.yz_sep.min(),
                            l.zx_sep.min()))
    if rmax==False:
        rmax = np.log10(max(l.xy_sep.max(),
                            l.yz_sep.max(),
                            l.zx_sep.max()))

    step   = (rmax-rmin)/bins

    rms_x       = np.array([])
    rms_x_error = np.array([])
    rms_y       = np.array([])
    rms_y_error = np.array([])
    rms_z       = np.array([])
    rms_z_error = np.array([])
    sep         = np.array([])
    samp_x      = np.array([])
    samp_y      = np.array([])
    samp_z      = np.array([])

    _dump,_dump,_dump,sep,samp_x,rms_x,rms_x_error = RMS_histogram_1D(l,rmin,rmax,bins,poiss,mc=montecarlo,
                                                                   seplabel='yz_sep',vellabel='dvx',velerrlabel='dvx_err')
    _dump,_dump,_dump,sep,samp_y,rms_y,rms_y_error = RMS_histogram_1D(l,rmin,rmax,bins,poiss,mc=montecarlo,
                                                                   seplabel='zx_sep',vellabel='dvy',velerrlabel='dvy_err')
    _dump,_dump,_dump,sep,samp_z,rms_z,rms_z_error = RMS_histogram_1D(l,rmin,rmax,bins,poiss,mc=montecarlo,
                                                                   seplabel='xy_sep',vellabel='dvz',velerrlabel='dvz_err')    
    
    return rmin,rmax,step,sep, samp_x,samp_y,samp_z, rms_x,rms_x_error,rms_y,rms_y_error,rms_z,rms_z_error, l

def make_histogram_spherical(table,rmin=False,rmax=False,rv_max=10,bins=10,poiss=False,montecarlo=True,N=999,
                             proy_sep = 'proy_sep_GDR2',
                             ra1 = 'ra2', ra2 = 'ra2',
                             dec1 = 'dec1', dec2 = 'dec2',
                             dist1 = 'dist1', dist2 = 'dist2',
                             dist1_err = 'dist1_err', dist2_err = 'dist2_err',
                             pmra1 = 'pmra1', pmra2 = 'pmra2',
                             pmra1_err = 'pmra1_err', pmra2_err = 'pmra2_err',
                             pmdec1 = 'pmdec1', pmdec2 = 'pmdec2',
                             pmdec1_err = 'pmdec1_err', pmdec2_err = 'pmdec2_err',
                             rv1 = 'rv1', rv2 = 'rv2',
                             
                             pmra2_corr = 'pmra2_corr',   pmra2_corr_err = 'pmra2_corr_err',
                             pmdec2_corr = 'pmdec2_corr', pmdec2_corr_err = 'pmdec2_corr_err',
                             rv2_corr = 'rv2_corr',
    
                             make_correction=False):
    """
        Table format:    Positions in (pc):     x1,y1,z1 for pri
                                                x2,y2,z2 for sec
                         Velocities in (km/s):  vx1,vy1,vz1 for pri
                                                vx2,vy2,vz2 for sec
                         And errors with suffix _err
    Returns:
        rmin = log10 of minimum separation,
        rmax = log10 of maximum separation,
        step,sep,samp_x,samp_y,samp_z, rms_x,rms_x_error,rms_y,rms_y_error,rms_z,rms_z_error, l
    """    
# Verificar que table contenga ra1,dec1,dist1,pmra1,pmdec1,rv1
#                              ra2,dec2,dist2,pmra2,pmdec2,rv2

    l = table.copy()

    if make_correction:
        l[pmra2_corr],l[pmdec2_corr],l[rv2_corr] = \
        spherical_corr(l[ra2],l[dec2],l[dist2], \
                       l[pmra2],l[pmdec2],l[rv2],
                       l[ra1],l[dec1])

    if proy_sep not in l.columns:
        l[proy_sep] = R3D(0.5*(l[dist1]+l[dist2]),0.5*(l[dist1]+l[dist2]),
                          cos_distang(l[ra1],l[dec1],l[ra2],l[dec2]))
    
    if rmin==False:
        rmin = np.log10(l[proy_sep].min())
    if rmax==False:
        rmax = np.log10(l[proy_sep].max())

    step   = (rmax-rmin)/bins

#    rms_ra        = np.array([])
#    rms_ra_error  = np.array([])
#    rms_dec       = np.array([])
#    rms_dec_error = np.array([])
#    sep           = np.array([])
#    samp_ra       = np.array([])
#    samp_dec      = np.array([])
    
    l['dv_ra']      = 1E-3*kappa*np.abs(l[pmra2_corr]*l[dist2]  - l[pmra1]*l[dist1])
    l['dv_dec']     = 1E-3*kappa*np.abs(l[pmdec2_corr]*l[dist2] - l[pmdec1]*l[dist1])
    l['dv_ra_err']  = 1E-3*kappa*np.abs(np.abs(l[pmra2_corr]*l[dist2_err]) +
                                        np.abs(l[pmra2_err]*l[dist2]) +
                                        np.abs(l[pmra1]*l[dist1_err]) +
                                        np.abs(l[pmra1_err]*l[dist1]))
    l['dv_dec_err'] = 1E-3*kappa*np.abs(np.abs(l[pmdec2_corr]*l[dist2_err]) +
                                        np.abs(l[pmdec2_err]*l[dist2]) +
                                        np.abs(l[pmdec1]*l[dist1_err]) +
                                        np.abs(l[pmdec1_err]*l[dist1]))
    if False:
        l['dv_ra']      = vel_tan_diff(l[pmra1],0,l[dist1],
                                       l[pmra2_corr],0,l[dist2])
        l['dv_dec']     = vel_tan_diff(0,l[pmdec1],l[dist1],
                                       0,l[pmdec2_corr],l[dist2])
        l['dv_ra_err']  = vel_tan_diff_err(l[pmra1],l[pmra1_err],1E-20,1E-30,l[dist1],l[dist1_err],
                                           l[pmra2_corr],l[pmra2_err],1E-20,1E-30,l[dist2],l[dist2_err])
        l['dv_dec_err'] = vel_tan_diff_err(1E-20,1E-30,l[pmdec1],l[pmdec1_err],l[dist1],l[dist1_err],
                                            1E-20,1E-30,l[pmdec2_corr],l[pmdec2_err],l[dist2],l[dist2_err])
        
    _dump,_dump,_dump,sep,samp_ra,rms_ra,rms_ra_error = \
    RMS_histogram_1D(l,rmin,rmax,bins,poiss,montecarlo,N,
                     seplabel=proy_sep, vellabel='dv_ra',velerrlabel='dv_ra_err')
    _dump,_dump,_dump,sep,samp_dec,rms_dec,rms_dec_error = \
    RMS_histogram_1D(l,rmin,rmax,bins,poiss,montecarlo,N,
                     seplabel=proy_sep,vellabel='dv_dec',velerrlabel='dv_dec_err')
    
    return rmin,rmax,step,sep, samp_ra,samp_dec, rms_ra,rms_ra_error,rms_dec,rms_dec_error, l

#Histogram maker
def RMS_histogram(sep,samp,step,rms_ra,rms_ra_error,rms_dec,rms_dec_error,
                  save='',figname='',fnt=10,
                  xlabl='r(AU)',ylabl='$<\Delta v^2>^{1/2}$(km s$^{-1}$)',
                  _xmin=10**(-3.5),
                  _xmax=10**(0.5),
                  _ymin=0.01,
                  _ymax=10.,
                  xsolid= 206000*rjacobi*10**(xspan),
                  ysolid=Omega_g_rj*10**(f2(xspan)),
                  showbins=False,
                  MONDline=7000,
                  showlegend=False):
    """
sep:           límite izquierdo del bin
samp:          miembros del bin
step:          separación entre bins
rms_ra:        RMS en ascención recta
rms_ra_error:  error del RMS en ascención recta
rms_dec:       RMS en declinacón
rms_dec_error: error del RMS en declinación
    """

#####    fig = plt.figure(figsize=(6,5),dpi=150)
    fig, (a0,a1) = plt.subplots(2, 1,
                                gridspec_kw={'height_ratios': [6, 1]},
                                dpi=150, figsize=(6,6))
    a0.loglog()
#    plt.semilogx()
    
    if figname:
        a0.set_title(figname)
    
    a0.set_xlabel(xlabl,fontsize=fnt)
    a0.set_ylabel(ylabl,fontsize=fnt)
    
    a0.set_xlim(_xmin,_xmax)
    a0.set_ylim(_ymin,_ymax)
    
    a0.plot(xsolid,ysolid,'r')
    
    a0.errorbar(10**(sep+0.01),rms_ra,
             rms_ra_error,
             (10**sep - 10**(sep-0.5*step),
              10**(sep+0.5*step) - 10**sep),
             fmt='g.',
             linewidth=0.5,
            label=r'$<\Delta v_\alpha^2>^{1/2}$')

    a0.errorbar(10**(sep-0.01),rms_dec,rms_dec_error,
             (10**sep - 10**(sep-0.5*step),
              10**(sep+0.5*step) - 10**(sep)),             
             fmt='b^',
             linewidth=0.5,
             label=r'$<\Delta v_\delta^2>^{1/2}$')
    
    if showbins:
        for i in range(len(sep)):
            a0.text(10**(sep[i]),1.1*rms_dec[i],str(int(samp[i])))
            
    a0.axvline(x=MONDline, ymin=_ymin, ymax=_ymax,
             linestyle='--',linewidth=0.9)
    if showlegend:
        a0.legend(loc='lower left')
    # SECOND PLOT
        a1.set_ylabel('freq',fontsize=fnt)
        a1.set_yscale('linear')
        a1.set_xlim(np.log10(_xmin),np.log10(_xmax))        
        a1.bar(sep,samp,width=step*0.5,
               color='white',edgecolor='black',tick_label='')
       
    if save:
        plt.savefig('%s.png'%save)
        plt.savefig('%s.jpg'%save)
        plt.savefig('%s.ps'%save)
        plt.savefig('%s.pdf'%save)
    return fig, (a0,a1)

