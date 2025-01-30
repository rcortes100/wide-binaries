import numpy as np


#index = int(sys.argv[1])

distance0 = (0.0, 25.0, 50.0)
distance1 = (25.0, 50.0, 100.0)
inner_radius = (0.01, 0.01, 0.005)
outer_radius = (20.0, 10.0, 5.0)
max_diff_dist = (20.0, 15.0, 20.0)
diff_pm_outer = (100.0, 75.0, 50.0)
deltamu_lim = np.array([35, 25, 10])
deltamu_inner = np.array([45, 35, 20])
deltamu_outer = np.array([100, 75, 50])
theta_inner = np.array([0.01, 0.01, 0.005])
theta_outer = np.array([20, 10, 5])

t_mu0 = (10.0, 3.44, 2.09)
t_alpha0 = (0.97, 0.93, 1.06)
t_deltamu1 = (7.13, 3.39, 2.28)
t_alpha1 = (0.97, 0.78, 1.15)
t_p_p = (0.5, 0.2, 0.07)

def sphericalToUnitVect(_a, _b):
   a = np.radians(_a)
   b = np.radians(_b)
   cosb = np.cos(b)
   return np.array([np.cos(a) * cosb,
                    np.sin(a) * cosb,
                    np.sin(b)])
   

#angular separation 
#parameters and result in degrees
def dsep (a1, b1, a2, b2):
   v1 = sphericalToUnitVect(a1, b1)
   v2 = sphericalToUnitVect(a2, b2)
   v1xv2 = np.cross (v1, v2)
   s = (v1xv2[0]**2 + v1xv2[1]**2 + v1xv2[2]**2)**0.5
   c = np.dot (v1, v2)
   if s != 0.0:
      return np.degrees(np.arctan2(s, c))
   else:
      return 0.0
   
def distance_2points_sph(a1, b1, d1, a2, b2, d2):
   v1 = sphericalToUnitVect(a1, b1)*d1
   v2 = sphericalToUnitVect(a2, b2)*d2
   dv = v2-v1
   return (dv[0]**2 + dv[1]**2 + dv[2]**2)**0.5

def compute_err_distance_2points_sph(a1, b1, d1, err_parallax_1, a2, b2, d2, err_parallax_2):
   v1 = sphericalToUnitVect(a1, b1)*d1
   v2 = sphericalToUnitVect(a2, b2)*d2
   r = distance_2points_sph(a1, b1, d1, a2, b2, d2)
   xp=v1[0]
   yp=v1[1]
   zp=v1[2]
   xc=v2[0]
   yc=v2[1]
   zc=v2[2]   
   xcomp = ((err_parallax_1*xp*xp)**2 + (err_parallax_2*xc*xc)**2)*((xp-xc)**2)
   ycomp = ((err_parallax_1*yp*yp)**2 + (err_parallax_2*yc*yc)**2)*((yp-yc)**2)
   zcomp = ((err_parallax_1*zp*zp)**2 + (err_parallax_2*zc*zc)**2)*((zp-zc)**2)
   return ((xcomp + ycomp + zcomp)**0.5)/r

#compute spatial velocity from proper motion
#units: pm_l, pm_b: mas/yr
#       vr: km/s
#        l: degrees
#        b: degrees
#        d: pc
#return: (vx, vy, vz) in km/s
def spatialVelFromPm (_pm_l, _pm_b, _vr, _l, _b, d):
   pm_l = _pm_l/1000.0
   pm_b = _pm_b/1000.0
   vr = _vr/4.74 #now vr is in AU/yr, 1 AU/yr=4.74km/s
   l=np.radians(_l)
   b=np.radians(_b)
   vx = vr*np.cos(b)*np.cos(l) - d*pm_l*np.sin(l) - d*pm_b*np.sin(b)*np.cos(l)
   vy = vr*np.cos(b)*np.sin(l) + d*pm_l*np.cos(l) - d*pm_b*np.sin(b)*np.sin(l)
   vz = vr*np.sin(b) + d*pm_b*np.cos(b)
   return vx*4.74, vy*4.74, vz*4.74

#compute proper motion from space velocity
#units: 
#        vx, vy, vz: km/s
#        l, b: degrees
#        d: pc
#return:
#        pm_l, pm_b: mas/yr
#        vr: km/s
def pmFromSpatialVel (_vx, _vy, _vz, _l, _b, d):
   vx = _vx/4.74
   vy = _vy/4.74
   vz = _vz/4.74
   l = np.radians(_l)
   b = np.radians(_b)
   pm_l = (-vx*np.sin(l)+vy*np.cos(l))/d
   pm_b = ((-vx*np.cos(l)-vy*np.sin(l))*np.sin(b)+vz*np.cos(b))/d
   vr = (vx*np.cos(l)+vy*np.sin(l))*np.cos(b)+vz*np.sin(b)
   return pm_l*1000, pm_b*1000, vr*4.74

def compute_pm_corr (l0, b0, pm_l0, pm_b0, l1, b1, vr, d0, d1):
   v = spatialVelFromPm(pm_l0, pm_b0, vr, l0, b0, d0)
   projected_pm = pmFromSpatialVel(v[0], v[1], v[2], l1, b1, d1)
   return (projected_pm[0]), projected_pm[1]

def compute_pm_diff (l0, b0, pm_l0, pm_b0, l1, b1, pm_l1, pm_b1, vr, d0, d1):
   v = spatialVelFromPm(pm_l0, pm_b0, vr, l0, b0, d0)
   projected_pm = pmFromSpatialVel(v[0], v[1], v[2], l1, b1, d1)
   return (projected_pm[0]-pm_l1), projected_pm[1]-pm_b1

#in solar masses
def estimate_mass (mag):
   mag_sun = 4.83
   #first, try 0.43Mo < M < 2Mo
   #a=4 in this case
   a=4
   mass = 10**(-0.4/a*(mag-mag_sun))
   if mass>0.43 and mass<=2:
      return mass
   else:
      a=3.5
      mass = 10**(-0.4/a*(mag-mag_sun))
      return mass

def appToAbsMag (m, d):
   return m - 5*(np.log10(d)-1)

   
