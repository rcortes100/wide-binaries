import MySQLdb
import numpy as np
import math
from binaries_common import *



sqlconn = MySQLdb.connect('localhost', 'database_username_here', '', 'astro_binaries')
sql = sqlconn.cursor()


for ndist in range(3):
   print 'NDIST ', ndist
   dmu_lim = deltamu_lim[ndist]
   nrows = sql.execute("""
                        select stars.source_id_1, 
                               stars.source_id_2, 
                               local_density.nstars,
                               ra_1,
                               dec_1,
                               parallax_1,
                               parallax_error_1,
                               pmra_1,
                               pmdec_1,
                               dist_1,
                               phot_g_mean_mag_1,
                               ra_2,
                               dec_2,
                               parallax_2,
                               parallax_error_2,
                               pmra_2,
                               pmdec_2,
                               dist_2,
                               phot_g_mean_mag_2,
                               pmra_error_1,
                               pmdec_error_1,
                               pmra_error_2,
                               pmdec_error_2
                        from stars left join local_density on stars.source_id_1=local_density.source_id
                        where sqrt((pmra_1-pmra_2)*(pmra_1-pmra_2)+(pmdec_1-pmdec_2)*(pmdec_1-pmdec_2))<%f and ndist=%d
                        
                       """ % (dmu_lim, ndist))
   sqlres=sql.fetchall()
   for i in xrange(len(sqlres)):
      source_id_1 = sqlres[i][0]
      source_id_2 = sqlres[i][1]
      N_f         = sqlres[i][2]
      ra_1        = sqlres[i][3]
      dec_1       = sqlres[i][4]
      parallax_1  = sqlres[i][5]
      parallax_error_1 = sqlres[i][6]
      pmra_1      = sqlres[i][7]
      pmdec_1     = sqlres[i][8]
      dist_1      = sqlres[i][9]
      phot_g_mean_mag_1 = sqlres[i][10]
      ra_2        = sqlres[i][11]
      dec_2       = sqlres[i][12]
      parallax_2  = sqlres[i][13]
      parallax_error_2 = sqlres[i][14]
      pmra_2      = sqlres[i][15]
      pmdec_2     = sqlres[i][16]
      dist_2      = sqlres[i][17]
      phot_g_mean_mag_2 = sqlres[i][18]
      err_pmra_1 = sqlres[i][19]
      err_pmdec_1 = sqlres[i][20]
      err_pmra_2 = sqlres[i][21]
      err_pmdec_2 = sqlres[i][22]
      
      
      pm_l0 = pmra_1
      pm_b0 = pmdec_1
      pm_l1 = pmra_2
      pm_b1 = pmdec_2
      l0 = ra_1
      b0 = dec_1
      l1 = ra_2
      b1 = dec_2
      d0 = dist_1
      d1 = dist_2
      vr = 0.0
      theta = dsep (ra_1, dec_1, ra_2, dec_2)
      parallax_primary = parallax_1
      err_parallax_primary = parallax_error_1
      err_parallax_secondary = parallax_error_2
      mag_primary = appToAbsMag(phot_g_mean_mag_1, d0)
      mag_companion = appToAbsMag(phot_g_mean_mag_2, d1)
      
      pos_star_1 = sphericalToUnitVect(ra_1, dec_1)*dist_1
      pos_star_2 = sphericalToUnitVect(ra_2, dec_2)*dist_2
      dpos = pos_star_2-pos_star_1
      delta_d = (dpos[0]**2 + dpos[1]**2 + dpos[2]**2)**0.5      
      
      mu_0 = t_mu0[ndist]
      alpha_0 = t_alpha0[ndist]
      delta_mu_1 = t_deltamu1[ndist]
      alpha_1 = t_alpha1[ndist]
      p_p = t_p_p[ndist]
      
      #compute the corrected parallel and ortogonal components
      #of the proper motion
      mod_mu1 = (pm_l1**2+pm_b1**2)**0.5
      mu1_unit = (pm_l1/mod_mu1, pm_b1/mod_mu1)
      mu0_corr = compute_pm_corr (l0, b0, pm_l0, pm_b0, l1, b1, vr, d0, d1)
      mod_mu_par = mu0_corr[0]*mu1_unit[0] + mu0_corr[1]*mu1_unit[1]
      mu_ort = (mu0_corr[0]-mod_mu_par*mu1_unit[0], mu0_corr[1]-mod_mu_par*mu1_unit[1])
      mod_mu_ort = (mu_ort[0]**2 + mu_ort[1]**2)**0.5
      delta_mu_par = np.abs(mod_mu_par-mod_mu1)
      
      #eq (8)
      p_mu_ort_c = np.exp(-((mod_mu_ort/mu_0)**alpha_0))
      
      #eq (9)
      p_mu_par_c = np.exp(-((delta_mu_par/delta_mu_1)**alpha_1))
      
      #eq (10)
      p_deltamu_c = 1 - (1-p_mu_ort_c)*(1-p_mu_par_c)
      
       
      #P(\Delta \mu, \Theta | field) 
      delta_mu = ((pm_l0-pm_l1)**2 + (pm_b0-pm_b1)**2)**0.5
      rho_f = N_f / (np.pi**2 *(deltamu_outer[ndist]**2 - deltamu_inner[ndist]**2)*(theta_outer[ndist]**2-theta_inner[ndist]**2)) #eq 17
      p_deltamu_theta_f = 1 - np.exp(-(np.pi**2)*rho_f*(theta**2)*(delta_mu**2)) #eq 18
      #print ''
      #print 'source_id_1', source_id_1
      #print 'source_id_2', source_id_2
      #print 'ndist', ndist
      #print 'N_f', N_f
      #print 'parallax_primary', parallax_primary
      #print 'parallax_primary**2', parallax_primary**2
      #print 'rho_f', rho_f
      #print 'theta', theta
      #print 'delta_mu', delta_mu
      #print 'theta**2', theta**2
      #print 'delta_mu**2', delta_mu**2
      #print 'exp', -(np.pi**2)*rho_f*(theta**2)*(delta_mu**2)
      
      #P(\Delta d | field)
      lambda_f = N_f/(2*max_diff_dist[ndist]) #eq (19)
      p_deltad_f = 1 - np.exp(-lambda_f*delta_d) #eq (20)
      
      #P(c)
      p_c = (1-p_deltad_f)*(1-p_deltamu_theta_f)*(p_p)/((1-p_deltad_f)*(1-p_deltamu_theta_f)*(p_p) + (p_deltad_f)*(p_deltamu_theta_f)*(1-p_p)) #eq 6
      #p_c = 0.5 #TODO test
      
      #print '(1-p_deltad_f)', (1-p_deltad_f)
      #print '(1-p_deltamu_theta_f)', (1-p_deltamu_theta_f)
      #print 'p_p', p_p
      #print 'den p_c', ((1-p_deltad_f)*(1-p_deltamu_theta_f)*(p_p) + (p_deltad_f)*(p_deltamu_theta_f)*(1-p_p))
      
      #P(r|c) eq (11)
      semi_axis = distance_2points_sph(l0, b0, d0, l1, b1, d1)*206265 #in AU
      err_semiaxis = compute_err_distance_2points_sph(l0, b0, d0, err_parallax_primary/1000, l1, b1, d1, err_parallax_secondary/1000)*206265 #in AU 
      mass_primary = estimate_mass(mag_primary)
      mass_companion = estimate_mass(mag_companion)
      period = ((semi_axis**3)/(mass_primary+mass_companion))**0.5 #this is in years
      e_log_p = 3*err_semiaxis/(2*np.log(10.0)*semi_axis)
      p_r_c = math.erfc ((np.log10(period*365)-4.8)/((2*(2.3**2 + e_log_p**2))**0.5)) #period is needed in days
      
      #p total eq (5)
      p_c_r_deltamu = (p_deltamu_c)*(p_r_c)*(p_c)/(p_deltamu_c*p_r_c*p_c+(1-p_deltamu_c)*(1-p_r_c)*(1-p_c)) 
      #print 'p_deltamu_c', p_deltamu_c
      #print 'p_r_c', p_r_c
      #print 'p_c', p_c
      #print '(p_deltamu_c)*(p_r_c)*(p_c)', (p_deltamu_c)*(p_r_c)*(p_c)
      #print '(p_deltamu_c*p_r_c*p_c+(1-p_deltamu_c)*(1-p_r_c)*(1-p_c))', (p_deltamu_c*p_r_c*p_c+(1-p_deltamu_c)*(1-p_r_c)*(1-p_c))
      if p_c_r_deltamu>0.5:
         print 'p_c_r_deltamu', p_c_r_deltamu, N_f, ndist
      sql.execute ("""
                   insert into probabilities (
                        source_id_1,
                        source_id_2,
                        ra_1,
                        dec_1,
                        d1,
                        mag_g_1,
                        parallax_1,
                        err_parallax_1,
                        pmra_1,
                        pmdec_1,
                        err_pmra_1,
                        err_pm_dec_1,
                        ra_2,
                        dec_2,
                        d2,
                        mag_g_2,
                        parallax_2,
                        err_parallax_2,
                        pmra_2,
                        pmdec_2,
                        err_pmra_2,
                        err_pm_dec_2,
                        angsep,
                        sep,
                        prob_companion) 
                   values (%d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f)
                   """ % (source_id_1, 
                          source_id_2,
                          ra_1,
                          dec_1,
                          dist_1,
                          phot_g_mean_mag_1,
                          parallax_1,
                          parallax_error_1,
                          pmra_1,
                          pmdec_1,
                          err_pmra_1,
                          err_pmdec_1,
                          ra_2,
                          dec_2,
                          dist_2,
                          phot_g_mean_mag_2,
                          parallax_2,
                          parallax_error_2,
                          pmra_2,
                          pmdec_2,
                          err_pmra_2,
                          err_pmdec_2,
                          theta,
                          semi_axis/206265,
                          p_c_r_deltamu))
   sqlconn.commit()
      
