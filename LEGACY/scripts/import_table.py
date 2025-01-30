#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parsec = 206265         #Astronomical Units
AU     = 1.49597871E11  #Meters

def DIST_ANGdDIST_ANG_rad(ra1,ra1_err,ra2,ra2_err,dec1,dec1_err,dec2,dec2_err):
   a    = np.abs(ra2-ra1)
   b    = np.abs((ra2+180)%360 - (ra1+180)%360) #minimum of a and b
   dra  = 0.5*(a + b - np.abs(a - b))
   ddec = np.abs(dec2-dec1)
   return dra*(ra2_err+ra1_err) + (dec2_err + dec1_err)*(ddec + 0.5*dra**2)

def sep(par1,par2,distang):
   dr = 1000*(1/par1-1/par2)
   return np.sqrt(dr**2 + (1000*np.radians(distang))**2/(par1*par2))

def dsep (r1,r1_err,r2,r2_err,
          ra1,ra1_err,ra2,ra2_err,
          dec1,dec1_err,dec2,dec2_err,sep,distang):
   # Primero hay que convertir a radianes. Coordenadas en grados, y errores en mas
   ra1  = np.radians(ra1)
   ra2  = np.radians(ra2)
   dec1 = np.radians(dec1)
   dec2 = np.radians(dec2)
   ra1_err  = np.radians(ra1_err/3.6E6)
   ra2_err  = np.radians(ra2_err/3.6E6)
   dec1_err = np.radians(dec1_err/3.6E6)
   dec2_err = np.radians(dec2_err/3.6E6)
   distang  = np.radians(distang)

   return abs(r2-r1)*(r1_err+r2_err)/sep + \
      (r1*r2/sep)*(DIST_ANGdDIST_ANG_rad(ra1,ra1_err,ra2,ra2_err,dec1,dec1_err,dec2,dec2_err)+ \
                   0.5*(r1_err/r1 + r2_err/r2)*distang**2)

def veltan(pmr1,pmd1,pmr2,pmd2,r1,r2):
   #Convertir de mas/yr a rad/s
   pmr1 = np.radians(pmr1/3.6E6)/(365*24*3600)
   pmr2 = np.radians(pmr2/3.6E6)/(365*24*3600)
   pmd1 = np.radians(pmd1/3.6E6)/(365*24*3600)
   pmd2 = np.radians(pmd2/3.6E6)/(365*24*3600)

   #Convertir de pc a km
   r1 = r1*parsec*AU/1000
   r2 = r2*parsec*AU/1000
   
   return np.sqrt((r2*pmr2-r1*pmr1)**2+(r2*pmd2-r1*pmd1)**2)

def dveltan(pmr1,pmr1_err,pmd1,pmd1_err,pmr2,pmr2_err,pmd2,pmd2_err,
            r1,r1_err,r2,r2_err):
   
   #Convertir de mas/yr a rad/s
   pmr1     = np.radians(pmr1/3.6E6)/(365*24*3600)
   pmr1_err = np.radians(pmr1_err/3.6E6)/(365*24*3600)
   pmr2     = np.radians(pmr2/3.6E6)/(365*24*3600)
   pmr2_err = np.radians(pmr2_err/3.6E6)/(365*24*3600)
   pmd1     = np.radians(pmd1/3.6E6)/(365*24*3600)
   pmd1_err = np.radians(pmd1_err/3.6E6)/(365*24*3600)
   pmd2     = np.radians(pmd2/3.6E6)/(365*24*3600)
   pmd2_err = np.radians(pmd2_err/3.6E6)/(365*24*3600)

   #Convertir de pc a km
   r1     = r1*parsec*AU/1000
   r1_err = r1_err*parsec*AU/1000
   r2     = r2*parsec*AU/1000
   r2_err = r2_err*parsec*AU/1000

   vtan = np.sqrt((r2*pmr2-r1*pmr1)**2+(r2*pmd2-r1*pmd1)**2)
   
   return (np.abs(r2*pmr2-r1*pmr1)*np.abs(r2*pmr2_err + pmr2*r2_err +
                                          r1*pmr1_err + pmr1*r1_err)+
           np.abs(r2*pmd2-r1*pmd1)*np.abs(r2*pmd2_err + pmd2*r2_err +
                                          r1*pmd1_err + pmd1*r1_err))/vtan

if __name__=='__main__':

# Remember: ra and dec are in degrees
#           ra_error and dec_error are in mas


   # Ver si se ha seleccionado el hemisferio norte, el sur, o todo

   north = False
   south = False
   
   for a in sys.argv:
      if a.startswith('-'):
         if a.count('n'):
            north = True
         if a.count('s'):
            south = True

   prefix = []
   if north:
      prefix.append('~/cdn.gea.esac.esa.int/results_wb/north_hemisphere')
   if south:
      prefix.append('~/cdn.gea.esac.esa.int/results_wb/south_hemisphere')

   exit_file = raw_input('Please insert the exit file name: ')

   tab = pd.DataFrame()

   for pre in prefix:
      if pre.count('north'):
#         declinations = [0.5*x for x in range(179)]
         declinations = [str(0.5*x).zfill(4) for x in range(179)]
      if pre.count('south'):
#         declinations = [-(0.5*x)-0.5 for x in range(179)]
         declinations = [str(-0.5*x-0.5).zfill(5) for x in range(179)]

      for dec in declinations:
         filename = '%s/dec%s/bins_dec_%s.tmp'%(pre,dec,dec)
         t1 = pd.read_csv(filename,dtype={'source_id_1':object,'source_id_2':object})

         #Corregir la separacion que está mal calculada
         t1 = t1.drop('separation',axis=1)

         #Limpiar calidad de datos
         t1 = t1[t1.parallax_1/t1.parallax_err_1 > 20]
         t1 = t1[t1.parallax_2/t1.parallax_err_2 > 20]

         #Crear nuevas columnas
         t1['r1'] = 1000/t1.parallax_1
         t1['r1_err'] = t1.r1*t1.parallax_err_1/t1.parallax_1
         t1['r2'] = 1000/t1.parallax_2
         t1['r2_err'] = t1.r2*t1.parallax_err_2/t1.parallax_2
         t1['separation']  = sep(t1.parallax_1,t1.parallax_2,t1.dtang)
         t1['dseparation'] = dsep (t1.r1,t1.r1_err,
                                   t1.r2,t1.r2_err,
                                   t1.ra_1,t1.ra_err_1,
                                   t1.ra_2,t1.ra_err_2,
                                   t1.dec_1,t1.dec_err_1,
                                   t1.dec_2,t1.dec_err_2,
                                   t1.separation,t1.dtang)
         t1['veltan']  = veltan(t1.pmra_1,t1.pmdec_1,t1.pmra_2,t1.pmdec_2,
                                t1.r1,t1.r2)
         t1['dveltan'] = dveltan(t1.pmra_1,t1.pmra_err_1,
                                 t1.pmdec_1,t1.pmdec_err_1,
                                 t1.pmra_2,t1.pmra_err_2,
                                 t1.pmdec_2,t1.pmdec_err_2,
                                 t1.r1,t1.r1_err,t1.r2,t1.r2_err)
         #Depurar por calidad de separacion radial y velocidad tangencial
         
         t1 = t1[t1.dseparation < 1]
         t1 = t1[t1.dveltan < 1]

         #Agregar tabla

         tab = tab.append(t1,ignore_index=True)
         print "File %s processed"%filename

      tab.to_csv(exit_file,index=False)
   print "Ready!"
   
"""
   t1['dsep'] = sepdsep(tabla.parallax_1,tabla.parallax_err_1,
                           tabla.parallax_2,tabla.parallax_err_2,
                           np.deg2rad(abs(tabla.ra_2-tabla.ra_1)),
                           np.deg2rad(tabla.ra_err_1+tabla.ra_err_2)/3.6E6,
                           np.deg2rad(tabla.dec_1),np.deg2rad(tabla.dec_err_1)/3.6E6,
                           np.deg2rad(tabla.dec_2),np.deg2rad(tabla.dec_err_2)/3.6E6)/tabla.separation
   
"""
