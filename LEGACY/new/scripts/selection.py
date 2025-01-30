#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import import_table as imp

cumulos = [[  5.66541667e+01 ,  2.41569444e+01],
           [  6.71333333e+01 ,  1.59455556e+01],
           [  1.30008333e+02 ,  1.96822222e+01],
           [  8.37458333e+01 , -4.60194444e+00],
           [  8.43791667e+01 , -1.18333333e-01],
           [  1.30070833e+02 ,  1.97372222e+01],
           [  1.30262500e+02 , -5.28647222e+01],
           [  2.71025000e+02 , -2.36788889e+01],
           [  2.70629167e+02 , -2.29927778e+01],
           [  2.71025000e+02 , -2.14722222e+01],
           [  2.46133333e+02 , -2.36163889e+01]]

if __name__=='__main__':
   archivo = raw_input('nombre del archivo: ')
   t = pd.read_csv(archivo,dtype={'source_id_1': object, 'source_id_2':object})
   h = pd.read_csv('../hipparcos_best_neighbours.csv',dtype={'source_id': object, 'original_ext_source_id': object})
   t = t.drop_duplicates(subset='source_id_1', keep=False)
   t = t[(t.dseparation/t.separation < 0.2) &
         (t.dveltan/t.veltan < 0.2) &
         (t.g_1 < 16) & (t.g_2 < 16)]

   #Tabla Shaya
   shaya = pd.read_csv('../shaya.csv',dtype={'Ind': object, 'Pri':object, 'Cmp':object})
   
   #Hacer una seleccion
#   t[(t.dseparation/t.separation < 0.2) & (t.dveltan/t.veltan < 0.2)].plot.scatter(x='separation',y='veltan',s=0.1)

"""   
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
         t1 = pd.read_csv(filename)

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

      tab.to_csv(exit_file)
   print "Ready!"
"""
