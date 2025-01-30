#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import requests
import time
import sys
import math

# PAGE ADRESSES
login_url = 'https://gea.esac.esa.int/tap-server/login'
logout_url = 'https://gea.esac.esa.int/tap-server/logout'
query_url = 'https://gea.esac.esa.int/tap-server/tap/async'

#TIME BETWEEN QUERIES
sleep_time = 3

phase = int(sys.argv[1])

#-------------------------------------
#Create job

# USER CREDENTIALS

if (phase==2 or phase==4):
   usname='rscarpa'
   pwd='X1mondo+='
else:
   usname='rcortesm'
   pwd='A$tronomyDomin3.'

print usname, phase
raw_input()

credentials = {'username': usname, 'password': pwd}
session = requests.Session()

request = session.post(login_url, params = credentials) 
cooks = request.cookies

print >> sys.stderr, 'Login done, cookies:\n', cooks

#Output of the job id, and the address to retrieve the data

# MAKE THE QUERY

index = 0;
if len(sys.argv) > 2:
   index = int(sys.argv[2]) #INDEX MEANS HOW DISTANT STARS WILL BE

distance0     = (0.0,     25.0, 50.0)    #Paralaje minima a la estrella en mas
distance1     = (999999., 50.0, 100.0)   #Paralaje maxima a la estrella en mas
inner_radius  = (0.01,    0.01, 0.005)   #Radio interno en grados
outer_radius  = (20.0,    10.0, 5.0)     #Radio externo en grados
max_diff_dist = (20.0,    15.0, 20.0)  #Separación maxima en componente radial
diff_pm_outer = (100.0,   75.0, 50.0)    #Max. vel. relativa de mov. prop.(mas/year)

print >> sys.stderr, 'Index of distance: ', index

#The query takes positions and velocities (including radial)
#We want to take the query with ra and dec centered
#
query = \
"""SELECT
  t1.source_id             AS source_id_1,
  t2.source_id             AS source_id_2,
  t1.ra                    AS ra_1,
  t1.ra_error              AS ra_err_1,
  t1.dec                   AS dec_1,
  t1.dec_error             AS dec_err_1,
  t2.ra                    AS ra_2,
  t2.ra_error              AS ra_err_2,
  t2.dec                   AS dec_2,
  t2.dec_error             AS dec_err_2,
  t1.parallax              AS parallax_1,
  t1.parallax_error        AS parallax_err_1,
  t2.parallax              AS parallax_2,
  t2.parallax_error        AS parallax_err_2,
  t1.pmra                  AS pmra_1,
  t1.pmra_error            AS pmra_err_1,
  t2.pmra                  AS pmra_2,
  t2.pmra_error            AS pmra_err_2,
  t1.pmdec                 AS pmdec_1,
  t1.pmdec_error           AS pmdec_err_1,
  t2.pmdec                 AS pmdec_2,
  t2.pmdec_error           AS pmdec_err_2,
  t1.radial_velocity       AS rvel_1,
  t1.radial_velocity_error AS rvel_err_1,
  t2.radial_velocity       AS rvel_2,
  t2.radial_velocity_error AS rvel_err_2,
  t1.phot_g_mean_mag       AS G_1,
  t2.phot_g_mean_mag       AS G_2,
  1000*ABS(t1.parallax-t2.parallax)/(t1.parallax*t2.parallax) AS dradial,
  DISTANCE(POINT('ICRS',t1.ra,t1.dec),POINT('ICRS',t2.ra,t2.dec)) AS dtang,
  1000*SQRT(POWER(1/t1.parallax-1/t2.parallax,2)+(COS(t1.dec*PI()/180)*COS(t2.dec*PI()/180)*POWER(t1.ra-t2.ra,2)+POWER(t1.dec-t2.dec,2))*POWER(PI()/180,2)/(t1.parallax*t2.parallax)) AS separation
	FROM
	gaiadr2.gaia_source AS t1
	INNER JOIN
	gaiadr2.gaia_source AS t2
	ON
	t1.source_id < t2.source_id AND
	t1.parallax_over_error > %f AND
	t2.parallax_over_error > %f AND
	POWER(%f,2)*t1.parallax*t2.parallax > POWER(17.45*DISTANCE(POINT('ICRS',t1.ra,t1.dec),POINT	('ICRS',t2.ra,t2.dec)),2) AND
	%f*t1.parallax*t2.parallax > 1000*ABS(t1.parallax-t2.parallax)
	WHERE	
	CONTAINS(POINT('ICRS',t1.ra,t1.dec),CIRCLE('ICRS',%f,%f,%f))=1 AND
	CONTAINS(POINT('ICRS',t2.ra,t2.dec),CIRCLE('ICRS',%f,%f,%f))=1 AND
	(t1.parallax BETWEEN %f AND %f) AND
	(t2.parallax BETWEEN %f AND %f)
"""
# parallax_over_error_minimum = 20
# maximum_separation R=1pc
# R
# Center of circle and radius: ra, dec, radius=0.5 degree
# min_parallax, max_parallax

print '#RA, JOBID'
print '#index='+str(index)

inicio = int(raw_input('Select initial range: '))

if phase==2:
   declinations = [x+0.5 for x in range(inicio,90)]
elif phase==3:
   declinations = range(-89,inicio+1)
   declinations.reverse()
elif phase==4:
   declinations = [-x-0.5 for x in range(inicio,90)]
else:
   declinations = range(inicio,90)

print declinations
raw_input()

for dec in declinations:
   #crear archivo
   file_name = './sessions/session'+time.strftime('%s')+'dec_'+str(dec)+'.txt'
   if (phase==2 or phase==4):
         file_name = './sessions/session2_'+time.strftime('%s')+'dec_'+str(dec)+'.txt'
   file_list = open(file_name,'w')
   ra_points = int(360*math.cos(math.radians(dec)))+1 #Distribucion de puntos
   
   for ra in range(ra_points):
      if (phase==2 or phase==4):
         ra = (ra+0.5)/math.cos(math.radians(dec)) #Distribucion de puntos
      else:
         ra = ra/math.cos(math.radians(dec)) #Distribucion de puntos
      print >> sys.stderr, 'Requesting RA ', ra, 'dec', dec
      
      
#      print query%(20,1,1,ra,dec,0.5,ra,dec,0.5,distance0[index],distance1[index],distance0[index],distance1[index])
#      raw_input()
      jobname = "bins_ra_"+str(ra)+"_dec_"+str(dec)
      credentials = {"REQUEST": "doQuery", \
                        "LANG": "ADQL", \
                      "FORMAT": "votable", \
                       "PHASE": "RUN", \
                     "JOBNAME": jobname, \
              "JOBDESCRIPTION": "Candidatos a estrellas binarias", \
                       "QUERY": query%(20,20,1,1,ra,dec,0.5,ra,dec,0.5,distance0[index],distance1[index],distance0[index],distance1[index])
	}	
      request = session.post(query_url, cookies = cooks, params = credentials)
      location = request.url
      #Escribir en el archivo
      file_list.write(jobname + '\t' + location + '\n')
      print 'location: '+str(location)
      jobid = location[location.rfind('/')+1:]
      print str(ra)+', '+jobid
      print >> sys.stderr, 'Done! ', jobid
      print >> sys.stderr, 'Sleeping for '+str(sleep_time)+' sec'
      time.sleep(sleep_time)
      
   # Cierro el archivo para abrirlo de nuevo más adelante
   file_list.close()     
request = session.post(logout_url, cookies = cooks, params = credentials) 
print 'done logout, ', request
