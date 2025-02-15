import requests
import time
import sys

login_url = 'https://gea.esac.esa.int/tap-server/login'
logout_url = 'https://gea.esac.esa.int/tap-server/logout'
query_url = 'https://gea.esac.esa.int/tap-server/tap/async'

#-------------------------------------
#Create job


p = {'username':'GAIA_ARCHIVE_USERNAME_HERE', 'password':'GAIA_ARCHIVE_PASS_HERE'}
s = requests.Session()

r = s.post(login_url, params = p) 
c = r.cookies

print >> sys.stderr, 'login done, cookies: ', c

#make query

index = int(sys.argv[1])
distance0 = (0.0, 25.0, 50.0)
distance1 = (25.0, 50.0, 100.0)
inner_radius = (0.01, 0.01, 0.005)
outer_radius = (20.0, 10.0, 5.0)
max_diff_dist = (20.0, 15.0, 20.0)
diff_pm_outer = (100.0, 75.0, 50.0)

sleep_time = 10

print >> sys.stderr, 'index of distance: ', index

query = """select t1.source_id as source_id_1, t1.ra as ra_1, t1.dec as dec_1, t1.parallax as parallax_1, t1.parallax_error as parallax_error_1, t1.pmra as pmra_1, t1.pmdec as pmdec_1, t1.pmra_error as pmra_error_1, t1.pmdec_error as pmdec_error_1, 1000.0/t1.parallax as dist_1, t1.phot_g_mean_mag as phot_g_mean_mag_1, t2.source_id as source_id_2, t2.ra as ra_2, t2.dec as dec_2, t2.parallax as parallax_2, t2.parallax_error as parallax_error_2, t2.pmra as pmra_2, t2.pmdec as pmdec_2, t2.pmra_error as pmra_error_2, t2.pmdec_error as pmdec_error_2, 1000.0/t2.parallax as dist_2, t2.phot_g_mean_mag as phot_g_mean_mag_2, abs(1000.0/t1.parallax-1000.0/t2.parallax) as ddist, sqrt(power(t1.pmra-t2.pmra, 2)+power(t1.pmdec-t2.pmdec,2)) as dpm 
          from gaiadr1.tgas_source as t1 inner join gaiadr1.tgas_source as t2 on
               t1.source_id <> t2.source_id and  
               abs(t1.parallax_error/t1.parallax)<0.1 and abs(t2.parallax_error/t2.parallax)<0.1 and
               0=CONTAINS(POINT('ICRS', t2.ra, t2.dec), CIRCLE('ICRS',t1.ra,t1.dec, %f)) and
               1=contains(POINT('ICRS', t2.ra, t2.dec), CIRCLE('ICRS',t1.ra,t1.dec, %f)) and 
               abs(1000.0/t1.parallax-1000.0/t2.parallax)<%f and 
               sqrt(power(t1.pmra-t2.pmra, 2)+power(t1.pmdec-t2.pmdec,2))<%f and 
               1000.0/t1.parallax between %f and %f                 
               where 
               t1.ra between %f and %f"""


print '#RA, JOBID'
print '#index='+str(index)
for ra in range (0, 360):
   print >> sys.stderr, 'Requesting RA ', ra
   p = {"REQUEST": "doQuery", \
	   "LANG":    "ADQL", \
      "FORMAT":  "votable", \
      "PHASE":  "RUN", \
      "JOBNAME": "binaries ra "+str(ra)+" index "+str(index), \
      "jobdescription": "una descripcion", \
      "QUERY": query%(inner_radius[index], outer_radius[index], max_diff_dist[index], diff_pm_outer[index], distance0[index], distance1[index], ra, ra+1)
	}	
   r = s.post(query_url, cookies = c, params = p)   
   location = r.url   
   jobid = location[location.rfind('/')+1:]
   print str(ra)+', '+jobid
   print >> sys.stderr, 'Done! ', jobid
   print >> sys.stderr, 'Sleeping for '+str(sleep_time)+' sec', time.time()
   time.sleep(sleep_time)


r = s.post(logout_url, cookies = c, params = p) 
print 'done logout, ', r


