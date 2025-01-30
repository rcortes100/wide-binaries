import requests
import time
import sys

login_url = 'https://gea.esac.esa.int/tap-server/login'
logout_url = 'https://gea.esac.esa.int/tap-server/logout'
download_url = 'https://gea.esac.esa.int/tap-server/tap/async/%s/results/result'

#-------------------------------------
#Create job

dirname=sys.argv[2]

p = {'username':'GAIA_ARCHIVE_USERNAME_HERE', 'password':'GAIA_ARCHIVE_PASS_HERE'}
s = requests.Session()

r = s.post(login_url, params = p) 
c = r.cookies

print >> sys.stderr, 'login done, cookies: ', c

nfile = sys.argv[1]
f=open(nfile)
fdat=f.read()
lines = fdat.split('\n')
for line in lines:
   if len(line)>0 and line[0]!='#':
    split = line.split(', ')
    if len(split)==2:
      jobid=split[1]
      print >> sys.stderr, 'Downloading result for RA '+split[0]+", jobid: "+jobid
      url = download_url %(jobid)
      r = s.get(url, cookies = c)   
      fres = open (dirname+'/'+jobid+'.xml', 'wb')
      data = r.content
      print >> sys.stderr, '... done, ', len(data), ' bytes'
      fres.write (data)
      fres.close()


r = s.post(logout_url, cookies = c, params = p) 
print 'done logout, ', r

