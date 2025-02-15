#!/usr/bin/python2.7

# -*- coding: utf-8 -*-
import getpass as gp
import requests
import time
import sys

login_url       = 'https://gea.esac.esa.int/tap-server/login'
logout_url      = 'https://gea.esac.esa.int/tap-server/logout'
download_suffix = '/results/result'

#-------------------------------------
#Create job

# USER CREDENTIALS
usname='rscarpa'
pwd='X1mondo+='

#usname = 'rcortesm'
#pwd    = 'A$tronomyDomin3.'

#usname = raw_input('Insert username:')
#pwd    = gp.getpass()

p = {'username': usname, 'password': pwd}
s = requests.Session()

r = s.post(login_url, params = p) 
c = r.cookies

print >> sys.stderr, 'login done, cookies: \n', c

# Enter file with the job IDs
nfile = sys.argv[1]
#nfile = './sessions/session1526381616dec_0.txt'

#Destination directory
dirname=sys.argv[2]
#dirname = './results'

# This must be performed for every value of declination

f     = open(nfile)      #Open file
fdat  = f.read()         #Read all the contents
lines = fdat.split('\n') #Split by lines

for line in lines:
   if len(line)>0 and line[0]!='#': #Ommit commented or empty lines
      split = line.split()          #Split line
   if len(split)==2:
      downlink     = split[1] + download_suffix #Take the second column as download page
      downlink_csv = split[1] + download_suffix + '?format=csv'
      print >> sys.stderr, 'Downloading result for ' + split[0] + \
         ' from\n' + downlink
      r     = s.get(downlink, cookies = c)     #Fetch the table from the url
      r_csv = s.get(downlink_csv, cookies = c) #Fetch the table from the url      
      fres = open (dirname+'/'+split[0]+'.xml', 'wb') #Open in destination directory
      fres_csv = open (dirname+'/'+split[0]+'.csv', 'w')
      data     = r.content
      data_csv = r_csv.content
      print >> sys.stderr, '... done, ', len(data), ' bytes', len(data_csv), 'bytes'
      fres.write(data)
      fres_csv.write(data_csv)
      fres.close()
      fres_csv.close()

r = s.post(logout_url, cookies = c, params = p) 
print 'done logout, ', r
