import sys
import MySQLdb
from astropy.io import fits
import numpy as np
import time
from astropy.io.votable import parse_single_table
import glob
import warnings
warnings.filterwarnings("ignore")

keyword_dict = {'DEC':'DECLINATION', 'RUN':'RUNID'}
types_dict = {bool:'bool', float:'double', int:'int', str:'varchar(255)', np.float32:'float'}

#Traduce usando el diccionario
def translate_key(k):
   if k in keyword_dict:
      return keyword_dict[k]
   else:
      return k.replace('-', '_')   

# Convierte los tipos de dato en SQL
def sql_quote(v):
   if type(v) is int:
      return str(v)
   elif type(v) is float:
      return str(v)
   elif type(v) is bool:
      return str(v)
   elif type(v) is np.float32:
      return str(v)   
   elif type(v) is str:
      return "'"+v+"'"
   else:
      print v, type(v)
      return None      

# Crea tabla
def create_tables (data):
   sql.execute ("drop table stars")
   query="""create table stars (
     id int not null auto_increment primary key,
     ndist int not null
     """ 

   for field in data.dtype.names:
      t=data.dtype[field]
      if t== np.dtype('int64'):
         query += ', '+field+' bigint not null'
      else:
         query += ', '+field+' float not null'
   query += ")"
   sql.execute(query)

# Aquí comienza el script

if __name__=='__main__':

   sqlconn = MySQLdb.connect('localhost', 'DATABASE_USERNAME', '', 'astro_binaries')
   sql = sqlconn.cursor()

   first_iter = True

   for dist in range(3):
      nfiles = glob.glob ('dist'+str(dist)+'/*.xml')
      for nfile in nfiles:
         print 'Parsing ', nfile
         table = parse_single_table(nfile)
         data = table.array
         if first_iter:
            create_tables(data)
            first_iter = False
            nfields = data.dtype.names
         for i in xrange(len(data)):
            query = 'insert into stars (ndist, '+ (', '.join(nfields)) +') values ('+str(dist)+',' + (', '.join([str(data[i][v]) for v in nfields]))+')'
            #print data[i], query
            sql.execute(query)
         sqlconn.commit()      
