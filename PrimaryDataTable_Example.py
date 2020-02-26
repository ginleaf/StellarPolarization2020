#!/usr/bin/env python
# coding: utf-8

# In[18]:


import h5py
import matplotlib.pyplot as plt
import numpy as np
import polarization_functions
import coordinate_conversions
import CDS
from importlib import reload
from astropy.coordinates import SkyCoord
import astropy.units as unit
from astropy.io import ascii
from astropy.table import Column


# # This example will walk us through getting from a table as given in a paper to a table in the format we need for this first step (the Primary Data Table)

# # Example1: Soam+2018

# ## No electronic format available: I copied the measurements to a csv file containing only the original columns in the paper

# ## Let's read in the file using astropy ascii

# In[174]:


data = ascii.read("Soam2018_table4.csv")


# In[129]:


data


# # This file so far has only what info was given in the paper. Note that I've named the relevant columns according to our convention.
# # We'll keep it simple and create a new dictionary that will hold the information of the columns and data we need. 
# ### (If you're familiar with astropy tables, feel free to use them instead.)

# In[130]:


# This is the number of sources in the table
Nsources = len(data['rawstarID'])

# Column names that we need in our final table - not those that need any amount of processing other than adding nans
colnames = ['starID','Name','rawstarID','RefID','Instrument','JD','FilterID',             'rawRA','rawDec','rawQ','rawE_Q','rawU','rawE_U','rawP','rawE_P',             'rawP_D','rawEVPA','rawE_EVPA','IntrPol','RA', 'Dec', 'q', 'e_q',             'u', 'e_u', 'p', 'e_p', 'evpa', 'e_evpa']

# Dictionary that will represent our table
my_table = {}

# Loop over all column names, add NaN columns where no info
for ii, cname in enumerate(colnames):
    if cname not in data.keys():
        if cname == 'IntrPol':
            NewArray = np.zeros(Nsources, dtype = int)
        elif cname == 'rawRA' or cname == 'rawDec':
            NewArray = np.array(['nan' for jj in range(Nsources)], dtype = 'U')
        else: 
            NewArray  = np.zeros(Nsources)+np.nan
        my_table[cname] = NewArray
    else:
        my_table[cname] = data[cname].data
        
        
# I'll add columns to the original table so we can visualize what happened        
for ii, cname in enumerate(colnames):
    if cname not in data.keys():
        NewColumn = Column(np.zeros(Nsources)+np.nan)
        data.add_column(NewColumn,index = ii+1, name = cname)
    else:
        continue


# In[131]:


data


# # Now let's create the minimally processed information

# ## Add metadata

# In[132]:


####JD

JD = 2456566.5
my_table['JD'] = np.zeros(Nsources)+JD

###Instrument

# I'm sleepy so this is stupid code :(
Instr = 'AIMPOL'
Il = []
for ii in range(Nsources):
    Il.append(Instr)
my_table['Instrument'] = np.array(Il, dtype = 'U')

### RefID
RefID = 14
my_table['RefID'] = np.zeros(Nsources, dtype = int)+RefID

## FilterID
FiltID = 0
my_table['FilterID'] = np.zeros(Nsources, dtype = int)+FiltID


# ## Convert rawstarID to string, if not already

# In[133]:


if type(my_table['rawstarID'][0])!=np.string_ or type(my_table['rawstarID'][0])!=np.str_:
    my_table['rawstarID'] = my_table['rawstarID'].astype('U')


# # Convert P to fraction

# In[134]:


my_table['p'] = my_table['rawP']/100.
my_table['e_p'] = my_table['rawE_P']/100.


# ## Rename EVPA, keep as is for this paper

# In[148]:


my_table['evpa'] = my_table['rawEVPA']
my_table['e_evpa'] = my_table['rawE_EVPA']


# # Get q, u from p, EVPA

# In[149]:


# remember to convert EVPA to radians when passing it into the function
my_table['q'], my_table['u'], my_table['e_q'], my_table['e_u'] =                 polarization_functions.qu_from_p_EVPA(my_table['p'], np.radians(my_table['rawEVPA']),                                                      my_table['e_p'], np.radians(my_table['rawE_EVPA']))


# # No coordinates were specified. Query Simbad to get RA, DEC

# In[141]:


from astroquery.simbad import Simbad
my_table['RA'] = np.zeros(len(data['Name'].data))+np.nan
my_table['Dec'] = np.zeros(len(data['Name'].data))+np.nan

for ii in range(len(data['Name'].data)):
    result_table = Simbad.query_object(data['Name'][ii])
    if result_table == None: # Then the query did not find a match, skip this source
        continue
    rahms = result_table['RA'][0].split()
    decdms = result_table['DEC'][0].split()
    star = SkyCoord('%sh%sm%ss'%(rahms[0],rahms[1],rahms[2]), '%sd%sm%ss'%(decdms[0],decdms[1],decdms[2]),frame='icrs')
    my_table['RA'][ii] = star.icrs.ra.deg
    my_table['Dec'][ii] = star.icrs.dec.deg



# # We are now ready to write our final ASCII (csv) file!

# In[171]:


filename = 'PrimaryDataTable.csv'
fop = open(filename,'w')

# Write the header
headerstr = ''
for cname in colnames:
    headerstr+=cname+','
fop.write(headerstr.strip(',')+'\n')

# Write the data row by row (I'm sure there's a better way..)
starIDcounter = 0 # the ID number you set for the first source from this paper
for ii in range(Nsources):
    writestr = '{0:d},{1:s},{2:s},{3:d},{4:s},{5:.4f},'.format(ii+starIDcounter, my_table['Name'][ii],                                                                my_table['rawstarID'][ii],my_table['RefID'][ii],                                                              my_table['Instrument'][ii],my_table['JD'][ii])
    writestr+= '{0:d},{1:s},{2:s},{3:.6f},{4:.6f},'.format(my_table['FilterID'][ii], my_table['rawRA'][ii],                                                         my_table['rawDec'][ii], my_table['rawQ'][ii], my_table['rawE_Q'][ii])
    writestr+= '{0:.6f},{1:.6f},{2:.6f},{3:.6f},{4:.6f},'.format(my_table['rawU'][ii], my_table['rawE_U'][ii],                                                         my_table['rawP'][ii], my_table['rawE_P'][ii], my_table['rawP_D'][ii])
    writestr+= '{0:.1f},{1:.1f},{2:d},{3:.6f},{4:.6f},'.format(my_table['rawEVPA'][ii], my_table['rawE_EVPA'][ii],                                                         my_table['IntrPol'][ii], my_table['RA'][ii], my_table['Dec'][ii])
    writestr+= '{0:.6f},{1:.6f},{2:.6f},{3:.6f},{4:.6f},'.format(my_table['q'][ii], my_table['e_q'][ii],                                                         my_table['u'][ii], my_table['e_u'][ii], my_table['p'][ii])
    writestr+= '{0:.6f},{1:.1f},{2:.1f}'.format(my_table['e_p'][ii], my_table['evpa'][ii],                                                         my_table['e_evpa'][ii])
    fop.write(writestr+'\n')
    
fop.close()


# # Let's visualize it

# In[172]:


PTable = ascii.read("PrimaryDataTable.csv")


# In[173]:


PTable


# In[ ]:




