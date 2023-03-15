"""
This is the configruation file for setting up iPOLARRIS. Herein, the specifics of the dataset need to be defined such as the experiment, location and type of reading, etc.
Written by Brenda Dolan (CSU) and Anthony Di Stefano (UBC)
Released: May 2017
Last Modified: June 2021
bdolan@atmos.colostate.edu
"""
from __future__ import print_function

import datetime
import glob
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import os
import pandas as pd
import re
import sys
import xarray as xr
import time
from copy import deepcopy
import RadarData
#import GeneralFunctions as GF
import RadarConfig
import plot_driver
from skewPy import SkewT
import re
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

def fix_my_data(ds):
    return(ds.drop(['VTZCS','CVECS']))

def find_dd_match(rdum,ddum,rdate,ddates):

    radlist=[]

    mdfiles = {}
    for v,cname in enumerate(rdum):
        #print cname
        base = os.path.basename(cname)
        dates = rdate[v]
        #print dates, etime,stime
            #print cname
        mval = match_dd(dates,ddates)
        #print( dates,ddates)

        if mval != -999:
            dfile = ddum[mval[0]]
            print('Found DD match!', dfile)
            mdfiles[cname] = dfile
        else:
            mdfiles[cname] = None
        
    return mdfiles

def match_dd(rdate,ddates):
    dum=abs(rdate-np.array(ddates))

    try:
        mval=np.argwhere(dum == np.min(dum))[0]
        diff=np.min(dum)
        if diff.total_seconds() < 600.:
            return mval
        else:
            return -999
    except ValueError: 
        print ('Something DD is not working here!')


def match_snd(rdate,sdates):
    dum=abs(rdate-np.array(sdates))

    try:
        mval=np.argwhere(dum == np.min(dum))
        diff=np.min(dum)
        #allow 12 hours between the radar obs and the sounding
        if diff.total_seconds() < 43200:
            return mval
        else:
            return 'no'
    except ValueError: 
        print ('Something sound is not working here!')


def find_snd_match(config):
    
    sdatetime = ''.join(re.findall(r'\d+', config['sdatetime']))
    sdatetime_format = '%Y%m%d%H%M'
    edatetime = ''.join(re.findall(r'\d+', config['edatetime']))
    edatetime_format = '%Y%m%d%H%M'
    sdt = datetime.datetime.strptime(sdatetime,sdatetime_format)
    edt = datetime.datetime.strptime(edatetime,edatetime_format)

    rdum =[]
    with open(config['rfiles']) as f: 
        for line in f:
            dat = (line)
            rdum.append(foo(dat))
    
    with open(config['sfiles']) as f:
        slist = f.read().splitlines()
    
    sdates=[]
    for v,sname in enumerate(slist):
        base = os.path.basename(sname)
        try: 
            snddate = re.search(r'\d{4}\d{2}\d{2}_\d{2}\d{2}\d{2}',base).group()
            sformat = '%Y%m%d_%H%M%S'
        except AttributeError: 
            snddate = re.search(r'\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}',base).group()
            sformat = '%Y-%m-%d_%H:%M:%S'
        except AttributeError: 
            snddate = re.search(r'\d{4}\d{2}\d{2}\d{2}',base).group()
            sformat = '%Y%m%d%H'
        dates = datetime.datetime.strptime('{r}'.format(r=snddate),sformat)
        sdates.append(dates)

    msfiles = {}

    for v,cname in enumerate(rdum):
        base = os.path.basename(cname)
        try: 
            radcdate = re.search(r'\d{4}\d{2}\d{2}_\d{2}\d{2}\d{2}',base).group()
            rformat = '%Y%m%d_%H%M%S'
        except AttributeError: 
            radcdate = re.search(r'\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}',base).group()
            rformat = '%Y-%m-%d_%H:%M:%S'
        dates = datetime.datetime.strptime('{r}'.format(r=radcdate),rformat)
        if (dates <= edt) and (dates >= sdt):
            mv = match_snd(dates,sdates)
            if mv != 'no':
                msfiles[cname] = np.array(slist)[mv[0][0]]
            else:
                return None
        
    return msfiles

def find_wrfpol_match(config):

    sdatetime = ''.join(re.findall(r'\d+', config['sdatetime']))
    sdatetime_format = '%Y%m%d%H%M'
    edatetime = ''.join(re.findall(r'\d+', config['edatetime']))
    edatetime_format = '%Y%m%d%H%M'
    sdt = datetime.datetime.strptime(sdatetime,sdatetime_format)
    edt = datetime.datetime.strptime(edatetime,edatetime_format)

    rdum =[]
    with open(config['rfiles']) as f: 
        for line in f:
            dat = (line)
            rdum.append(foo(dat))
    
    with open(config['wfiles']) as f:
        slist = f.read().splitlines()
    
    wdates=[]
    for v,sname in enumerate(slist):

        base = os.path.basename(sname)
        wcdate = re.search(r'\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}',base).group()
        wformat = '%Y-%m-%d_%H:%M:%S'
        dates=datetime.datetime.strptime('{r}'.format(r=wcdate),wformat)
        wdates.append(dates)

    msfiles = {}

    for v,cname in enumerate(rdum):
        base = os.path.basename(cname)
        try: 
            radcdate = re.search(r'\d{4}\d{2}\d{2}_\d{2}\d{2}\d{2}',base).group()
            rformat = '%Y%m%d_%H%M%S'
        except AttributeError: 
            radcdate = re.search(r'\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}',base).group()
            rformat = '%Y-%m-%d_%H:%M:%S'
        dates = datetime.datetime.strptime('{r}'.format(r=radcdate),rformat)
        if (dates <= edt) and (dates >= sdt):
            mv = match_snd(dates,wdates)
            if mv != 'no':
                msfiles[cname] = np.array(slist)[mv[0][0]]
            else:
                return None
    
    return msfiles

def foo(s1):
    return '{}'.format(s1.rstrip())

def reduce_dim(ds):
    try:
        t1= ds['time'][0].values
    except KeyError as ke:
        #print(f"{ke} skipping preprocessing")
        return(ds)
    for v in ds.data_vars.keys():
        try:
            ds[v]=ds[v].sel(time=t1).drop('time')
        except KeyError as k:
        #except ValueError as e:
            pass
#            print(e)
#            print(v)
    return(ds)
    
from matplotlib.dates import DateFormatter,HourLocator
dayFormatter = DateFormatter('%H%M')      # e.g., 12
hourFormatter = DateFormatter('%H')      # e.g., 12


def hasNumbers(inputString):
     return any(char.isdigit() for char in inputString)


def polarris_driver(configfile):
    # =====
    # (1) Read in config file line by line.
    # =====

    config = {} # Load variable for config file data

    print('\nReading '+str(configfile)+'...')
    with open(configfile) as f:
        lines1 = [mm for mm in (line.replace('\t',' ') for line in f) if mm]
        lines2 = [nn for nn in (line.strip() for line in lines1) if nn] # NEW! Allow new lines in config file - can be skipped over!
        for line in lines2: #f:
            if not line.startswith("#"):
                key, val, comment = line.split('==')
                vval = val.replace(" ","")
                numck = hasNumbers(vval)
                #if key.replace(" ", "") == 'exper' or key.replace(" ", "") == 'dz_name' or key.replace(" ", "") == 'drop_vars' or key.replace(" ", "") == 'dr_name' or key.replace(" ", "") == 'kd_name' or key.replace(" ", "") == 'rh_name' or key.replace(" ", "") == 'vr_name' or key.replace(" ", "") == 'data' or key.replace(" ", "") == 'xname' or key.replace(" ", "") == 'yname' or key.replace(" ", "") == 'zname' or key.replace(" ", "") == 'latname' or key.replace(" ", "") == 'lonname':
                if key.replace(" ", "") == 'exper' or key.replace(" ", "") == 'drop_vars' or key.replace(" ", "") == 'data':
                    numck = False
                if key.replace(" ", "") == 'exper': # or key.replace(" ", "") == 'ptype':
                    vval = vval.strip("''")
                if key.replace(" ", "") == 'image_dir' or key.replace(" ", "") == 'rfiles':
                    numck = True

                if numck is True or vval == 'None' or vval == 'True' or vval == 'False':
                    try:
                        config[(key.replace(" ", ""))] = eval(vval)
                    except:
                        if "datetime" in vval:
                            config[(key.replace(" ", ""))] = vval
                else:
                    config[(key.replace(" ", ""))] = vval

                if key.replace(" ", "") == 'mix_vars' or key.replace(" ", "") == 'cfad_vars' or key.replace(" ", "") == 'rhi_vars' or key.replace(" ", "") == 'cappi_vars' or key.replace(" ", "") == 'cfad_compare_vars':
                    config[(key.replace(" ", ""))] = val

    print('Read-in complete.\n')

    # =====
    # (2) Find input radar files and concatenate the data. Rename x, y, z variables. 
    # =====

    print('Station/experiment: '+config['exper'])
    if config['type'].startswith('wrf'): print('MP Scheme: '+config['data'].upper())
    else: print('Data Source: '+config['data'].upper())
    print('Start: '+config['sdatetime'])
    print('End: '+config['edatetime'])
    print('')
    time.sleep(3)

    drop_vars = config['drop_vars']
    sdatetime = ''.join(re.findall(r'\d+', config['sdatetime']))
    edatetime = ''.join(re.findall(r'\d+', config['edatetime']))
    rfiles = []
    with open(config['rfiles'], 'r') as f:
        allrfiles = f.read().splitlines()
        for rfile in allrfiles:
            fullname = os.path.basename(rfile)
            print(fullname)
            try: filedatestr = re.search(r'\d{4}\d{2}\d{2}_\d{2}\d{2}\d{2}',fullname).group()
            except AttributeError: filedatestr = re.search(r'\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}',fullname).group()
            filedatestr = filedatestr.replace('_','').replace(':','').replace('-','')
            filedate = filedatestr[0:-2]
            if int(filedate) >= int(sdatetime) and int(filedate) <= int(edatetime):
                rfiles.append(rfile)
    
    print('')
    if rfiles == []:
        print("\nOops! There is no radar data for the dates given in your config file. Exiting...\n")
        sys.exit(1)

    if config['exper'] == 'MC3E' and config['data'] == 'obs':
        print("special handling for ",config['exper'])

        file = open(config['rfiles'], "r")
        rf1=[]
        rf2=[]
        for line in file:
            if re.search('vtzms', line):
                rf1.append(line.rstrip('\n'))
            else:
                rf2.append(line.rstrip('\n'))
        if not rf2:
            rvar = xr.open_mfdataset(rf1,autoclose=True,combine='nested',concat_dim='d',preprocess=fix_my_data)
        else:
            rvar1 = xr.open_mfdataset(rf1,autoclose=True,combine='nested',compat='override',preprocess=fix_my_data)
            rvar2 = xr.open_mfdataset(rf2,autoclose=True,concat_dim='d')
            rvar = xr.concat((rvar1,rvar2),dim='d')
            rfiles = list(np.append(rf1,rf2))         
    else:
        rvar = xr.open_mfdataset(rfiles,autoclose=True,combine='nested',concat_dim='d',preprocess=reduce_dim)
        #rvar = xr.open_mfdataset(rfiles,autoclose=True,concat_dim='d',preprocess=reduce_dim,combine='by_coords')
        #rvar = xr.open_mfdataset(rfiles,autoclose=True,concat_dim='d',preprocess=reduce_dim)
    
    if config['data'].startswith('nexrad'):
      
        rvar = rvar.rename({'x0':'x'})
        rvar = rvar.rename({'y0':'y'}) 
        rvar = rvar.rename({'z0':'z'})
        
        currlat = deepcopy(np.squeeze(rvar['lat0'].values))
        newlat = xr.DataArray(currlat, dims=['y','x'], name='XLAT')
        rvar['XLAT'] = newlat

        currlon = deepcopy(np.squeeze(rvar['lon0'].values))
        newlon = xr.DataArray(currlon, dims=['y','x'], name='XLONG')
        rvar['XLONG'] = newlon
       
        rvar = rvar.reset_coords(['lat0','lon0'])
        rvar = rvar.set_coords(['XLAT','XLONG'])

    elif config['type'].startswith('wrf'):
       
        currx = deepcopy(rvar['x'].values)
        newx = xr.DataArray(currx-np.mean(currx), dims=['x'], name='x')
        rvar['x'] = newx
       
        curry = deepcopy(rvar['y'].values)
        newy = xr.DataArray(curry-np.mean(curry), dims=['y'], name='y')
        rvar['y'] = newy

        if 'd' in rvar['hgt'].dims: currz = rvar['hgt'].values[0,:]
        else: currz = rvar['hgt'].values
        newz = xr.DataArray(currz, dims=['z'], name='z')
        rvar['z'] = newz

        currlat = deepcopy(np.squeeze(rvar['latitude'].values[0,:,:]))
        newlat = xr.DataArray(currlat, dims=['y','x'], name='XLAT')
        rvar['XLAT'] = newlat

        currlon = deepcopy(np.squeeze(rvar['longitude'].values[0,:,:]))
        newlon = xr.DataArray(currlon, dims=['y','x'], name='XLONG')
        rvar['XLONG'] = newlon
        
        rvar = rvar.set_coords(['XLAT','XLONG'])
        
        refvals = deepcopy(rvar['zhh01'].values)
        refvals = np.where(refvals < float(config['refthresh']), np.nan, refvals)
        
        elevs = deepcopy(rvar['elev01'].values) 
        refvals = np.where(np.logical_or(elevs < float(config['mincosthresh']),elevs > float(config['maxcosthresh'])), np.nan, refvals)
        
        newref = xr.DataArray(refvals, dims=['d','z','y','x'], name='zhh01')
        rvar['zhh01'] = newref
        
        parsers = ['zdr01','kdp01','rhohv01','vrad01']
        for v in parsers:
            newvals = deepcopy(rvar[v].values)
            newvals = np.where(newvals == -999.0,np.nan,newvals)
            newvals = np.where(newvals == 0.0,np.nan,newvals)
            newvals = np.where(np.isnan(refvals),np.nan,newvals)
            newvar = xr.DataArray(newvals, dims=['d','z','y','x'], name=v)
            rvar[v] = newvar
 
        qvars = ['qc', 'qr', 'qi', 'qs', 'qg', 'qh']
        for v in qvars:
            qvals = deepcopy(rvar[v].values)
            qvals = np.where(np.logical_or(elevs < float(config['mincosthresh']),elevs > float(config['maxcosthresh'])), np.nan, qvals)
            newq = xr.DataArray(qvals, dims=['d','z','y','x'], name=v)
            rvar[v] = newq
  
    
    # =====
    # (3) Get datetime objects from radar file names.
    # =====

    #if not config['type'].startswith('wrf'): print(rvar['REF'].coords['lat0'])
    #else: print(rvar['zhh01'].coords['lat'])
    #input()

    tm = []
    for d in rfiles:
        base = os.path.basename(d)
        try: 
            radcdate = re.search(r'\d{4}\d{2}\d{2}_\d{2}\d{2}\d{2}',base).group()
            dformat = '%Y%m%d_%H%M%S'
        except AttributeError: 
            radcdate = re.search(r'\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}',base).group()
            dformat = '%Y-%m-%d_%H:%M:%S'
        date=datetime.datetime.strptime(radcdate,dformat)
        tm.append(date)
   
    if drop_vars:
        print("dropping extra variables for memory!")
        rvar= rvar.drop(['vrad03','vdop02','elev03','elev02','vdop03','vang02','vang03','vrad02','zhh02','zhh03','zdr02','zdr03','kdp02','kdp03','rhohv02','rhohv03'])
 
    print('Radar files ready.')
    time.sleep(3)

    # =====
    # (4) 
    # =====

    if config['dd_on']:
        
        print('\nIn your config file, dd_on is set to True.')
        time.sleep(3)
        with open(config['dfiles'], 'r') as g:
            dfiles1 = g.read().splitlines()
        tmd = []
        for d in dfiles1:
            base = os.path.basename(d)
            try: 
                radcdate = re.search(r'\d{4}\d{2}\d{2}_\d{2}\d{2}\d{2}',base).group()
                dformat = '%Y%m%d_%H%M%S'
            except AttributeError: 
                radcdate = re.search(r'\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}',base).group()
                dformat = '%Y-%m-%d_%H:%M:%S'
            except AttributeError:
                radcdate = re.search(r'\d{4}\d{2}\d{2}_\d{2}\d{2}',base).group()
                dformat = '%Y%m%d_%H%M'
            dat2=datetime.datetime.strptime(radcdate,dformat)
            tmd.append(dat2)

        print('Matching Dual-Doppler')
        dmatch = find_dd_match(rfiles,dfiles1,tm,tmd)
        try:
            dvar = xr.open_mfdataset(dfiles1,concat_dim='d')
        except ValueError as ve:
            print('Trying nested instead of concat_dim to read DD files')
            #dvar = xr.open_mfdataset(dfiles1,combine='nested',concat_dim='d')
            dvar = xr.open_mfdataset(dfiles1,autoclose=True,combine='nested',concat_dim='d',preprocess=reduce_dim)

        if config['type'].startswith('obs'):
 
            if 'eastward_wind' in dvar.keys(): dvar = dvar.rename({'eastward_wind':'u'})
            if 'northward_wind' in dvar.keys(): dvar = dvar.rename({'northward_wind':'v'})
            if 'upward_air_velocity' in dvar.keys(): dvar = dvar.rename({'upward_air_velocity':'w'})
 
            if np.array_equal(dvar.variables['x'].values, 1000.0*rvar.variables['x'].values): dvar['x'] = rvar['x']
            if np.array_equal(dvar.variables['y'].values, 1000.0*rvar.variables['y'].values): dvar['y'] = rvar['y']
            if np.array_equal(dvar.variables['z'].values, 1000.0*rvar.variables['z'].values): dvar['z'] = rvar['z']

        elif config['type'].startswith('wrf'):
            currx = deepcopy(dvar['x'].values)
            newx = xr.DataArray(currx-np.mean(currx), dims=['x'], name='x')
            dvar['x'] = newx
            
            curry = deepcopy(dvar['y'].values)
            newy = xr.DataArray(curry-np.mean(curry), dims=['y'], name='y')
            dvar['y'] = newy

            if 'd' in rvar['hgt'].dims: currz = dvar['hgt'].values[0,:]
            else: currz = dvar['hgt'].values
            newz = xr.DataArray(currz, dims=['z'], name='z')
            dvar['z'] = newz
        
        unew = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        unew.fill(np.nan)

        vnew = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        vnew.fill(np.nan)

        wnew = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        wnew.fill(np.nan)

        conv = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        conv.fill(np.nan)
        
        xsubmin = np.where(rvar.variables['x']==np.min(dvar.variables['x']))[0][0]
        xsubmax = np.where(rvar.variables['x']==np.max(dvar.variables['x']))[0][0]

        ysubmin = np.where(rvar.variables['y']==np.min(dvar.variables['y']))[0][0]
        ysubmax = np.where(rvar.variables['y']==np.max(dvar.variables['y']))[0][0]

        zsubmin = np.where(rvar.variables['z']==np.min(dvar.variables['z']))[0][0]
        zsubmax = np.where(rvar.variables['z']==np.max(dvar.variables['z']))[0][0]        

        for q,d in enumerate(dmatch.keys()):
            if dmatch[d] is not None:
                dfile = dmatch[d]
                if dfile in dfiles1:
                    i = dfiles1.index(dfile)
                    unew[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar['u'][i,:,:,:]
                    vnew[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar['v'][i,:,:,:]
                    wnew[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar['w'][i,:,:,:]
      
        if config['type'].startswith('wrf'):
            unew = np.where(np.logical_or(elevs < float(config['mincosthresh']),elevs > float(config['maxcosthresh'])), np.nan, unew)
            #unew = np.where(unew == -999.0, np.nan, unew)
            vnew = np.where(np.logical_or(elevs < float(config['mincosthresh']),elevs > float(config['maxcosthresh'])), np.nan, vnew)
            #vnew = np.where(vnew == -999.0, np.nan, vnew)
            wnew = np.where(np.logical_or(elevs < float(config['mincosthresh']),elevs > float(config['maxcosthresh'])), np.nan, wnew)
            #wnew = np.where(wnew == -999.0, np.nan, wnew)
    
        print(np.nanmin(wnew))
        print(np.nanmax(wnew))
        input()
        rvar['u'] = (['d','z','y','x'],unew)
        rvar['v'] = (['d','z','y','x'],vnew)
        rvar['w'] = (['d','z','y','x'],wnew)
        rvar['conv'] = (['d','z','y','x'],conv)

    print('\nSending data to RadarData...')
 
    if config['wrft_on']:

        print('In your config file, wrft_on is set to True.')
        time.sleep(3)
        if not 't_air' in list(rvar.keys()):
            wmatch = find_wrfpol_match(config)
            if len(wmatch) > 0:
                print('Found POLARRIS-f files!')
                try:
                    tvar = xr.open_mfdataset(list(wmatch.values()),concat_dim='d')
                except ValueError as ve:
                    tvar = xr.open_mfdataset(list(wmatch.values()),combine='nested',concat_dim='d')
                rvar['T'] = tvar['t_air']-273.15
        else:
            rvar['t_air'].values = deepcopy(rvar['t_air'])-273.15
        
        hnew = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        for ii in range(len(rvar['z'])):
            hnew[:,ii,:,:] = rvar['z'][ii]
        newh = xr.DataArray(hnew, dims=['d','z','y','x'], name='z3d')
        rvar['z3d'] = newh

    elif config['snd_on']:

        print('In your config file, snd_on is set to True.')
        time.sleep(3)
        smatch = find_snd_match(config)
        if len(smatch) > 0:
            print ('Found soundings!')
            hnew = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
            for ii in range(len(rvar['z'])):
                hnew[:,ii,:,:] = rvar['z'][ii]
            newh = xr.DataArray(hnew, dims=['d','z','y','x'], name='z3d')
            rvar['z3d'] = newh

            tnew = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
            tnew.fill(np.nan)
            for jj in range(len(tm)):
                sfile = smatch[rfiles[jj]]
                snd = SkewT.Sounding(sfile)
                #rdata.add_sounding_object(snd) # this will add the sounding object to the radar object and then will take the heights and temps
                #rdata.interp_sounding()
                tnew[jj,:,:,:] = np.interp(hnew[jj,:,:,:],snd.data['hght']/1000.0,snd.data['temp'])
            newt = xr.DataArray(tnew, dims=['d','z','y','x'], name='T')
            rvar['T'] = newt
        else:
            print('No soundings available for this study period. Exiting gracefully.')
            sys.exit(1)

   
    if config['type'].startswith('obs') or config['type'].startswith('wrf'):
        rdata = RadarData.RadarData(rvar,tm,ddata = None,band = config['band'],lat_r=config['lat'],lon_r=config['lon'],lat_0=config['lat'],lon_0=config['lon'],exper=config['exper'],rtype=config['type'],rsrc=config['data'],z_thresh=0,conv_types=config['conv_types'],strat_types=config['strat_types'],color_blind=config['cb_friendly'],dd_on=config['dd_on'],rr_on=config['rr_on'],hid_on=config['hid_on'],hid_cats=config['hid_cols'])
    else:
        rdata = RadarData.RadarData(rvar,tm,ddata = None,dz=config['dz_name'],zdr=config['dr_name'],kdp=config['kd_name'],rho=config['rh_name'],temp=config['t_name'],u=Uname,v=Vname,w=Wname,conv=config['convname'],rr=config['rr_name'],band = config['band'],vr = config['vr_name'],lat_r=config['lat'],lon_r=config['lon'],lat=config['latname'], lon=config['lonname'],lat_0=config['lat'],lon_0=config['lon'],exper=config['exper'],rtype=config['type'],rsrc=config['data'],z_thresh=0,conv_types=config['conv_types'],strat_types=config['strat_types'],color_blind=config['cb_friendly'],hid_cats=config['hid_cols'])

    if config['mask_model']:
        print('masking model data')
        rdata.mask_model()
   
    hid_on,qr_on,rr_on = False if config['hid_on'] == '' else True, False if config['qr_on'] == '' else True, False if config['rr_on'] == '' else True
    rdata.calc_pol_analysis(tm,hid_on,qr_on,rr_on,rr_dir=config['rr_dir'],mode=config['type'],classify=config['hid_cols']) # HCA, RR, QR
    
    if not config['cs_z'] == '':
        rdata.calc_cs_shy(cs_z=config['cs_z'])
        rdata.raintype=rdata.data['CSS'].values
    
    if not config['comb_vicr'] == '' and hid_on:
        whvi = np.where(rdata.hid == 6)
        rdata.hid[whvi] = 3
 
#Do some quick masking of the data####
#     mask = np.zeros([rdata.data.dims['d'],rdata.data.dims['z'],rdata.data.dims['y'],rdata.data.dims['x']])
#     whbad = np.logical_or(np.logical_or(np.logical_or(np.logical_or(rdata.data[rdata.dz_name].values>-20.,rdata.data[rdata.zdr_name].values>-2.),rdata.data[rdata.kdp_name].values<10.),rdata.data[rdata.zdr_name].values<10.),rdata.data[rdata.dz_name].values<70.)
#     whbad2= np.where(~whbad)
#     mask[whbad] = 1
#     if np.nanmin(rdata.data['CSS'].values)<1.:
#         mask[rdata.data['CSS'].values<=0] = 0
#     else:
#         mask[np.isnan(rdata.data['CSS'].values)]=0
#
#     rdata.data['CSS'] = rdata.data['CSS'].where(mask ==1)
#     rdata.data[rdata.dz_name].values[whbad2] = np.nan
#     rdata.data[rdata.zdr_name].values[whbad2] = np.nan
#     rdata.data[rdata.kdp_name].values[whbad2] = np.nan
#     rdata.data[rdata.rho_name].values[whbad2] = np.nan
#     rdata.data[rdata.w_name].values[whbad2] = np.nan

    return rdata, config
