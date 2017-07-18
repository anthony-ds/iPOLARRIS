"""
This is the configruation file for setting up iPOLARRIS. Herein, the speicifcs of the dataset need to be defined such as the experiment, location and type of reading, etc.
Written by Brenda Dolan
May 2017
bdolan@atmos.colostate.edu
"""

import glob
import os
import sys
from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import RadarData
import datetime
import matplotlib.pyplot as plt
import numpy as np
import GeneralFunctions as GF
import RadarConfig
import plot_driver
from skewPy import SkewT



###############
def get_data(config, tm, rfile, dmatch,smatch):

#    print tm
    rdata = RadarData.RadarData(rfile,tm,ddata = dmatch,dz = config['dz_name'],zdr=config['dr_name'],
                                              kdp=config['kd_name'],rho=config['rh_name'],temp=config['t_name'],
                                              u=config['uname'],v=config['vname'],w=config['wname'],x=config['xname'],
                                              rr=config['rr_name'],band = config['band'],
                                              y=config['yname'],z=config['zname'],lat=config['lat'], lon=config['lon'],
                                              exper=config['exper'],mphys=config['mphys'],radar_name = config['radarname'],
                                              z_thresh=config['zthresh'],cs_z = config['cs_z'],zconv=config['zconv'],
                                              zdr_offset=config['zdr_offset'], conv_types = config['conv_types'],strat_types = config['strat_types'],
                                              mixed_types = config['mixed_types'])
                                              
    if dmatch is not None:
         wvardum = np.zeros_like(rdata.data[rdata.dz_name].data)
#          wvardum = np.ma.masked_equal(wvardum,0)
#          rdata.data[rdata.w_name].data = wvardum
#          rdata.wvar = wvardum
         rdata.w_name = config['wname']

    if smatch is not None:
        print 'Smatch',smatch
        snd = SkewT.Sounding(smatch)
        rdata.add_sounding_object(snd) # this will add the sounding object to the radar object
                    # and then will take the heights and temps
        rdata.interp_sounding()

    #vvar.close()
    #rdata.radar_name = radarname
    if config['convert_Tk_Tc'] == True:
        rdata.convert_t()
    #print 'Calculating polarimetric fields like HID and rain...'
    if config['pol_on'] == True:
        rdata.calc_pol_analysis()

    return(rdata)
    
def foo(s1):
    return '{}'.format(s1.rstrip())
#######################################################
#######################################################
def get_time(tp,te,filename,dformat):
#    print tp,te
    base = os.path.basename(filename)
    radar_name = base[:3]
    radcdate=np.str(base[tp:tp+te])
#    print radcdate
    date=datetime.datetime.strptime(radcdate,dformat)
    
    return date

def run_exper(config, dmatch = None, smatch=None,interactive=False):

    dum =[]
    with open(config['radar_files']) as f: 
        for line in f:
            dat = (line)
            dum.append(foo(dat))

    dzcfad_a=[]
    drcfad_a=[]
    kdcfad_a=[]
    rhcfad_a=[]
    hidcfad_a=[]

    dzcfadc_a=[]
    drcfadc_a=[]
    kdcfadc_a=[]
    rhcfadc_a=[]
    hidcfadc_a=[]

    dzcfads_a=[]
    drcfads_a=[]
    kdcfads_a=[]
    rhcfads_a=[]
    hidcfads_a=[]


    hts_a=[]
    dzbins_a=[]
    drbins_a=[]
    kdbins_a=[]
    rhbins_a=[]

    times_a=[]
    warea=[]
    wareat=[]
    histzzdr_a = []
    edg_a = []

    histkdz_a = []
    edgkdz_a = []

    histzw_a = []
    edgzw_a = []

    histwr_a = []
    edgwr_a = []


    water_vert = []
    graup_vert =[]
    hail_vert = []
    snow_vert = []
    hid_hts = []
    
    wcfad_a =[]
    wcfadc_a =[]
    wcfads_a =[]
    wbins_a =[]

    hid_cfad = []
    hids_cfad = []
    hidc_cfad = []
    runrr = False
    runw = False

    for d in dum:
        print d
    #        dn = 
        tp = config['time_parse'][0]
        te = config['time_parse'][1]
        tm=get_time(tp,te,d,config['wdate_format'])
#        print tm
        rvar = xr.open_dataset(d,autoclose=True)
#        print dmatch
        if dmatch is not None:
            print dmatch[d]
#            print dmatch.keys()
            if dmatch[d] is not None:
                ddata = xr.open_dataset(dmatch[d],autoclose=True)
            else:
                ddata = None
        else:
            ddata = None
        if smatch is not None:
#            print 'sounding file', smatch[d]
            sfile = smatch[d]
        else:
            sfile = None

        rdat = get_data(config,tm,rvar,ddata,sfile)

        flags = {}
        for k in eval(config['ks']):
            flags[k]=config[k]

        if any(flags.values()) == True:
            plot_driver.make_single_pplots(rdat,flags,config)

        
        bad = np.less(rdat.data[rdat.dz_name].data,config['zthresh'])
        rdat.data[rdat.dz_name].data[bad]=np.nan
    
        whconv= np.where(rdat.raintype != 2)
        whstrat= np.where(rdat.raintype != 1)
        dzcfad, hts, dzbins = GF.cfad(data = rdat.data[rdat.dz_name].data,hts = rdat.data[rdat.z_name][:].data,value_bins = config['dzbins'],
                            ret_z=1,ret_bin = 1,thresh=config['zthresh'])

        dzcfad_a.append(dzcfad)
        dzbins_a.append(dzbins)

        dzcfadc, hts, dzbins = GF.cfad(data = rdat.data[rdat.dz_name].data,mask = whconv,hts = rdat.data[rdat.z_name][:].data,value_bins = config['dzbins'],
                            ret_z=1,ret_bin = 1,thresh=config['zthresh'])
        dzcfadc_a.append(dzcfadc)


        dzcfads, hts, dzbins = GF.cfad(data = rdat.data[rdat.dz_name].data,mask = whstrat,hts = rdat.data[rdat.z_name][:].data,value_bins = config['dzbins'],
                            ret_z=1,ret_bin = 1,thresh=config['zthresh'])
        dzcfads_a.append(dzcfads)


        rdat.data[rdat.zdr_name].data[bad]=np.nan
    
        drcfad, hts, drbins = GF.cfad(data = rdat.data[rdat.zdr_name].data,hts = rdat.data[rdat.z_name][:].data,value_bins = config['drbins'],
                            ret_z=1,ret_bin = 1,thresh=config['zthresh'])
        drcfad_a.append(drcfad)

        drcfadc, hts, drbins = GF.cfad(data = rdat.data[rdat.zdr_name].data,mask = whconv,hts = rdat.data[rdat.z_name][:].data,value_bins = config['drbins'],
                            ret_z=1,ret_bin = 1,thresh=config['zthresh'])
        drcfadc_a.append(drcfadc)

        drcfads, hts, drbins = GF.cfad(data = rdat.data[rdat.zdr_name].data,mask = whstrat,hts = rdat.data[rdat.z_name][:].data,value_bins = config['drbins'],
                            ret_z=1,ret_bin = 1,thresh=config['zthresh'])
        drcfads_a.append(drcfads)
        drbins_a.append(drbins)

        kdcfad, hts, kdbins = GF.cfad(data = rdat.data[rdat.kdp_name].data,hts = rdat.data[rdat.z_name][:].data,value_bins = config['kdbins'],
                            ret_z=1,ret_bin = 1,thresh=config['zthresh'])
        kdcfad_a.append(kdcfad)

        kdcfadc, hts, kdbins = GF.cfad(data = rdat.data[rdat.kdp_name].data,mask = whconv,hts = rdat.data[rdat.z_name][:].data,value_bins = config['kdbins'],
                            ret_z=1,ret_bin = 1,thresh=config['zthresh'])
        kdcfadc_a.append(kdcfadc)

        kdcfads, hts, kdbins = GF.cfad(data = rdat.data[rdat.kdp_name].data,mask = whstrat,hts = rdat.data[rdat.z_name][:].data,value_bins = config['kdbins'],
                            ret_z=1,ret_bin = 1,thresh=config['zthresh'])
        kdcfads_a.append(kdcfads)
        kdbins_a.append(kdbins)
#        if config['wname']  == True:
    
        if config['wname'] in rdat.data.variables.keys():
            print '{c} is good!'.format(c=config['wname'])
            print np.nanmax(rdat.data[rdat.w_name].data)
            wcfad, hts, wbins = GF.cfad(data = rdat.data[rdat.w_name].data,hts = rdat.data[rdat.z_name][:].data,value_bins = config['wbins'],
                                ret_z=1,ret_bin = 1,thresh=-100)
            print np.nanmax(wcfad)
            wcfad_a.append(wcfad)

            wcfadc, hts, wbins = GF.cfad(data = rdat.data[rdat.w_name].data,mask = whconv,hts = rdat.data[rdat.z_name][:].data,value_bins = config['wbins'],
                                ret_z=1,ret_bin = 1,thresh=-100)
            wcfadc_a.append(wcfadc)

            wcfads, hts, wbins = GF.cfad(data = rdat.data[rdat.w_name].data,mask = whstrat,hts = rdat.data[rdat.z_name][:].data,value_bins = config['wbins'],
                                ret_z=1,ret_bin = 1,thresh=-100)
            wcfads_a.append(wcfads)

            wbins_a.append(wbins)

            histzw,edgzw = GF.hist2d(rdat.data[rdat.w_name].data,rdat.data[rdat.dz_name].data,binsx = config['wbins'],binsy=config['dzbins'])
    #        edgx = 
            histzw_a.append(histzw)
            edgzw_a.append(edgzw)
            runw = True
        else:
            wbins = config['wbins']

        if config['rr_name'] in rdat.data.variables.keys():
                histwr,edgwr = GF.hist2d(rdat.data[rdat.rr_name].data,rdat.data[rdat.dz_name].data,binsx = config['rrbins'],binsy=config['dzbins'])
    #        edgx = 
                histwr_a.append(histwr)
                edgwr_a.append(edgwr)
                runrr = True
        else:
            edgwr_a.append(None)

        hts_a.append(hts)
        times_a.append(tm)
    
        hist,edg = GF.hist2d(rdat.data[rdat.zdr_name].data,rdat.data[rdat.dz_name].data,binsx = config['drbins'],binsy=config['dzbins'])
#        edgx = 
        histzzdr_a.append(hist)
        edg_a.append(edg)

        histkdz,edgkdz = GF.hist2d(rdat.data[rdat.kdp_name].data,rdat.data[rdat.dz_name].data,binsx = config['kdbins'],binsy=config['dzbins'])
#        edgx = 
        histkdz_a.append(histkdz)
        edgkdz_a.append(edgkdz)

        rvar.close()
        if dmatch is not None:
            if dmatch[d] is not None:
                ddata.close()
        
        hidhts, mwrf_water_vert = GF.hid_vertical_fraction(rdat.hid,rdat.data[rdat.z_name][:].data,config['hidwater'],rdat.species,z_resolution =config['z_resolution'],z_ind=1)
        hidhts, mwrf_graup_vert = GF.hid_vertical_fraction(rdat.hid,rdat.data[rdat.z_name][:].data,config['hidgraup'],rdat.species,z_resolution =config['z_resolution'],z_ind=1)
        hidhts, mwrf_hail_vert = GF.hid_vertical_fraction(rdat.hid,rdat.data[rdat.z_name][:].data,config['hidhail'],rdat.species,z_resolution =config['z_resolution'],z_ind=1)
        hidhts, mwrf_snow_vert = GF.hid_vertical_fraction(rdat.hid,rdat.data[rdat.z_name][:].data,config['hidsnow'],rdat.species,z_resolution =config['z_resolution'],z_ind=1)

        water_vert.append(mwrf_water_vert)
        graup_vert.append(mwrf_graup_vert)
        hail_vert.append(mwrf_hail_vert)
        snow_vert.append(mwrf_snow_vert)
    
        hidhts,hiddum =GF.hid_cdf(rdat.hid, rdat.data[rdat.z_name][:].data,rdat.species,z_resolution=config['z_resolution'], pick=None,z_ind =0)
        hid_cfad.append(hiddum)

        hidhts,hiddumc =GF.hid_cdf(rdat.hid, rdat.data[rdat.z_name][:].data,rdat.species,z_resolution=config['z_resolution'], pick=None,z_ind =0,mask=whconv)
        hidc_cfad.append(hiddumc)

        hidhts,hiddums =GF.hid_cdf(rdat.hid, rdat.data[rdat.z_name][:].data,rdat.species,z_resolution=config['z_resolution'], pick=None,z_ind =0,mask = whstrat)
        hids_cfad.append(hiddums)

#        print np.shape(hiddum)
        hid_hts.append(hidhts)
        
        
        if config['wname'] in rdat.data.variables.keys():
#            print np.shape(rdat.T)
            tmp, m_warea_wrf = GF.updraft_width_profile(rdat.data[rdat.w_name].data,rdat.data[rdat.z_name].data,thresh=config['wthresh'], temps=config['trange'],\
                z_ind=0,tcoord = True,temp = rdat.T[:,0,0])
            warea_wrf = m_warea_wrf*rdat.dx*rdat.dy/rdat.ntimes
            if rdat.data[rdat.x_name].units == "[deg]":
                ware_wrf=warea_wrf*110.*110.


            warea.append(warea_wrf)
            wareat.append(tmp)
        
        rconf = RadarConfig.RadarConfig()
        rconf.date = times_a
        rconf.mphys = rdat.mphys
        rconf.exper = rdat.exper
        rconf.radar_name = rdat.radar_name
        
        if interactive is True:
            return rdat
            
        del rdat


    dat = {'rconf': rconf,
          'dzcfad':dzcfad_a,
          'drcfad':drcfad_a,
          'kdcfad':kdcfad_a,
          'wcfad':wcfad_a,

          'dzccfad':dzcfadc_a,
          'drccfad':drcfadc_a,
          'kdccfad':kdcfadc_a,
          'wccfad':wcfadc_a,

          'dzscfad':dzcfads_a,
          'drscfad':drcfads_a,
          'kdscfad':kdcfads_a,
          'wscfad':wcfads_a,

          'hts':hts_a,
          'histzzdr':histzzdr_a,
          'edgzzdr':edg_a,
          'histkdz':histkdz_a,
          'edgkdz':edgkdz_a,
          'histzw':histzw_a,
          'edgzw':edgzw_a,
          'histwr':histwr_a,
          'edgwr':edgwr_a,
          'hidhts':hid_hts,
          'water_vert':water_vert,
          'graup_vert':graup_vert,
          'hail_vert':hail_vert,
          'snow_vert':snow_vert,
           'warea':warea,
           'wareat':wareat,
          'time':times_a,
          'dzbins':dzbins,
          'drbins':drbins,
          'kdbins':kdbins,
          'wbins':wbins,
          'dzcbins':dzbins,
          'drcbins':drbins,
          'kdcbins':kdbins,
          'wcbins':wbins,

          'dzsbins':dzbins,
          'drsbins':drbins,
          'kdsbins':kdbins,
          'wsbins':wbins,

          'hidcfad':hid_cfad,
          'hidccfad':hidc_cfad,
          'hidscfad':hids_cfad,
          'runrr':runrr,
          'runw':runw}
    return dat
