# Brody Fuchs, CSU, October 2016
# brfuchs@atmos.colostate.edu

# Classes and function that read in polarimetric/DD files
# and plot them and do things with them

# Also included here is a Cell object to do some cell calculations
# which will be used for my PhD stuff

# 10/13/16: This is hooked into my old version of things like HID, will need to get 
#             csuradartools into this at some point

# 2/22/2016: Trying to modularize to allow for integration with larger datasets such as 
# WRF output and other functions.

from __future__ import division
from __future__ import print_function
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#plt.style.use('./presentation.mplstyle')
import os
import sys
#from pyhid import beta_functions, cdf_fhc, radar_calculations
from matplotlib import colors
import matplotlib.ticker as plticker
#from general_tools import make_cmap
import scipy.ndimage as ND
import matplotlib as mpl
from copy import deepcopy
import datetime
import scipy.ndimage as ndi
import xarray as xr
import raintype as rt
import steiner_houze_yuter_cs as shy

from pyproj import Proj, transform

from tqdm import tqdm

#import analysis_tools as AT
#import lightning_tools as LT
#from CSU_RadarTools.csu_radartools import csu_fhc
import csu_fhc
#import general_tools as gentools
import RadarConfig
# Up here are just some general functions

np.set_printoptions(threshold=sys.maxsize)

class RadarData(RadarConfig.RadarConfig): 

    def __init__(self, data,times, ddata = None,dz='DZ', zdr='DR', kdp='KD', ldr='LH', rho='RH', hid='HID',conv='Con',temp='T', x='x', y='y', z='z', u='U', v='V', w='Wvar', rr='RR',vr='VR',lat=None, lon=None, band='C',exper='CASE',lat_r=None,lon_r=None,mphys=None,dd_data = None,z_thresh=-10.0,cs_z = 2.0,zconv = 41.,zdr_offset=0, remove_diffatt = False,lat_0 = 0.0,lon_0=90.0,conv_types = ['CONVECTIVE'],strat_types = ['STRATIFORM'],mixed_types = ['UNCERTAIN'],mixr=['qr','qs','qc','qi','qh','qg'],return_scores=False,color_blind=False): 

        super(RadarData, self).__init__(dz=dz, zdr=zdr, kdp=kdp, ldr=ldr, rho=rho, hid=hid, conv=conv,temp=temp, x=x, y=y,lat_0=lat_0,lon_0=lon_0,lat_r=lat_r,lon_r=lon_r, z=z, u=u, v=v, w=w,rr=rr,vr=vr,mphys=mphys,exper=exper,lat=lat,lon=lon,tm = times,color_blind=color_blind)

        # ********** initialize the data *********************
#        self.data = {} 
        #self.date = None
#        self.radar_file = radar_file
#        self.dd_file = dd_file

        # now if radar_file and dd_file have been set, go and try to read them in, else do nothing
#        if self.radar_file is not None:
#            self.time_parse = time_parse
#            self.dformat = dformat
        self.conv_types = conv_types
        self.strat_types = strat_types
        self.mixed_types = mixed_types
        self.return_scores=return_scores
        self.band = band
        self.mixr=mixr
        self.lat_0 = lat_0
        self.lon_0 = lon_0
        if ddata is not None:
            self.data= data.combine_first(ddata)
        else:
            self.data = data
        self.zdr_offset = zdr_offset
        if zdr_offset != 0:
            print ('fixing Zdr!')
            self.zdr_correct()
        self.cs_z =cs_z
        self.zconv =zconv
        self.t_name = temp
        self.z_thresh=z_thresh
        self.badval = -999.0
        self.zind = 1
        self.yind = 2
        self.xind = 3
        self.ntimes =1
        if 'd' in self.data[self.z_name].dims:
            self.nhgts = len(self.data[self.z_name].values[0])
        else:
            self.nhgts = len(self.data[self.z_name].values)
            #self.data[self.z_name].values = self.data[self.z_name].values[0]
#        if 'd' in self.data[self.z_name].dims:
#            self.nhgts = np.shape(self.data[self.z_name].values)[self.zind][0]
#        else:
#            self.nhgts = np.shape(self.data[self.z_name].values)
#        print(self.data[self.z_name].values)
#        print(self.data[self.dz_name].sizes[self.z_name])
#        self.nhgts = self.data[self.dz_name].sizes[self.z_name]
#            self.read_data_from_nc(self.radar_file)
        print ('calculating deltas')
        self.calc_deltas()
        print ('calculating rain area')
        self.radar_area()
        self.rr_name = rr
        if self.rr_name == None:
            self.rr_name = 'RR'

       #print 'masking data'
        #print('masking data')
        #self.mask_dat()
        if remove_diffatt == True:
            self.corr_zdr()
#        self.raintype_calc()
#        self.calc_pol_analysis()
        #print 'going to HID'
#        self.set_hid(use_temp = 'True',band=self.band,zthresh = self.z_thresh)

#        else:
#            pass

        # down here goes some zdr checking
        self.zdr_offset = 0 # initialize as 0
        try:
            self.top_index = np.where(self.data[self.z_name]== np.max(self.data[self.z_name]))[self.zind]
        except:
            self.top_index = np.where(self.data[self.z_name]== np.max(self.data[self.z_name]))[0]
        try:
            self.bot_index = np.where(self.data[self.z_name]== np.min(self.data[self.z_name]))[self.zind]
        except:
            self.bot_index = np.where(self.data[self.z_name]== np.min(self.data[self.z_name]))[0]

#         if hasattr(self.nc,'Latitude_deg') == True:
#             print 'Getting attribute!'
#             try:
#                 self.radar_lat = getattr(self.nc,'Latitude_deg')[0]
#             except:
#                 self.radar_lat = getattr(self.nc,'Latitude_deg')
#             
#         if hasattr(self.nc,'Longitude_deg') == True:
#             print 'Getting attribute!'
#             try:
#                 self.radar_lat = getattr(self.nc,'Longitude_deg')[0]
#             except:
#                 self.radar_lat = getattr(self.nc,'Longitude_deg')
        try:
            self.T = self.data[temp]
        except:
            print ('')
#        self.set_masks()
#        self.remove_nc()
    # Have some options of how to pass it files? 
    # Nc gridded files
    # pyart?
    # some other format??

    def read_data_from_nc(self, filename, format = 'grid', date = None,dd_flag =0):
    # will need to change time_parse if filename is different, also assumes 4 digits for year
    # this will read in a netcdf file calling functions below assuming a certain filename format
    # RADyyyymmdd_hhmmss*.ext
        tp = self.time_parse[0]
        te = self.time_parse[1]
        self.filename = filename
        self.base = os.path.basename(filename)
        self.radar_name = self.base[:3]
        if date is None:
            try:
                radcdate=np.str(self.base[tp:tp+te])
                #print radcdate
                self.date=datetime.datetime.strptime(radcdate,self.dformat)
            except ValueError:
                pass

        # run some basic functions so not needed to be typed out
        self.nc4_read()
        self.nc4_varnames()
        
#        self.nc4_dict()

#############################################################################################################

    def nc4_read(self):
        "Read netcdf4 file with netCDF4"
        self.nc = Dataset(self.filename)
    #############################################################################################################

    def nc4_varnames(self):
        "Returns all of the variables of the netcdf4 file"
        self.ncvariables = self.nc.variables.keys()

    #############################################################################################################

    def raintype_calc(self):
        print ('Setting up default C/S parameters. Can change in RadarData.py')
        minZdiff = 20 
        deepcoszero = 40
        shallowconvmin = 28
        #truncZconvthres = 43;
        #dBZformaxconvradius = 46;
        truncZconvthres = self.zconv
        dBZformaxconvradius = self.zconv+3;
        weakechothres = 7
        backgrndradius = 5       #(in km)
        maxConvRadius = 10       #(in km)
        minsize = 8              #(in km^2)
        startslope = 50          #(in km^2)
        maxsize = 2000           #(in km^2)

        if self.x_name == 'longitude':
            dx = self.dx*110.
        else:
            dx = self.dx
#        print dx
        zlev = self.get_ind(self.cs_z,self.data[self.z_name].values)
#        print zlev

        print ('Unfortunatley have to run the raintype per file. Might take a minute....{n}'.format(n=self.data.dims['d']))
        rntypetot = []
        for q in tqdm(range(self.data.dims['d'])):        
        
            refl = np.squeeze(self.data[self.dz_name].sel(d=q,z=slice(zlev,zlev)).values)
            #print 'refl shape',np.shape(refl)
            if len(np.shape(refl)) > 3:
                refl = np.squeeze(self.data[self.dz_name].sel(z=slice(zlev,zlev+1)).values)
            #print 'refl shape',np.shape(refl)
            #print type(refl)
            refl.setflags(write=1)
            refl[(np.isnan(refl))] = -99999.9
            refl_missing_val = -99999.9
            #print self.data[self.dz_name].min()
            raintype,self.rtypes = rt.raintype(refl, refl_missing_val=refl_missing_val,
                                       refl_dx=dx, minZdiff=minZdiff, deepcoszero=deepcoszero,
                                       shallowconvmin=shallowconvmin,truncZconvthres=truncZconvthres,
                                       dBZformaxconvradius=dBZformaxconvradius,
                                       weakechothres=weakechothres, backgrndradius=backgrndradius,
                                       maxConvRadius=maxConvRadius,minsize=minsize,
                                       startslope=startslope, maxsize=maxsize)
            nlevs = np.shape(self.data[self.z_name].data)[0]
    #        print nlevs
            rpt = np.tile(raintype,(nlevs,1,1))
            rntypetot.append(rpt)
        self.raintype= np.array(rntypetot)
        #self.def_convstrat()
        self.add_field((self.data[self.dz_name].dims,self.raintype,), 'RRT')    
        
#############################################################################################################

    def check_size(self,ddata):
        if self.x_name in self.data:
            ochksz = np.shape(self.data[self.x_name])
#            print self.nc.variables[self.x_name]
            ddsz = np.shape(ddata.variables[self.x_name])
            
            if ochksz != ddsz:
                "Create a new array of the dimensions for the pol file to put DD data in"
                xvals = ddata.variables[self.x_name][:]
                yvals = ddata.variables[self.y_name][:]
                zvals = ddata.variables[self.z_name][:]
                
                whx = np.nonzero(np.in1d(self.data[self.x_name], xvals))[0]
                why = np.nonzero(np.in1d(self.data[self.y_name], yvals))[0]
                whz = np.nonzero(np.in1d(self.data[self.z_name], zvals))[0]

                
                for key in ddata.variables.keys():
                    if ddata.variables[key].ndim == 3:
                        test =np.zeros_like(self.data[self.dz_name])
                        test = np.ma.masked_where(test==0,test)
                        test[whz[0]:whz[-1]+1,why[0]:why[-1]+1,whx[0]:whx[-1]+1] = ddata.variables[key][:]
    #                    for i,v in enumerate(whz):
    #                        for j,u in enumerate(why):
    #                            for k,w in enumerate(whx):
    #                                test[v,u,w]=self.nc.variables[key][i,j,k]
#                        print np.shape(test)
                        if not key in self.data.variables.keys():
#                            print key
                            self.data[key]=((self.data.dims['z'],self.data.dims['y'],self.data.dims['x']),test)
                
                
    #############################################################################################################

    def nc4_dict(self):
        "Goes into the netcdf file and grabs all of the data and puts into a dictionary"
        #out = {}
        for key in self.nc.variables.keys():
            self.data[key] = self.nc.variables[key][:]

    #self.data[data_type] = out
    #############################################################################################################

    def get_data_from_dict(self, dict):    
        # if want to pass a dictionary already
        self.data = dict
#############################################################################################################
    def corr_zdr(self):
        #whzdr = np.where(np.logical_and(self.data[self.zdr_name].data<-1., self.data[self.kdp_name].data >.0))
        whzdr = np.where(self.data[self.zdr_name].data<-1.)
        self.data[self.zdr_name].data[whzdr]= np.nan
        

    def mask_dat(self):    
        # if want to pass a dictionary already
        for k in self.data.variables.keys():
#            print k, np.shape(self.data[self.dz_name].data), np.shape(self.data[k].data)
            if k != 'TIME' and k != self.t_name:
                #print(k)
#                print np.shape(self.data[self.dz_name].data), np.shape(self.data[k].data)
                try:
                    whbad = np.where(np.less_equal(self.data[self.dz_name].values, self.z_thresh))
                    self.data[k].data[whbad]= np.nan
                    whbad2 = np.where(np.isnan(self.data[self.dz_name].values))
                    self.data[k].values[whbad2]= np.nan
#                    self.data[k].data = np.ma.masked_where(self.data[self.dz_name].data < self.z_thresh,self.data[k].data)

                    whbad3 = np.where(np.less_equal(self.data[k].values, self.badval))
                    self.data[k].values[whbad3] = np.nan
                    whbad4 = np.where(np.isnan(self.data[k].values))
                    self.data[k].values[whbad4] = np.nan
                except:
                    pass
#                    print np.shape(self.data[self.dz_name].data), np.shape(self.data[k].data), 'wrong shapes!'
#############################################################################################################

    def convert_t(self):    
        # if want to pass a dictionary already
        #self.T = np.ma.masked_less(self.T, 0)
        #self.T.values=np.ma.masked_less(self.T.values,0)
        self.T = self.T #- 273.15
        
        #self.T.values = self.T.values-273.15
        
#############################################################################################################
    def set_masks(self):
        for k in self.data.keys():
            self.data[k] = np.ma.masked_where(self.data[k] < -998.0,self.data[k])

    def radar_area(self):
        #Define the radar coverage area. We can do this from radial velocity in the model
        if self.mphys == 'obs':
            self.radar_area = len(self.data[self.x_name].values)*len(self.data[self.y_name].values)*self.dx*self.dy  
        else:
            vrcomp = self.data[self.vr_name].sel(d=0).max(axis=0).values
            whmask = np.where(vrcomp > -50)
            
            if 'd' in self.data[self.x_name].dims:
                x,y = self.convert_ll_to_xy(self.data[self.lat_name].sel(d=0),self.data[self.lon_name].sel(d=0))
                dx = np.average(np.diff(x.sel(y=0)))
                dy = np.average(np.diff(y.sel(x=0)))
            else:
                x,y = self.convert_ll_to_xy(self.data[self.lat_name],self.data[self.lon_name])
                dx = np.average(np.diff(x))
                dy = np.average(np.diff(y))

            self.x = x.values
            self.y = y.values
            self.dx = dx
            self.dy = dy
            dummy=np.zeros_like(vrcomp)
            dummy[whmask] = 1
            radar_area = np.count_nonzero(dummy)*dy*dx
            self.radar_area = radar_area
            print('radar _area',radar_area)
            
            dummy = 0.0
            vrcomp = 0.0

    def mask_model(self):
     
        whmask2= np.where(np.logical_or(self.data[self.vr_name] < -50.,self.data[self.dz_name]<=0))

        mask_dat=[self.dz_name,self.zdr_name,self.vr_name,self.rr_name,self.kdp_name,self.w_name,self.u_name,self.v_name]
        for k in mask_dat:
            print(k)
            dathold = self.data[k].values
            dathold[whmask2]=np.nan
            
            self.data[k].values=dathold
            dathold = 0.0

    def valid_vars(self):
        return np.intersect1d(self.pol_vars, self.data.keys())
#############################################################################################################


    def calc_deltas(self): # get grid sizes for x, y, z
        
        if 'd' in self.data[self.x_name].dims:
            print,'calc deltas'
            self.dx = np.average(abs(np.diff(self.data[self.x_name].sel(d=0,y=0))))
            self.dy = np.average(abs(np.diff(self.data[self.y_name].sel(d=0,x=0))))        
            self.dz = np.average(abs(np.diff(self.data[self.z_name].sel(d=0))))
        else:
            self.dx = np.average(abs(np.diff(self.data[self.x_name].values)))
            self.dy = np.average(abs(np.diff(self.data[self.y_name].values)))
            self.dz = np.average(abs(np.diff(self.data[self.z_name].values)))
            


#############################################################################################################


    def add_field(self, in_array, name, pol_var=False):
        "Add some field to the radar object and give it a name"
#        self.data[name] = (['d','z','y','x'],in_array)
        self.data[name] = in_array
        if pol_var:
            self.add_pol_var(name)
#############################################################################################################

    def add_pol_var(self, name):
        self.pol_vars = np.append(self.pol_vars, name)
#############################################################################################################

    def get_names(self):
        pass

#############################################################################################################

    ##### some private worker bee functions #####

    def _get_ab_incides(self, above=None, below=None):

        if above is not None:
            bot_index = np.argsort(abs(self.data[self.z_name].values - above))[0]
        else:
            bot_index = deepcopy(self.bot_index)
        if below is not None:
            top_index = np.argsort(abs(self.data[self.z_name].values - below))[0]
        else:
            top_index = deepcopy(self.top_index)

        return int(bot_index), int(top_index)


    def _pick_data(self, data, pick):

        if pick is None:
            return data
        else:
            # want to mask the values that are false, hence the logical not
            return np.ma.masked_where(np.logical_not(pick), data) 


    def convert_xy_to_ll(self, x_data, y_data):
    # code to use the radar lat/lon to convert x and y coords to lat lon
        lat = self.radar_lat + y_data/111.0
        lon = self.radar_lon + x_data/(111.0*np.cos(self.radar_lat*np.pi/180.0))

        pass


    def convert_ll_to_xy(self, lat_data, lon_data):
        # code to use the radar lat/lon to convert lat/lon data to x/y radar relative coords
        # just going to assume a flat plane for now
        if (self.radar_lat is not None) and (self.radar_lon is not None):
            dlat = lat_data - self.radar_lat
            dlon = lon_data - self.radar_lon    

            x = dlon*111.0*np.cos(self.radar_lat*np.pi/180.0)
            y = dlat*111.0

            return x, y

        else:
            print ('Need to set the radar lat/lon')
            return None, None
    


#############################################################################################################

    ################### Some quality control stuff ###############################################

    def zdr_check(self, bins=np.arange(-3, 3, 0.15), thresh=0.3):
            # This will do a binned histogram of the ZDR everywhere and makes sure that it peaks around 0
        hist, edges = np.histogram(np.isfinite(self.data[self.zdr_name].data), bins=bins)
        dbin = np.diff(bins)[0]/2.0
        maxarg = np.argmax(hist)

        offset = edges[maxarg]+dbin
#        print hist, edges
        self.zdr_offset = deepcopy(offset)
        if abs(offset) >= thresh:
            self.zdr_correct()

    def zdr_correct(self):
        self.data[self.zdr_name].data -= self.zdr_offset




    
    def find_melting_level(self):
        # this sets the levels that are above and below the melting level
        mlevel = self.get_T_height(0)

        # now need to 
        self.above_melting = np.zeros((self.data[self.dz_name].shape), dtype=bool)
        self.below_melting = np.zeros((self.data[self.dz_name].shape), dtype=bool)

        # alright, now need to loop thru the vertical
        for iz in range(len(self.data[self.z_name])):
            if self.data[self.z_name][iz] < mlevel:
                self.below_melting[iz, :, :] = True
            elif self.data[self.z_name][iz] > mlevel:
                self.above_melting[iz, :, :] = True
            else:
                pass # this only happens if they're equal I think


#############################################################################################################

    def get_max_loc(self, variable, smooth=0, above=None, below=None):
        # this will find the maximum updraft and return the x, y, z position
        
        #print smooth
        bot_index, top_index = self._get_ab_incides(above=above, below=below)

        if smooth > 0:
            arr = ndi.filters.gaussian_filter(self.data[variable], smooth)[bot_index:top_index]
        else:
            arr = self.data[variable][bot_index:top_index]

        max_ind = np.where(arr == arr.max())
        zloc = self.data[self.z_name][bot_index:top_index][max_ind[0]][0]
        yloc = self.data[self.y_name][max_ind[1]][0]
        xloc = self.data[self.x_name][max_ind[2]][0]

        return zloc, yloc, xloc

    def remove_nc(self):
        del self.nc


    ##################### Here is some sounding stuff ####################################

    def add_sounding_object(self, sobj):
        self.sobj = sobj
        # add the object to be used later or for other purposes
        # now call add_sounding from the sounding object info
        self.add_sounding(height=sobj.data['hght']/1000.0, temp=sobj.data['temp']) # have to flip it to make right side up




    def add_sounding(self, height=None, temp=None, interp=False):
        """This will add sounding data and interp to the radar heights, if you leave data=None, it will put in a basic
        default sounding"""
        if (height is not None) and (temp is not None):
            self.snd_temp = deepcopy(temp)
            self.snd_height = deepcopy(height)
        else: # do a default temp structure
            t0 = 25.0 # base temp
            z0 = 0.5 # base height
            gamma = 6.5 # lapse rate
            nz = 20

            self.snd_height = np.arange(20)+z0
            self.snd_temp = -1*gamma*np.arange(20)+t0
            self.snd_temp[self.snd_temp <= -60] = -60 # make it -60 as the coldest, don't think this will matter too much

#############################################################################################################


    def interp_sounding(self):

        """This will take the sounding data and interpolate it to the radar coordinates"""
        self.gridded_height = np.zeros(self.data[self.dz_name].shape)
        #print(self.x)
        if 'd' in self.data[self.z_name].dims:
            for i in range(self.data[self.z_name].shape[1]):
                self.gridded_height[:,i,:,:] = self.data[self.z_name][0,i]
        else:
            for i in range(self.data[self.z_name].shape[0]):
                #self.gridded_height[:,:,i,...] = self.data[self.z_name][i]
                self.gridded_height[:,i,:,:] = self.data[self.z_name][i]

        self.T = np.interp(self.gridded_height, self.snd_height, self.snd_temp)
        self.add_field((self.data[self.dz_name].dims,self.T,), self.temp_name)

#############################################################################################################


    def get_T_height(self, temp, interp=False):
        temp_index = np.argmin(abs(self.T[:,0,0] - temp))
        return self.gridded_height[0,:,0,0][temp_index]


#############################################################################################################
   ############ Here is calling CSU radartools for HID, RR, etc... ############################
#############################################################################################################
    def calc_pol_analysis(self,tm,hid_on,qr_on,rr_on,rr_dir=None,**kwargs):
        
        import time
        
        if hid_on:

            start = time.time()
            self.set_hid(use_temp = 'True',band=self.band,zthresh = self.z_thresh,return_scores=self.return_scores)
            print('HID runtime = '+str(time.time()-start))
        
        if self.mphys == 'obs':
            
            if qr_on:

                print("running pol qr")
                start = time.time()
                self.calc_qr_pol()
                print('QR runtime = '+str(time.time()-start))

            if rr_on:
                
                print("running pol rr")
                start = time.time()
                self.calc_rr_pol(tm,rr_dir,**kwargs)
                print('RR runtime = '+str(time.time()-start))


#############################################################################################################
   ############ Here is calling CSU radartools for HID, RR, etc... ############################
#############################################################################################################


   # Just a wrapper on the CSU radartools HID function
    def set_hid(self, band=None, use_temp=False, name='HID',zthresh = -9999.0,return_scores=False):
        
        #        print zthresh
        #        print self.dz_name
        #       bad = self.data[self.dz_name].where(self.data[self.dz_name].values<zthresh)
        #       print type(self.data[self.dz_name])
        if band is None:
           self.hid_band = self.band
        else:
           self.hid_band = band
        
        hid = []
        scores=[]
        print ("Unfortunately need to run HID by time")
        for v in tqdm(range(len(self.data[self.dz_name]))):
            dzhold =np.squeeze(self.data[self.dz_name].sel(d=v)).values
#             drhold =np.squeeze(self.data[self.zdr_name].sel(d=v)).values
#             kdhold = np.squeeze(self.data[self.kdp_name].sel(d=v)).values
#             rhhold = np.squeeze(self.data[self.rho_name].sel(d=v)).values
#            print('shape holds',np.shape(dzhold))
            if use_temp and hasattr(self, 'T'):
               #print ('Using T!')
               #tdum = self.T[v,...]
               tdum = self.T[v,:,:,:]
               
               #print(type(tdum),'tdum is')
               #print('T:',np.shape(tdum))
            else:
               tdum = None
            #print('You have entered hid band:',self.hid_band, 'ln 624 RadarData')
            hiddum = csu_fhc.csu_fhc_summer(dz=dzhold, zdr=np.squeeze(self.data[self.zdr_name].sel(d=v)).values, rho=np.squeeze(self.data[self.rho_name].sel(d=v)).values, kdp=np.squeeze(self.data[self.kdp_name].sel(d=v)).values, band=self.hid_band, use_temp=True, T=tdum, return_scores=self.return_scores)
#            scores.append(scoresdum)
            #hiddum = np.argmax(scoresdum,axis=0)+1
#            print(np.shape(tdum),'tdum')
            #whbad = np.where(np.logical_and(hiddum ==1,tdum <-5.0))
            if tdum.any() == None: whbad = np.where(np.logical_and(hiddum == 1,tdum == None))
            else: whbad = np.where(np.logical_and(hiddum == 1,tdum < -5.0))
            hiddum[whbad] = -1
            dzmask = np.where(np.isnan(dzhold))
            hiddum = np.array(hiddum,dtype='float64')
            hiddum[dzmask] =np.nan
            hid.append(hiddum)
#        print "Returned to RadarData"
#        self.scores=np.array(scores)
        #print 'np.shape self.scores',np.shape(self.scores)
#         #       self.data[self.dz_name].values[bad] = np.nan
#        dzmask = np.where(np.isnan(self.data[self.dz_name].values))
#            # set the hid
#        self.hid = np.argmax(scores, axis=1)+1
#        print('dzmask shape',np.shape(dzmask))
 #       hid[dzmask] = np.nan
        self.hid = np.array(hid)
#        dzmask=0.0
#        print ('setting bad hid with Drizzle, T<-5')

#       = try:
#         #           print 'Trying to mask HID!'
#        self.hid = np.ma.masked_where(dzmask==True,self.hid)
        
#         self.hid = np.ma.masked_where(self.T.mask,self.hid)
#         self.hid = np.ma.masked_where(self.data[self.dz_name].mask,self.hid)
#        except:
#            print 'Problem trying to threshold HID'
#        try:
#         scmask = np.isnan(self.scores[:,0,...])
# #         #            print 'trying to mask via scores!'
#         self.hid[scmask] =-1
#        except:
#             print 'Cant threshold on scores!'
#         #           self.hid=np.ma.masked_where(self.data[self.dz_name].data < self.z_thresh,self.hid)
        self.add_field((self.data[self.dz_name].dims,self.hid,), name)

#############################################################################################################
        ###Calculate the mixing ratios and add to the radar object

    def calc_qr_pol(self):

        import time
    ### This is where I do mixing ratio calculations ### 
        #print self.band
        #print self.expr
#        print self.band, self.exper
        if self.band == 'C' or self.band == 'c':
 #           print 'got C!'
#            print (self.exper)
            if self.exper == 'MC3E':
                k_c= 1.4195
                k_m=0.7489

                z_c=0.0014297
                z_m=0.6729
            elif self.exper == 'TWPICE': 
                k_c= 0.89191
                k_m=0.6529

                z_c=0.00194
                z_m=0.5957
            else:
                k_c= 1.4195
                k_m=0.7489

                z_c=0.0014297
                z_m=0.6729
           
        elif self.band == 'S' or self.band == 's':
            if self.exper == 'MC3E':
                k_c= 1.4489
                k_m=0.6589

                z_c=0.0019393
                z_m=0.592032
            elif self.exper == 'TWPICE': 
                k_c= 2.559
                k_m=0.76687

                z_c=0.0014377
                z_m=0.66743
            else:
                print("No ice-mass coefficients for your project. Using MC3E relations")
                k_c= 1.4489
                k_m=0.6589
                z_c=0.0014297
                z_m=0.6729
        
        else:
            print('Problem in ice-mass. No criteria (wavelength or project)')
            #print ('Your wavelength has not been run yet! Please return to fundamentals.')
            return
            
        dbzz = np.ma.masked_less(self.data[self.dz_name].values,35)
        kdpp = np.ma.masked_where(dbzz.mask,self.data[self.kdp_name].values)

        kdppp = np.ma.masked_less(kdpp,0.2)
        dbzzz = np.ma.masked_where(kdppp.mask,dbzz)

        M_Kdp = k_c*(kdppp**k_m)
        
        kdpp =0.0
        kdppp = 0.0
        dbzz = 0.0
        dbzzz = 0.0
        
        print('Masking data')
        start = time.time()
        kdp_notmet = np.ma.masked_greater_equal(self.data[self.kdp_name].values,0.2)
        dbz_notmet = np.ma.masked_where(kdp_notmet.mask,self.data[self.dz_name].values)

        dbz_notmett = np.ma.masked_greater_equal(dbz_notmet,35.0)
        kdp_notmett = np.ma.masked_where(dbz_notmett.mask,kdp_notmet)

        kdp_notmettt = np.ma.masked_less(kdp_notmett,-100.)
        dbz_notmettt = np.ma.masked_where(kdp_notmettt.mask,dbz_notmett)

        linz = 10.**(dbz_notmettt/10.)  # for method w/o HID
        M_Z = z_c*(linz**z_m)
    
        dbz_notmett = 0.0
        kdp_notmett = 0.0
        dbz_notmettt=0.0
        kdp_notmettt=0.0
        print(time.time()-start)
        print('got Kdp')
        start = time.time() 
        M_Z.set_fill_value(0.0)
        M_Kdp.set_fill_value(0.0)
        lwccc = M_Z.filled()+M_Kdp.filled()
        qrr = lwccc/1.225
        qrr = np.ma.masked_less_equal(qrr,0.0)
        print(time.time()-start)
        print("Made through calculations, saving data")
        start = time.time()
        self.add_field((self.data[self.dz_name].dims,qrr,), 'rqr')
        print(time.time()-start)
        print('saved data')


    def calc_rr_pol(self,tm,rr_dir,band=None):

#        import pydisdrometer as pyd
#        import pytmatrix as pyt
        import csu_blended_rain_julie
        import time

        ### This is where I do mixing ratio calculations ### 
        if band == None:
            band = self.band
        
        if band == 'C' or band == 'c':
            if self.expr == 'MC3E':
                k_c=20.59
                k_m=0.75

                z_c=0.0179
                z_m=0.6855
            
                azzdr_coeff=0.00842
                bzzdr_coeff=0.894
                czzdr_coeff=-3.4029
            
                a=azdrk_coeff=33.0097
                b=bzdrk_coeff=0.9148
                c=czdrk_coeff=-0.92778
            if self.expr == 'TWPICE': 
                k_c=28.0669
                k_m=0.8049

                z_c=0.0152889
                z_m=0.738
            
                azzdr_coeff=0.008097
                bzzdr_coeff=0.87966
                czzdr_coeff=-2.6350
            
                a=azdrk_coeff=41.0249
                b=bzdrk_coeff=0.882359
                c=czdrk_coeff=-1.30948
        elif band == 'S' or band == 's':
            if self.expr == 'MC3E':
                k_c=35.77975
                k_m=0.767

                z_c=0.017850
                z_m=0.6817
            
                azzdr_coeff=0.00463
                bzzdr_coeff=0.9539
                czzdr_coeff=-3.9243
            
                a=azdrk_coeff=96.8418
                b=bzdrk_coeff=0.965216
                c=czdrk_coeff=-1.9881
            if self.expr == 'TWPICE': 
                k_c=51.7402
                k_m=0.81788

                z_c=0.015364
                z_m=0.73256
            
                azzdr_coeff=0.006482
                bzzdr_coeff=0.92738
                czzdr_coeff=-4.0552
            
                a=azdrk_coeff=94.1082
                b=bzdrk_coeff=0.9337
                c=czdrk_coeff=-1.9350
            else:
                print('Specific rain rate coefficients not avialable for your project. Using MC3E')
                k_c=35.77975
                k_m=0.767

                z_c=0.017850
                z_m=0.6817
            
                azzdr_coeff=0.00463
                bzzdr_coeff=0.9539
                czzdr_coeff=-3.9243
            
                a=azdrk_coeff=96.8418
                b=bzdrk_coeff=0.965216
                c=czdrk_coeff=-1.9881
                
        else:
            print ('Sorry, your wavelength has not been run yet! Return to first principles!')
            return
        
        start = time.time()
        rr_arr,rm = csu_blended_rain_julie.csu_hidro_rain(self.data[self.dz_name].values,self.data[self.zdr_name].values,self.data[self.kdp_name].values,z_c,z_m,k_c,k_m,azdrk_coeff,bzdrk_coeff, czdrk_coeff,azzdr_coeff,bzzdr_coeff,czzdr_coeff,band=band,fhc=self.hid)
        print(time.time()-start)
        #rr_arr,rm = csu_blended_rain_julie.calc_blended_rain(self.data[self.dz_name].values,self.data[self.zdr_name].values,self.data[self.kdp_name].values,z_c,z_m,k_c,k_m,None,None,azdrk_coeff,bzdrk_coeff, czdrk_coeff,azzdr_coeff,bzzdr_coeff,czzdr_coeff,band=band)
#         mask=np.where(np.isnan(self.data[self.dz_name].values))
#         rr[mask]=np.nan
#         rm[mask]=-1

        self.add_field((self.data[self.dz_name].dims,rr_arr,),self.rr_name)
        self.add_field((self.data[self.dz_name].dims,rm,),'RRM')
        #whbad = np.where(np.isnan(self.data[self.dz_name])) #<<< attempt to speed up by commenting
        #self.data[self.rr_name].values[whbad]=np.nan #<<<
        #self.data['RRM'].values[whbad]=-1 #<<<
        
        it = self.data[self.rr_name].shape[0]
        print('\nSaving RR to files...')
        for ii in tqdm(range(it)):
            filepath = rr_dir+'/RR_hidro_'+self.exper+'_'+str(tm[ii].strftime('%Y%m%d_%H%M%S'))+'.nc'
            if not os.path.isfile(filepath): self.data[self.rr_name].sel(d=ii).to_netcdf(path=filepath, mode='w')  
        print('')
        
        return
#############################################################################################################

    def score_xsec_plot(self, y=None, title_flag=False, *args, **kwargs):

        "Cross section showing scores for all species, used as HID analysis tool"
        # first, get the appropriate y index from the y that's wanted
        if y is None:
            y_ind = int(len(self.data[self.y_name])/2.0)
        else:
            y_ind = np.argmin(abs(y - self.data[self.y_name]))

        fig, ax = plt.subplots(5,2, figsize = (9,9))
        
        axf = ax.flatten()
        fig.subplots_adjust(top = 0.90, bottom = 0.12, left = 0.05, right = 0.92, hspace = 0.4)

        for c in range(10): # looping thru each number corresponding to each species

            dummy = axf[c].pcolormesh(self.data[self.x_name], self.data[self.z_name], self.scores[c,:,y_ind,:], cmap = plt.cm.gist_heat_r,
                     vmin = 0, vmax = 1.0)


        axf[c].set_title(self.species[c])

        cbar_ax = fig.add_axes([0.15, 0.06, 0.7, 0.03])
        cb = fig.colorbar(dummy, cax = cbar_ax, orientation = 'horizontal')
        cb.set_label(r'$\mu$'+' score')

        if title_flag:
            fig.suptitle('%s %s Cross Section HID scores' %(self.print_date(), self.band+'-band'), fontsize = 14)


        return fig, ax

#############################################################################################################

    def score_cappi_plot(self, z=1.0, title_flag=False, *args, **kwargs):

        "CAPPI plot showing scores for all species, used as HID analysis tool"

        # first, get the appropriate z index from the z that's wanted in altitude
        z_ind = np.argmin(abs(z - self.data[self.z_name]))

        fig, ax = plt.subplots(5,2, figsize = (9,9))
        
        axf = ax.flatten()
        fig.subplots_adjust(top = 0.90, bottom = 0.12, left = 0.05, right = 0.92, hspace = 0.4)

        for c in range(10): # looping thru each number corresponding to each species

            dummy = axf[c].pcolormesh(self.data[self.x_name], self.data[self.y_name], self.scores[c,z_ind,:,:], cmap = plt.cm.gist_heat_r,
                     vmin = 0, vmax = 1.0)

        axf[c].set_title(self.species[c])

        cbar_ax = fig.add_axes([0.15, 0.06, 0.7, 0.03])
        cb = fig.colorbar(dummy, cax = cbar_ax, orientation = 'horizontal')
        cb.set_label(r'$\mu$'+' score')

        if title_flag:
            fig.suptitle('%s %s CAPPI HID scores' %(self.print_date(), self.band+'-band'), fontsize = 14)


        return fig, ax


    def xsec(self, var, y=None, xlim=None, zmax=None, cbar=1, ts = None,varlist=None, ax=None, title_flag=False, vectors=None, cblabel=None, res=2.0,cbpad=0.03, labels=True, xlab=False, ylab=False, latlon=False, lblsz=16, lblpad=15, **kwargs):
        
        if ts is not None: 
            try:
                tmind = np.where(np.array(self.date) == ts)[0][0]
            except IndexError as e:
                tmind = np.where(np.array(self.date) == ts)[0]
                
        if not xlim:
            xmin, xmax = np.floor(self.data[self.x_name].values.min()), np.ceil(self.data[self.x_name].values.max())
        else:
            xmin, xmax = xlim[0], xlim[1]
        if not zmax:
            zmax = self.data[self.z_name].values.max()

        if 'y' in self.data[self.x_name]:
            xdataset = self.data[self.x_name].sel(x=slice(xmin,xmax),y=slice(ymin,ymax))
        else:
            xdataset = self.data[self.x_name].sel(x=slice(xmin,xmax))

        zdataset = self.data[self.z_name].sel(z=slice(0,zmax))

        if 'd' in xdataset:
            xdat = np.squeeze(xdataset.sel(d=tmind).values)
            zdat = np.squeeze(zdataset.sel(d=tmind).values)
        else:
            xdat = np.squeeze(xdataset.values)
            zdat = np.squeeze(zdataset.values)

        dataset = self.data[var].sel(z=slice(0,zmax),y=y,x=slice(xmin,xmax))
        data = np.squeeze(dataset.sel(d=tmind).values)
        data = np.ma.masked_where(~np.isfinite(data),data)
           
        if var.startswith('HID'):
            data = np.ma.masked_where(data < 1,data)

        if ax is None:
            fig = plt.figure(figsize=(10,8))
            ax = fig.add_subplot(111)
        else:
            fig = ax.get_figure()

        if var in self.lims.keys():
            dummy = ax.pcolormesh(xdat, zdat, data, vmin = self.lims[var][0], vmax = self.lims[var][1], cmap = self.cmaps[var], **kwargs)
        else:
            dat = self.data[var].values
            dat[dat<-900.0]=np.nan
            range_lim  = np.nanmax(dat) - np.nanmin(dat)
            dummy = ax.pcolormesh(xdat,zdat, data,
                vmin = np.nanmin(dat), vmax = np.nanmax(dat),cmap = plt.cm.gist_ncar, **kwargs)
        
        ####### plotting limits getting set here ######
        if self.x_name == 'longitude':
            if labels:
                ax.set_xlabel('Distance E of Radar (km)', fontsize=lblsz)
                ax.set_ylabel('Altitude (km MSL)', fontsize=lblsz)
                ax.set_xlim([xmin,xmax])
                ax.tick_params(axis='both', which='major', labelsize=lblsz)
            else:
                if xlab:
                    ax.set_xlabel('Distance E of Radar (km)', fontsize=lblsz)
                    ax.tick_params(axis='x', which='major', labelsize=lblsz)
                else:
                    ax.set_xticklabels([])
                    ax.xaxis.set_ticks_position('none')
                if ylab:
                    ax.set_ylabel('Altitude (km MSL)', fontsize=lblsz)
                    ax.tick_params(axis='y', which='major', labelsize=lblsz)
                else:
                    ax.set_yticklabels([])
                    ax.yaxis.set_ticks_position('none')
        
        else:
            if labels:
                ax.set_xlabel('Distance E of radar (km)',fontsize=lblsz)
                ax.set_ylabel('Altitude (km MSL)',fontsize=lblsz)
                ax.tick_params(axis='both', which='major', labelsize=lblsz)
            else:
                if xlab:
                    ax.set_xlabel('Distance E of radar (km)',fontsize=lblsz)
                    ax.tick_params(axis='x', which='major', labelsize=lblsz)
                else: 
                    ax.set_xticklabels([])
                    ax.xaxis.set_ticks_position('none')
                if ylab:
                    ax.set_ylabel('Altitude (km MSL)',fontsize=lblsz)
                    ax.tick_params(axis='y', which='major', labelsize=lblsz)
                else:
                    ax.set_yticklabels([])
                    ax.yaxis.set_ticks_position('none')
 
        ax.set_xlim([xmin,xmax])
        ax.set_xticks(np.linspace(xmin,xmax,5))
        ax.set_ylim([0,zmax]) 
        ax.set_yticks(np.linspace(0,zmax,6))
        ax.grid(color='grey', linestyle='-', linewidth=1)
    
        if cbar == 1:  # call separate HID colorbar function for bar plots
            if var.startswith('HID'): 
                lur,bur,wur,hur = ax.get_position().bounds
                cbar_ax_dims = [lur+wur+0.015,bur,0.03,hur]
                cbt = self.HID_barplot_colorbar(fig,cbar_ax_dims)
                cbt.ax.tick_params(labelsize=16)
                cbt.set_label(self.names_uc[var]+' '+self.units[var], fontsize=16, rotation=270, labelpad=20)
            else:
                self.mycbar(fig,ax,dummy,self.longnames[var]+' '+self.units[var])

        if cbar == 2 and var.startswith('HID'):
            lur,bur,wur,hur = ax.get_position().bounds
            cbar_ax_dims = [lur,bur-0.125,wur,0.03]
            self.HID_barplot_colorbar(fig,cbar_ax_dims,orientation='horizontal',names='longnames', lblsz=12)
 
        if vectors:
#                 try:
#                     print( zmax,xlim,ts,res)
                self.xsec_vector(ax=ax, y=y,zmax=zmax,xlim=xlim,ts=ts,res=res)
#                 except Exception as e:
#                     print ('Error trying to plot xsec vectors: {}'.format(e))

        if title_flag:
            ax.set_title('%s %s Cross Section' %(ts, self.band+'-band'), fontsize = 14)
        else:
            dummy = fig
#        print type(dummy),dummy

        return dummy, ax

    
    def xsec_multiplot(self, y=0.5, xlim=[], zmax=[], ts=None, res = 2.0, varlist=None, vectors=None, latlon=False, **kwargs):
    
        if ts is not None:
            try:
                tmind = np.where(np.array(self.date)==ts)[0][0]
            except:
                tmind = np.where(np.array(self.date)==ts)[0]

        if varlist is not None:
            good_vars = varlist
        else:
            good_vars = self.valid_vars()

        nvars = len(good_vars)
        if 'scores' in good_vars:
            print('needing to add + to HID, ln 1257 RadarData')
            if hasattr(self, 'scores'):
                nvars += 1
        if nvars <= 3:
            ncols = 1
            nrows = deepcopy(nvars)
            figx = 5
            figy = 3*nrows
        elif (nvars > 3 and nvars < 7):
            nrows = 2
            ncols = int(np.ceil(nvars/2))
        else:
            ncols = 2
            nrows = int(np.ceil(nvars/2))
            figx= 18
            figy = 18

        if 'HID' in good_vars:
            fig, ax = plt.subplots(nrows,ncols,figsize=(16,8),gridspec_kw={'wspace': 0.45, 'hspace': 0.07, 'top': 1., 'bottom': 0., 'left': 0., 'right': 1.})
        else:
            fig, ax = plt.subplots(nrows,ncols,figsize=(14,8),gridspec_kw={'wspace': 0.3, 'hspace': 0.07, 'top': 1., 'bottom': 0., 'left': 0., 'right': 1.})

        if not isinstance(ax, np.ndarray) or not isinstance(ax, list): 
            ax = np.array([ax])
        axf = ax.flatten()

        for i in range(len(good_vars),ncols*nrows):
            fig.delaxes(axf[i])
 
        for i, var in enumerate(good_vars):
            if vectors is not None:
                vect = vectors[i]
            else:
                vect = None
            botpanels = np.arange(nvars-ncols,nvars)
            xlabbool = True if i in botpanels else False
            lspanels = [ncols*n for n in range(0,nrows)]
            ylabbool = True if i in lspanels else False
            #dummy = self.xsec(var, ts=ts, y=y, vectors=vect, xlim=xlim, zmax=zmax, ax=axf[i],res=res,xlab=xlabbool,ylab=ylabbool,labels=False,lblsz=28,lblpad=28,**kwargs)
            dummy = self.xsec(var, ts=ts, y=y, xlim=xlim, zmax=zmax, ax=axf[i],res=res,xlab=xlabbool,ylab=ylabbool,labels=False,latlon=latlon,**kwargs)

        axf[0].text(0, 1, '{e} {r}'.format(e=self.exper,r=self.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=20, color='k', zorder=10, weight='bold', transform=axf[0].transAxes) # (a) Top-left

        axf[ncols-1].text(1, 1, '{d:%Y-%m-%d %H:%M:%S} UTC'.format(d=ts), horizontalalignment='right', verticalalignment='bottom', size=20, color='k', zorder=10, weight='bold', transform=axf[ncols-1].transAxes) # (a) Top-left
        
        axf[ncols-1].text(0.99, 0.99, 'y = {a} km'.format(a=y), horizontalalignment='right',verticalalignment='top', size=20, color='k', zorder=10, weight='bold', transform=axf[ncols-1].transAxes, bbox=dict(facecolor='w', edgecolor='none', pad=0.0))
 
        return fig

    def get_ind(self,val,dat):
        dum = abs(val - dat)
        wh_t = np.squeeze(np.where(dum == np.min(dum)))
        try:
            t = (len(wh_t))
            if t==1:
                return wh_t[0]
            elif t==2:
                try:
                    return wh_t[0][0]
                except:
                    return wh_t[0]
            elif t==3:
                return wh_t[0][0][0]
            else:
                return wh_t
        except:
                 return wh_t
    
######################### Here is the 4 stuff ##############################

    def cappi(self, var, z=1.0, xlim=[], ylim=[], latlon=False, ax=None, ts = None, title_flag=False, vectors=None, cblabel=None, labels=True, xlab=1, ylab=1, cbar=1, res = 2.0, thresh_dz=False,contours = None,statpt=False, dattype='obs', **kwargs):

        from copy import deepcopy
        import cartopy.crs as ccrs
        import xarray as xr

        if ts is not None:
            try:
                tmind = np.where(np.array(self.date)==ts)[0][0]
            except IndexError as e:
                tmind = np.where(np.array(self.date)==ts)[0]
     
        if latlon:

            if 'd' in self.data[self.lon_name].dims:
                lons = self.data[self.lon_name].sel(d=tmind).values
                lats = self.data[self.lat_name].sel(d=tmind).values
            else:
                lons = self.data[self.lon_name].values
                lats = self.data[self.lat_name].values
            
            if not xlim:
                xmin, xmax = lons[(0,0)], lons[(0,-1)]
                ymin, ymax = lats[(0,0)], lats[(-1,0)]
            else:
                xmin, xmax = np.min(xlim), np.max(xlim)
                ymin, ymax = np.min(ylim), np.max(ylim)
            
            xshift = 0.5*(self.data[self.x_name].shape[0]-1)
            xmin = np.where(lons[0,:] >= xmin)[0][0]-xshift
            xmax = np.where(lons[0,:] <= xmax)[0][-1]-xshift

            yshift = 0.5*(self.data[self.y_name].shape[0]-1)
            ymin = np.where(lats[:,0] >= ymin)[0][0]-yshift
            ymax = np.where(lats[:,0] <= ymax)[0][-1]-yshift
            
            xdataset = self.data[self.lon_name].sel(x=slice(xmin,xmax),y=slice(ymin,ymax))
            ydataset = self.data[self.lat_name].sel(x=slice(xmin,xmax),y=slice(ymin,ymax))
        
        else:

            if not xlim: xmin, xmax = self.data[self.x_name].values.min(), self.data[self.x_name].values.max()
            else: xmin, xmax = np.min(xlim), np.max(xlim)
            
            if not ylim: ymin, ymax = self.data[self.x_name].values.min(), self.data[self.x_name].values.max()
            else: ymin, ymax = np.min(ylim), np.max(ylim)

            if 'y' in self.data[self.x_name].dims: xdataset = self.data[self.x_name].sel(x=slice(xmin,xmax),y=slice(ymin,ymax))
            else: xdataset = self.data[self.x_name].sel(x=slice(xmin,xmax))

            if 'x' in self.data[self.y_name].dims: ydataset = self.data[self.y_name].sel(x=slice(xmin,xmax),y=slice(ymin,ymax))
            else: ydataset = self.data[self.y_name].sel(y=slice(ymin,ymax))
                
        if 'd' in xdataset.dims:
            xdat = np.squeeze(xdataset.sel(d=tmind).values)
            ydat = np.squeeze(ydataset.sel(d=tmind).values)
        else:
            xdat = np.squeeze(xdataset.values)
            ydat = np.squeeze(ydataset.values)
           
        dataset = self.data[var].sel(z=z,x=slice(xmin,xmax),y=slice(ymin,ymax))
        
        data = np.squeeze(dataset.sel(d=tmind).values)
        data = np.ma.masked_where(~np.isfinite(data),data)
        if var.startswith('HID'): 
            data = np.ma.masked_where(data < 1,data)
        elif var.startswith('RR'):
            data = np.ma.masked_where(data < 0.5,data)

        if ax is None:
            if latlon: 
                fig = plt.figure(figsize=(10,8))
                ax = fig.add_subplot(111, projection=ccrs.Mercator())
            else:
                fig = plt.figure(figsize=(10,8))
                ax = fig.add_subplot(111)      
        else:
            fig = ax.get_figure()
        
        if var in self.lims.keys():
            if latlon:
                dummy = ax.pcolormesh(xdat,ydat, data, vmin = self.lims[var][0], vmax = self.lims[var][1], cmap = self.cmaps[var], transform=ccrs.PlateCarree(), **kwargs)
            else:
                dummy = ax.pcolormesh(xdat,ydat, data, vmin = self.lims[var][0], vmax = self.lims[var][1], cmap = self.cmaps[var], **kwargs)
                
        else:
            dat = self.data[var].data
            dat[dat<-900.0]=np.nan
            range_lim  = np.nanmax(dat) - np.nanmin(dat)
            dummy = ax.pcolormesh(xdat,ydat, data,
                vmin = np.nanmin(dat), vmax = np.nanmax(dat), cmap = plt.cm.gist_ncar,**kwargs)

        if contours is not None:
            if contours == 'CS':
                csvals = np.squeeze(self.data[self.cs_name].sel(d=tmind,z=z,x=slice(xmin,xmax),y=slice(ymin,ymax)).values)
                if latlon:
                    cb = ax.contour(xdat, ydat, csvals, levels=[1,2], colors=['k'], linewidths=[3], alpha=0.8, zorder=10, transform=ccrs.PlateCarree())
                else:
                    cb = ax.contour(xdat, ydat, csvals, levels=[1,2], colors=['k'], linewidths=[3], alpha=0.8, zorder=10)
        
        if not xlim:
          
            maxxs = []
            minxs = []
            maxys = []
            minys = []
            
            for ii in range(dataset.shape[0]):

                if 'd' in xdataset.dims:
                    xdat_masked = deepcopy(xdataset.sel(d=ii).values)
                    ydat_masked = deepcopy(ydataset.sel(d=ii).values)
                else:
                    xdat_masked = deepcopy(xdataset.values)
                    ydat_masked = deepcopy(ydataset.values)

            if not latlon and xdat_masked.ndim == 1: xdat_masked,ydat_masked = np.meshgrid(xdat_masked,ydat_masked) 
            
            xdat_masked = np.where(np.isnan(dataset.sel(d=ii).values),np.nan,xdat_masked)
            ydat_masked = np.where(np.isnan(dataset.sel(d=ii).values),np.nan,ydat_masked)

            minxs.append(np.nanmin(xdat_masked))
            maxxs.append(np.nanmax(xdat_masked))
            minys.append(np.nanmin(ydat_masked))
            maxys.append(np.nanmax(ydat_masked))
            
            if latlon:
            
                minx = np.round(np.min(minxs)-0.1,1)
                maxx = np.round(self.lon_0+(self.lon_0-np.min(minxs))+0.1,1)
                miny = np.round(np.min(minys)-0.1,1)
                maxy = np.round(np.max(maxys)+0.1,1)
                
            else:
                
                minx = np.round(np.min(minxs),1)
                maxx = np.round(np.abs(np.min(minxs)),1)
                miny = np.round(np.min(minys),1)
                maxy = np.round(np.max(maxys),1)
 
        else:

            minx = np.min(xlim)
            maxx = np.max(xlim)
            miny = np.min(ylim)
            maxy = np.max(ylim)
           
        self.co_gridlines(fig,ax,latlonyn=latlon,minx=minx,maxx=maxx,miny=miny,maxy=maxy,xlab=xlab,ylab=ylab)
        
        if statpt: 
            if latlon: ax.plot(self.lon_0,self.lat_0,markersize=16,marker='^',color='k',transform=ccrs.PlateCarree())
            else: ax.plot(0,0,markersize=16,marker='^',color='k')

        if cbar == 1:
            if var.startswith('HID'): 
                lur,bur,wur,hur = ax.get_position().bounds
                cbar_ax_dims = [lur+wur+0.015,bur,0.03,hur]
                cbt = self.HID_barplot_colorbar(fig,cbar_ax_dims) 
                cbt.ax.tick_params(labelsize=16)
                cbt.set_label(self.names_uc[var]+' '+self.units[var], fontsize=16, rotation=270, labelpad=20)
            else:
                self.mycbar(fig,ax,dummy,self.longnames[var]+' '+self.units[var])
        
        elif cbar == 2 and var.startswith('HID'):
            lur,bur,wur,hur = ax.get_position().bounds
            if latlon: cbar_ax_dims = [lur,bur-0.075,wur,0.03]
            else: cbar_ax_dims = [lur,bur-0.125,wur,0.03]
            self.HID_barplot_colorbar(fig,cbar_ax_dims,orientation='horizontal',names='longnames', lblsz=14)
 
        # Now check for the vectors flag, if it's there then plot it over the radar stuff
        if vectors is not None:
#            try:
#                print 'RadarDAta 1177:', res,z,ax
                self.plan_vector(ax=ax, ts=ts,z=z,res=res,thresh_dz=thresh_dz,xlim=xlim,ylim=ylim)
#            except Exception, e:
#                print 'Error trying to plot vectors: {}'.format(e)
        if 'd' in self.data[self.z_name].dims:
            hts = self.data[self.z_name].sel(d=tmind).values
        else:
            hts = self.data[self.z_name].values

        if title_flag:
            ax.set_title('%s %s CAPPI %.1f km MSL' %(ts, self.band+'-band', \
                    hts[z_ind]), fontsize = 14)
#        print type(dummy),dummy
        return dummy,ax

#############################################################################################################

    def cappi_multiplot(self, z=1.0, xlim=[], ylim=[], ts=None,res = 2, varlist=None, vectors=None, contours = None,thresh_dz = False, latlon=False, statpt=False, dattype='obs', **kwargs):

        import cartopy.crs as ccrs

        if ts is not None:
            try:
                tmind = np.where(np.array(self.date)==ts)[0][0]
            except:    
                tmind = np.where(np.array(self.date)==ts)[0]

        if varlist is not None:
            good_vars = varlist
        else:
            good_vars = self.valid_vars()
            
        nvars = len(good_vars)
        if 'scores' in good_vars:
            if hasattr(self, 'scores'):
                nvars += 1
        if nvars <= 3:
            ncols = 1
            nrows = deepcopy(nvars)
            figx = 5
            figy = 3*nrows
        else:
            nrows = 2
            ncols = int(np.ceil(nvars/2))

        if latlon:
            fig, ax = plt.subplots(nrows,ncols,figsize=(16,8),subplot_kw={'projection': ccrs.Mercator()},gridspec_kw={'wspace': 0.45, 'hspace': 0.075, 'top': 1., 'bottom': 0., 'left': 0., 'right': 1.})
        else:
            fig, ax = plt.subplots(nrows,ncols,figsize=(16,8),gridspec_kw={'wspace': 0.45, 'hspace': 0.07, 'top': 1., 'bottom': 0., 'left': 0., 'right': 1.})
        
        if not isinstance(ax, np.ndarray) or not isinstance(ax, list): 
            ax = np.array([ax], **kwargs)
        axf = ax.flatten()
        
        for i, var in enumerate((good_vars)):
            if contours is not None:
                vcont = contours[i]
            else:
                vcont = None
            if vectors is not None:
                vect = vectors[i]
            else:
                vect = None
            botpanels = np.arange(nvars-ncols,nvars)
            xlabbool = True if i in botpanels else False
            lspanels = [ncols*n for n in range(0,nrows)]
            ylabbool = True if i in lspanels else False
            #dummy = self.cappi(var, z=z, ax=axf[i], xlim=xlim, ylim=ylim,ts = ts, vectors=vect,res=res,contour=vcont,thresh_dz =thresh_dz,xlab=xlabbool,ylab=ylabbool,cbar=1,labels=False,statpt=statpt,latlon=False)
            dummy = self.cappi(var, z=z, ax=axf[i], xlim=xlim, ylim=ylim,ts = ts, res=res, thresh_dz=thresh_dz,xlab=xlabbool,ylab=ylabbool,cbar=1,labels=False,statpt=statpt,latlon=latlon,dattype=dattype)
        
        for i in range(len(good_vars),ncols*nrows):
            fig.delaxes(axf[i])
         
        axf[0].text(0, 1, '{e} {r}'.format(e=self.exper,r=self.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=20, color='k', zorder=10, weight='bold', transform=axf[0].transAxes) # (a) Top-left
        
        axf[ncols-1].text(1, 1, '{d:%Y-%m-%d %H:%M:%S} UTC'.format(d=ts), horizontalalignment='right', verticalalignment='bottom', size=20, color='k', zorder=10, weight='bold', transform=axf[ncols-1].transAxes) # (a) Top-left
        
        axf[ncols-1].text(0.99, 0.99, 'z = {a} km'.format(a=z), horizontalalignment='right',verticalalignment='top', size=20, color='k', zorder=10, weight='bold', transform=axf[ncols-1].transAxes,bbox=dict(facecolor='w', edgecolor='none', pad=0.0))
        
        return fig


#############################################################################################################
    # Down here is dual doppler plotting stuff

    def xsec_vector(self, y=None, xlim=None,zmax=None,ts=None,ax=None, res=2.0, ht_offset=0.2, **kwargs):
        if ts is None:
            print('xsec_vector got no time')
            ts=self.date[0]
            tsi = 0
        try:
            tmind = np.where(np.array(self.date)==ts)[0][0]
        except:    
            tmind = np.where(np.array(self.date)==ts)[0]
#         else:
#             print 'Check your dates!', ts
#        print 'ts:',ts

        if zmax is None:
            if 'd' in self.data[self.z_name].dims:
                zmax=len(self.data[self.z_name].sel(d=0))
            else:
                zmax=len(self.data[self.z_name])
            
        if 'd' in self.data[self.z_name].dims:
#                print('getting z-ind')
            zmini = self.get_ind(0,np.squeeze(self.data[self.z_name].sel(d=tmind).values))
            zmaxi = self.get_ind(zmax,np.squeeze(self.data[self.z_name].sel(d=tmind).values))
        else:
            zmini = self.get_ind(0,np.squeeze(self.data[self.z_name].values))
            zmaxi = self.get_ind(zmax,np.squeeze(self.data[self.z_name].values))

#        print('xlims 1203',xlim,tmind)
        if self.x_name == 'longitude':
            if 'd' in self.data[self.x_name].dims:
                if 'y' in self.data[self.x_name].dims:
                    xmini = self.get_ind(xlim[0],self.data[self.x_name].sel(d=tmind,y=0).values)
                    xmaxi = self.get_ind(xlim[1],self.data[self.x_name].sel(d=tmind,y=0).values)
                else:
                    xmini = self.get_ind(xlim[0],self.data[self.x_name].sel(d=tmind).values)
                    xmaxi = self.get_ind(xlim[1],self.data[self.x_name].sel(d=tmind).values)
            
            else:                
                if 'y' in self.data[self.x_name].dims:
                    xmini = self.get_ind(xlim[0],self.data[self.x_name].sel(y=0).values)
                    xmaxi = self.get_ind(xlim[1],self.data[self.x_name].sel(y=0).values)
                else:
                    xmini = self.get_ind(xlim[0],self.data[self.x_name].values)
                    xmaxi = self.get_ind(xlim[1],self.data[self.x_name].values)
        else:
            xmini, xmaxi = xlim
        
        xmin, xmax = xlim
        if np.size(res) == 1:
            resx = res
            resz = res
        else:
            resx = res[0]
            resz = res[1]
        
#         if self.data[self.x_name].units == "[deg]":
#             skip = int(np.round(resx/self.dx))
#             xskip = int(np.round(skip)#/110.))
# #            print '1245', xskip
#             
#         else:
        skip = int(np.round(resx/self.dx))
        xskip = skip
        if self.data[self.z_name].units == "[deg]":
            skip = int(np.round(resz/self.dz))
            zskip = skip/110.
        else:
            skip = int(np.round(resz/self.dz))
            zskip = 1
#        print('resx, dx',resx,self.dx)
#        print 'zskip 1257',zskip
            
        #print skip,xskip, zskip
        #print skip

    # first, get the appropriate y index from the y that's wanted
#        print('y is ',y)
        if y is None:
            y_ind = int(len(self.data[self.y_name].data[0])/2.0)
        else:
            if 'd' in self.data[self.y_name].dims:
                if 'x' in self.data[self.y_name].dims:
                    y_ind = self.get_ind(y,self.data[self.y_name].sel(d=tmind,x=0).values)
                else:
                    y_ind = self.get_ind(y,self.data[self.y_name].sel(d=tmind).values)
            else:
                if 'x' in self.data[self.y_name].dims:
                    y_ind = self.get_ind(y,self.data[self.y_name].sel(x=0).values)
                else:
                    y_ind = self.get_ind(y,self.data[self.y_name].values)            

#        print 'y in xsec', y

        if ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig = ax.get_figure()

#        print 'ln 1274', xmini,xmaxi,y_ind,zmini,zmaxi,skip,xlim,y

        if self.u_name in self.data.variables.keys():

            try:
                if self.y_name == 'latitude':
                    if 'd' in self.data[self.u_name].dims:
                        udat= np.squeeze(np.squeeze(self.data[self.u_name]).sel(d=tmind,x=slice(xmini,xmaxi),y=y_ind,z=slice(zmini,zmaxi)).values)
                        wdat= np.squeeze(np.squeeze(self.data[self.w_name]).sel(d=tmind,x=slice(xmini,xmaxi),y=y_ind,z=slice(zmini,zmaxi)).values)
                        xdat= np.squeeze(np.squeeze(self.data[self.x_name]).sel(d=tmind,x=slice(xmini,xmaxi)).values)

                        if 'y' in self.data[self.z_name].dims:
                        
                            zdat= np.squeeze(np.squeeze(self.data[self.z_name]).sel(d=tmind,z=slice(zmini,zmaxi),y=y_ind).values)
                        else:
            
                            zdat= np.squeeze(np.squeeze(self.data[self.z_name]).sel(d=tmind,z=slice(zmini,zmaxi)).values)

                    else:
                        udat= np.squeeze(np.squeeze(self.data[self.u_name]).sel(x=slice(xmini,xmaxi),y=slice(y,y+1),z=slice(zmini,zmaxi)).values)
                        wdat= np.squeeze(np.squeeze(self.data[self.w_name]).sel(x=slice(xmini,xmaxi),y=slice(y,y+1),z=slice(zmini,zmaxi)).values)
                        xdat= np.squeeze(np.squeeze(self.data[self.x_name]).sel(x=slice(xmini,xmaxi)).values)
                        if 'y' in self.data[self.z_name].dims:
                            zdat= np.squeeze(np.squeeze(self.data[self.z_name]).sel(z=slice(zmini,zmaxi),y=y_ind).values)
                        else:
                            zdat= np.squeeze(np.squeeze(self.data[self.z_name]).sel(z=slice(zmini,zmaxi)).values)
                        
                        

                else:
    #                print np.shape(xdat), np.shape(zdat)
                    #xdat = np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi+1),y=y_ind).values) 
                    if 'd' in self.data[self.u_name].dims:
#                        print('made line 1936')
                        xdat = np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi+1)).values)
                        zdat = np.squeeze(self.data[self.z_name].sel(z=slice(zmini,zmaxi+1)).values)
                        #udat = np.squeeze(self.data[self.u_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=y_ind).values)
                        udat = np.squeeze(self.data[self.u_name].sel(d=tmind,z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y,y+1)).values)
                        #wdat = np.squeeze(self.data[self.w_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=y_ind).values)
                        wdat = np.squeeze(self.data[self.w_name].sel(d=tmind,z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y,y+1)).values)                  
                        if len(np.shape(udat))>2:
#                            print('Udat shape is too many!')
                            udat = np.squeeze(self.data[self.u_name].sel(d=tmind,z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y,y)).values)
                            #wdat = np.squeeze(self.data[self.w_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=y_ind).values)
                            wdat = np.squeeze(self.data[self.w_name].sel(d=tmind,z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y,y)).values)
                    else: 
#                        print('made 1949')
                        xdat = np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi+1)).values)
                        zdat = np.squeeze(self.data[self.z_name].sel(z=slice(zmini,zmaxi+1)).values)
                        #udat = np.squeeze(self.data[self.u_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=y_ind).values)
                        udat = np.squeeze(self.data[self.u_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y,y+1)).values)
                        #wdat = np.squeeze(self.data[self.w_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=y_ind).values)
                        wdat = np.squeeze(self.data[self.w_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y,y+1)).values)

            except:
#                print 'uh-oh, exception'
                xdat = np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi+1)).values)
                zdat = np.squeeze(self.data[self.z_name].sel(z=slice(zmini,zmaxi+1)).values)
    #            print np.shape(xdat),np.shape(zdat),np.shape(self.data[self.u_name].data)
                #udat = np.squeeze(self.data[self.u_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=y_ind).values)
                udat = np.squeeze(self.data[self.u_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y,y+1)).values)
                #wdat = np.squeeze(self.data[self.w_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=y_ind).values)
                wdat = np.squeeze(self.data[self.w_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y,y+1)).values)
                if len(np.shape(udat))>2:
 #                   print('Udat shape is too many!')
                    udat = np.squeeze(self.data[self.u_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y,y)).values)
                    #wdat = np.squeeze(self.data[self.w_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=y_ind).values)
                    wdat = np.squeeze(self.data[self.w_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y,y)).values)
                
#            print('w is',np.nanmax(wdat[::zskip,::xskip]))

            if self.y_name == 'latitude':
                q_handle = ax.quiver(xdat[::xskip], zdat[::zskip]+ht_offset, \
                    udat[::zskip, ::xskip], wdat[::zskip, ::xskip], \
                        scale=60, scale_units='inches', pivot='middle', width=0.0025, headwidth=6, **kwargs)
            else:
                q_handle = ax.quiver(xdat[::xskip], zdat[::zskip]+ht_offset, \
                    udat[::zskip, ::xskip], wdat[::zskip, ::xskip], \
                        scale=70, scale_units='inches', pivot='middle', width=0.0025, headwidth=6, **kwargs)

            qk = ax.quiverkey(q_handle, 0.85, 0.85, 20, r'20 $\frac{m}{s}$',coordinates='axes', \
                        fontproperties={'weight':'bold','size':14})
        else:
            print ('No vectors. No {t}'.format(t=self.u_name))
            q_handle = fig
            return q_handle
        return q_handle

#############################################################################################################

    def plan_vector(self, z=1.0, ax=None, xlim=None,ylim = None,ts=None,res=2.0, ht_offset=0.2,thresh_dz=False, **kwargs):


        if np.size(res) == 1:
            resx = res
            resy = res
        else:
            resx = res[0]
            resy = res[1]

#             xskip = int(np.round(resx/(self.dx*110.)))
#             yskip = int(np.round(resy/(self.dy*110.)))
#             #yskip = xskip
#             #yskip = int(np.round(res/(self.dy*110.)))
#             #skip = int(np.round(skip/110.))
#         else:
        xskip = int(np.round(resx/self.dx))
        yskip = int(np.round(resy/self.dy))
#        print xskip,yskip
        if ts is not None:
            tmind = np.where(np.array(self.date)==ts)[0][0]
        else:
            tmind = 0
        
        if z is None:
            z_ind = 2
            
        else:
            if 'd' in self.data[self.z_name].dims:
#                print('getting z-ind')
                z_ind = self.get_ind(z,np.squeeze(self.data[self.z_name].sel(d=tmind).values))
            else:
                z_ind = self.get_ind(z,self.data[self.z_name].values)

#        print('xlims 1203',xlim,tmind)
        if self.x_name == 'longitude':
            if 'd' in self.data[self.x_name].dims:
                if 'y' in self.data[self.x_name].dims:
                    xmini = self.get_ind(xlim[0],self.data[self.x_name].sel(d=tmind,y=0).values)
                    xmaxi = self.get_ind(xlim[1],self.data[self.x_name].sel(d=tmind,y=0).values)
                else:
                    xmini = self.get_ind(xlim[0],self.data[self.x_name].sel(d=tmind).values)
                    xmaxi = self.get_ind(xlim[1],self.data[self.x_name].sel(d=tmind).values)
            
            else:                
                if 'y' in self.data[self.x_name].dims:
                    xmini = self.get_ind(xlim[0],self.data[self.x_name].sel(y=0).values)
                    xmaxi = self.get_ind(xlim[1],self.data[self.x_name].sel(y=0).values)
                else:
                    xmini = self.get_ind(xlim[0],self.data[self.x_name].values)
                    xmaxi = self.get_ind(xlim[1],self.data[self.x_name].values)
        else:
        
            xmini, xmaxi = xlim
        
        xmin, xmax = xlim

        if ylim is None:
 #           print 'ylim is None!'
            if self.y_name == 'latitude':
                ymint, ymaxt = self.data[self.y_name].values.min(), self.data[self.y_name].values.max()
                ymini = self.get_ind(ymint,self.data[self.y_name].sel(d=tmind).values[:,0])
                ymaxi = self.get_ind(ymaxt,self.data[self.y_name].sel(d=tmind).values[:,0])
                ymin =ymint
                ymax = ymaxt

            else:
                ymini, ymaxi = self.data[self.y_name].values.min(), self.data[self.y_name].values.max()
        else:
            if self.y_name == 'latitude':
#                print 'trying to get indices'
                ymini = self.get_ind(ylim[0],self.data[self.y_name].sel(d=tmind).values[:,0])
                ymaxi = self.get_ind(ylim[1],self.data[self.y_name].sel(d=tmind).values[:,0])
                ymin = ylim[0]
                ymax = ylim[1]

            else:
                ymini, ymaxi = ylim
                ymin = ylim[0]
                ymax = ylim[1]
             
        if ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig = ax.get_figure()
#        print xmini,xmaxi,ymini,ymaxi

        if self.u_name in self.data.variables.keys():
            if 'd' in self.data[self.x_name].dims:
                if 'y' in self.data[self.x_name].dims:
                    xdat = np.squeeze(np.squeeze(self.data[self.x_name].sel(d=tmind,x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).data))
                    ydat = np.squeeze(np.squeeze(self.data[self.y_name].sel(d=tmind,y=slice(ymini,ymaxi+1),x=slice(xmini,xmaxi+1)).data))
                
                else:
                    xdat = np.squeeze(np.squeeze(self.data[self.x_name].sel(d=tmind,x=slice(xmini,xmaxi+1)).data))
                    ydat = np.squeeze(np.squeeze(self.data[self.y_name].sel(d=tmind,y=slice(ymini,ymaxi+1)).data))
                
            else:
                if 'y' in self.data[self.x_name].dims:
                    xdat = np.squeeze(np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).data))
                    ydat = np.squeeze(np.squeeze(self.data[self.y_name].sel(y=slice(ymini,ymaxi+1),x=slice(xmini,xmaxi+1)).data))
                else:
#                    print(np.shape(self.data['x']), self.data['x'].dims)
#                    print(np.shape(self.data['y']))
                    xdat = np.squeeze(np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi+1)).data))
                    ydat = np.squeeze(np.squeeze(self.data[self.y_name].sel(y=slice(ymini,ymaxi+1)).data))
#                    print(np.shape(xdat),np.shape(ydat))
            if 'd' in self.data[self.u_name].dims:
                #udat = np.squeeze(np.squeeze(self.data[self.u_name].sel(d=tmind,z=z_ind,x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).values))
                #vdat = np.squeeze(np.squeeze(self.data[self.v_name].sel(d=tmind,z=z_ind,x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).values))
                udat = np.squeeze(np.squeeze(self.data[self.u_name].sel(d=tmind,z=slice(z_ind,z_ind+1),x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).values))
                vdat = np.squeeze(np.squeeze(self.data[self.v_name].sel(d=tmind,z=slice(z_ind,z_ind+1),x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).values))
                if np.shape(udat)[0] == 2:
                    udat = udat[0]
                    vdat = vdat[0]
                
                #print(np.shape(udat),np.shape(vdat))
                
            else:
                #udat = np.squeeze(np.squeeze(self.data[self.u_name].sel(z=z_ind,x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).values))
                #vdat = np.squeeze(np.squeeze(self.data[self.v_name].sel(z=z_ind,x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).values)) 
                udat = np.squeeze(np.squeeze(self.data[self.u_name].sel(z=slice(z_ind,z_ind+1),x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).values))
                vdat = np.squeeze(np.squeeze(self.data[self.v_name].sel(z=slice(z_ind,z_ind+1),x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).values)) 
            
            #z_ind = self.get_ind(z_ind,self.data[self.z_name].data)

            if thresh_dz == True:
                dzdat = np.squeeze(self.data[self.dz_name].sel(d=tmind,z=z_ind,x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).values)
                print ('trying to threshold...',np.shape(vdat),np.shape(dzdat))
                msk = np.less(dzdat,self.z_thresh)
                msk2 =np.where(np.logical_or(udat< -1000.,vdat < -1000.))

                vdat[msk] = np.nan
                udat[msk] = np.nan
#                xdat= np.ma.masked_where(msk,xdat)
                vdat[msk2] = np.nan
                udat[msk2] = np.nan
#                ydat= np.ma.masked_where(msk,ydat)
                #print type(vdat)
            #print np.max(vdat)
            #print 'vect shp',np.shape(udat),np.shape(vdat),np.shape(xdat),min(ydat),max(ydat)
#            print type(xdat),type(ydat),type(udat)
#            print('yskip','xskip',yskip,xskip)
            if xdat.ndim > 1:
#                print('xdat.ndim is > 1')
                try:
                    print(np.shape(xdat))
                    xdatskip = xdat[::yskip,::xskip]
                    ydatskip = ydat[::yskip,::xskip]
                    udatskip = udat[::yskip,::xskip]
                    vdatskip = vdat[::yskip,::xskip]
                except:
                    xdatskip = xdat[::yskip][::xskip]
                    ydatskip = ydat[::yskip][::xskip]
                    udatskip = udat[::yskip][::xskip]
                    vdatskip = vdat[::yskip][::xskip]

            else:
                xdatskip = xdat[::yskip]#,::xskip]
                ydatskip = ydat[::yskip]#,::xskip]
                udatskip = np.squeeze(udat[::yskip,::xskip])
                vdatskip = np.squeeze(vdat[::yskip,::xskip])
                xdatskipmesh,ydatskipmesh = np.meshgrid(xdatskip,ydatskip)
#            print np.shape(xdatskip),np.shape(ydatskip),np.shape(udatskip),np.shape(vdatskip)
#            print ('RadarData 1516:', xskip, yskip,np.shape(xdat),np.shape(ydat),np.shape(udat),np.shape(vdat))
#            print('RadarData line 2111',np.shape(xdatskipmesh),np.shape(ydatskipmesh),np.shape(udatskip),np.shape(vdatskip))
            q_handle = ax.quiver(xdatskipmesh, ydatskipmesh, \
                udatskip, vdatskip, \
                    scale=100, scale_units='inches', pivot='middle', width=0.0025, headwidth=4, **kwargs)

            qk = ax.quiverkey(q_handle, 0.08, 0.05, 20, r'20 $\frac{m}{s}$',coordinates='axes', \
                        fontproperties={'weight':'bold','size':14})
        else:
            return
        return q_handle

#############################################################################################################


    def cfad(self, var, value_bins=None, above=None, below=15.0, pick=None, ret_z=0,z_resolution=1.0,cs_over = None,cscfad = None):
    # pick a variable and do a CFAD for the cell

        if value_bins is None: # set a default if nothing is there
            value_bins = np.linspace(self.lims[var][0], self.lims[var][1], 20)
        else:
            pass
        #print('value bins',value_bins)
        nbins = value_bins.shape[0]
        if above is not None:
            bot_index, top_index = self._get_ab_incides(above=above, below=below)
        #print self.data.keys()
#         data = self._pick_data(self.data[var], pick)
#         if 't' in self.data.dims.keys():
#             data = np.mean(data,axis=0)
#             
#         print np.shape(data),type(data)

        if np.mod(z_resolution, self.dz) != 0:
            print('Need even multiple of vertical resolution: %.1f'%self.dz)
            return
#        print('ln 1592')
        multiple = int(z_resolution/self.dz)
        if 'd' in self.data[self.z_name].dims:
            sz=np.shape(self.data[self.z_name].sel(d=0).values)[0]
            hts = np.squeeze(self.data[self.z_name].sel(d=0).values)
        else:
            sz=np.shape(self.data[self.z_name].values)[0]
            hts=np.squeeze(self.data[self.z_name].values)
        #print( np.shape(sz),sz, multiple)
        looped = np.arange(0, sz, multiple)
        #cfad_out = np.zeros((sz//multiple, nbins-1))
        cfad_out = np.zeros((len(looped), nbins-1))
        #print (np.shape(cfad_out),'cfad shape')
        #print looped
        #print cfad_out.shape, multiple
#        if tspan == None:
        tsi = 0
        tei = self.ntimes-1
#         else:
#             ts=tspan[0]
#             te=tspan[1]
#             tsi = self.get_ind(ts,np.array(self.date))
#             tei = self.get_ind(te,np.array(self.date))
# 
#        print 'cscfad',cscfad
        
        holddat = deepcopy(self.data[var].values)
        if cs_over:
            if cscfad == 'convective':
                #mask = np.where(self.raintype != 2)
               mask= np.where(self.raintype != 2)
               self.data[var].values[mask] = np.nan
    #           print ('in conv'
            elif cscfad == 'stratiform':
               mask = np.where(self.raintype != 1)
               self.data[var].values[mask] = np.nan
    #           print ('in strat',type(self.data[var].data))
            else:
               mask = np.where(self.raintype > 100)
         #      print('entering deep copy')
               holddat2 = deepcopy(self.data[var].values)
               holddat2[mask] = np.nan
               #self.data[var].values[mask] = np.nan
               self.data[var].values = holddat2
            #print('ready to go in loop!')
            # if left blank, check the whole thing
            
            
        for ivl, vl in (enumerate(tqdm(looped[:-1]))):
            #print ivl, vl
#             try:
            #cfad, ed = (np.histogram(mc3e_wrf.data['w'].sel(z=slice(3,15))))            
            #print ts,te,'ts,te'
            v = hts[vl]
            v2 = hts[vl+multiple]
#            print v, v2, 'line1527'
            try:
                dum = self.data[var].sel(z=slice(v,v2)).where(np.isfinite)
            except:
#                v = self.get_ind(vl,self.data[self.z_name].data)
#                v2 = self.get_ind(vl+multiple,self.data[self.z_name].data])
#                print v,v2
                dum = self.data[var].sel(z=slice(vl,vl+multiple)).where(np.isfinite)
            #print np.max(dum)
#            dum2 = np.ma.masked_less(dum,-900.0)
#             dum2 = np.where(np.isfinite(dum))
#             print(type(dum[dum2]))
            #print dum[dum2]
            lev_hist, edges = np.histogram(dum, bins=value_bins, density=True) 
            '''
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.hist(dum.values.flatten(), bins=value_bins, density=True)
            if self.mphys.startswith('obs'):
                plt.savefig('hist1x_'+var+'_'+str(ivl).zfill(2)+'.png')
            else:
                plt.savefig('hist2x_'+var+'_'+str(ivl).zfill(2)+'.png')
            plt.close()
            '''
#           except:
#            lev_hist, edges = np.histogram(data[vl:vl+multiple].ravel(), bins=value_bins, density=True) 
            #print lev_hist, edges
                    # this keeps it general, can make more elaborate calls in other functions
            #print(np.shape(lev_hist))
            
            lev_hist = 100.0*lev_hist/np.sum(lev_hist)
            if np.max(lev_hist) > 0:
                #print(looped[ivl])
                #print('Hello')
                #print(max(lev_hist))
                cfad_out[ivl, :] = lev_hist
            #else:
            #    print('Negative!')
            #    print(looped[ivl])
            #    print(max(lev_hist))
        #print np.shape(cfad_out)
#        if cscfad == 'convective' or cscfad == 'stratiform':
#            print 'setting data back'
        self.data[var].values = holddat

        if ret_z == 1:
        
            return cfad_out,value_bins,hts[looped]
        else:
            return cfad_out

#############################################################################################################

    def cfad_multiplot(self, varlist=None, z_resolution=1.0, cs_over=False, zmax=None, **kwargs):

        st = self.date[0].strftime('%Y%m%d_%H%M%S')
        en = self.date[-1].strftime('%Y%m%d_%H%M%S')
        if st.startswith(en): dtlab = st[0:4]+'-'+st[4:6]+'-'+st[6:8]+' '+st[9:11]+':'+st[11:13]+' UTC'
        else:
            if st[0:8].startswith(en[0:8]): dtlab = st[0:4]+'-'+st[4:6]+'-'+st[6:8]+' '+st[9:11]+':'+st[11:13]+'-'+en[9:11]+':'+en[11:13]+' UTC'
            else: dtlab = st[0:4]+'-'+st[4:6]+'-'+st[6:8]+' '+st[9:11]+':'+st[11:13]+' - '+en[0:4]+'-'+en[4:6]+'-'+en[6:8]+' '+en[9:11]+':'+en[11:13]+' UTC'

        if varlist is not None:
            good_vars = varlist
        else:
            good_vars = self.valid_vars()
            
        nvars = len(good_vars)
        if 'scores' in good_vars:
            if hasattr(self, 'scores'):
                nvars += 1
        if nvars <= 3:
            ncols = 1
            nrows = deepcopy(nvars)
            figx = 5
            figy = 3*nrows
        else:
            nrows = 2
            ncols = int(np.ceil(nvars/2))

        fig, ax = plt.subplots(nrows,ncols,figsize=(14,8),gridspec_kw={'wspace': 0.1, 'hspace': 0.2, 'top': 1., 'bottom': 0., 'left': 0., 'right': 1.})
        if not isinstance(ax, np.ndarray) or not isinstance(ax, list): 
            ax = np.array([ax], **kwargs)
        axf = ax.flatten()
 
        for i in range(len(good_vars),ncols*nrows):
            fig.delaxes(axf[i])
        
        for i, var in enumerate((good_vars)):
            lspanels = [ncols*n for n in range(0,nrows)]
            ylabbool = True if i in lspanels else False
            if 'HID' in good_vars:
                cbpanels = [ncols*n+2 for n in range(0,nrows)]
                cbbool = True if i in cbpanels else False
            else:
                cbbool = 0
            if not var.startswith('HID'):
                cfad_ma, hts, pc, figout, axout = self.cfad_plot(var,ax=axf[i],ylab=ylabbool,cbar=cbbool,bins=self.cfbins[var],levels=1,zmax=zmax,z_resolution=z_resolution,cs_over=cs_over)
            else:
                self.plot_hid_cdf(ax=axf[i],ylab=ylabbool,zmax=zmax,z_resolution=z_resolution,cs_over=cs_over)
 
        if 'HID' not in good_vars:
            cbar_ax_dims = [1.02,0.0,0.03,1.0]
            cbar_ax = fig.add_axes(cbar_ax_dims)
            cbt = plt.colorbar(pc,cax=cbar_ax)
            cbt.ax.tick_params(labelsize=16)
            cbt.set_label('Frequency (%)', fontsize=16, rotation=270, labelpad=20)
            cbt.set_ticks(self.cfad_levs)
              
        axf[0].text(0, 1, '{e} {r}'.format(e=self.exper,r=self.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=20, color='k', zorder=10, weight='bold', transform=axf[0].transAxes) # (a) Top-left
        axf[ncols-1].text(1, 1, dtlab, horizontalalignment='right', verticalalignment='bottom', size=20, color='k', zorder=10, weight='bold', transform=axf[ncols-1].transAxes) # (a) Top-left
 
        return fig, ax

#############################################################################################################

    def cfad_plot(self, var, cfad=None, hts=None, nbins=20, ax=None, maxval=10.0, above=None, below=15.0, bins=None, log=False, diff=False, pick=None, z_resolution=1.0, levels=None, cont=False, cs_over=False, cscfad=False, cbar=1, xlab=1, ylab=1, zmax=8,**kwargs):

        from matplotlib.colors import from_levels_and_colors

        if bins is None:
            bins = np.linspace(self.lims[var][0], self.lims[var][1], nbins)
        else:
            pass

        multiple = int(z_resolution/self.dz)

        if np.nanstd(self.data[var].values) < abs(bins[0]-bins[-1]):
            if self.names_uc[var].startswith('RHO'):
                bins = np.arange(0.95,bins[-1],0.001)
        
        if not diff:
            cfad,value_bins,hts = self.cfad(var, value_bins=bins, above=above, below=below, pick=pick, z_resolution=z_resolution, cs_over=cs_over, cscfad=cscfad,ret_z=1)

        if above is not None:
            bot_index, top_index = self._get_ab_incides(above=above, below=below)

        if ax is None:
            fig = plt.figure(figsize=(10,8))
            ax = fig.add_subplot(111)
        else:
            fig = ax.get_figure()

        cfad_ma = np.ma.masked_where(cfad==0,cfad)
        if log:
            norm = colors.LogNorm(vmin=1e-5, vmax=1e2)
        elif diff:
            cmap = plt.set_cmap('bwr')
            vmin,vmax = -np.nanmax(cfad_ma),np.nanmax(cfad_ma)
            norm = colors.Normalize(vmin=vmin,vmax=vmax) 
        else:
            cmap, norm = from_levels_and_colors(self.cfad_levs,self.cfad_cols)
        
        if cont is True:
            pc = ax.contourf(bins[0:-1],hts,(cfad_ma),levels=self.cfad_levs,color=self.cfad_cols,cmap=cmap,norm=norm,extend='both')
        else:
            if levels is not None:
                pc = ax.pcolormesh(bins[0:-1], hts, cfad_ma, norm=norm, cmap=cmap)
            else:
                pc = ax.pcolormesh(bins, hts, cfad_ma, vmin=0, vmax=maxval, norm=norm,cmap=cmap, **kwargs)
      
        if xlab:
            ax.set_xlabel('%s %s' %(self.names_uc[var], self.units[var]),fontsize=16)
            ax.tick_params(axis='x',labelsize=16)
        else:
            ax.tick_params(axis='x',labelsize=0,left=False)
        if cbar:
            cbar = self.mycbar(fig,ax,pc,'Frequency (%)')
            cbar.set_ticks(self.cfad_levs)
        if ylab: 
            ax.set_ylabel('Height (km MSL)',fontsize=16)
            ax.tick_params(axis='y',labelsize=16)
        else: 
            ax.tick_params(axis='y',labelsize=0,left=False)
        
        ax.set_xlim(np.min(bins),np.max(bins))
        ax.set_ylim(0,zmax)
        ax.grid(color='grey', linestyle='-', linewidth=1)

        return cfad_ma, hts, pc, fig, ax

#############################################################################################################
    
    def hist2d(self, varx=None, vary=None, binsx=None, binsy=None, above=None, below=None, pick=None,xthr = -900.0,ythr=-900.0):
        # This will just look at the whole volume
#        if above is None:
            
 #       bot_index, top_index = self._get_ab_incides(above=above, below=below)
        if above is None:
            above = 0
        if below is None:
            below = self.nhgts#[0]
        
        datax = self.data[varx].sel(z=slice(above,below)).values
        datay = self.data[vary].sel(z=slice(above,below)).values
                    
        dum = np.less(datax,xthr)
        dumy = np.less(datay,ythr)
        datax[dum]=np.nan
        datay[dumy]=np.nan

        hist, edges = gentools.hist2d(varx=datax, vary=datay, 
                binsx=binsx, binsy=binsy)

#         hist, edges = gentools.hist2d(varx=datax[bot_index:top_index], vary=datay[bot_index:top_index], 
#                 binsx=binsx, binsy=binsy)

        return hist, edges

#############################################################################################################
#############################################################################################################

    def plot_2dhist(self, hist,edge,ax=None,cbon = True):
        if ax is None:
            fig, ax = plt.subplots()
        else:
        # ax has been passed in, do nothing to ax, but need to get the parent fig
            fig = ax.get_figure()
    
        cb = ax.contourf(edge[0][:-1],edge[1][:-1],hist.T,norm=colors.Normalize(vmin=0, vmax=max(hist)),levels=np.arange(0.01,max(hist),0.01))
        if cbon == True:
            #print ' making colorbar'
            col = plt.colorbar(cb,ax=ax)
            col.ax.tick_params(labelsize=24)

        ax.tick_params(axis='both',labelsize=24)
        
        # This will just look at the whole volume
#        if above is None:
        return fig, ax


#############################################################################################################
    def percentile(self,wup =False,wdown=False,use_frzht=False,frzhgt=None):

        if wup is not False:
            wdat = deepcopy(self.data[self.w_name].values)
            wdat[wdat<=0] = np.nan
            t1=99
            t2=90
            t3=50
        elif wdown is not False:
            wdat = deepcopy(self.data[self.w_name].values)
            wdat[wdat>=0] = np.nan
            t1=1
            t2=10
            t3=50
        else:
            wdat =deepcopy(self.data[self.w_name].values)

            t1 = 99
            t2 = 90
            t3 = 50
    #    print np.ma.min(np.ma.compressed(winds))
        # Time to make the percentile lists
        if 'd' in self.data[self.z_name].dims:
            height = self.data[self.z_name].sel(d=0).values
        else:
            height =self.data[self.z_name].values
        percentile_50 = np.zeros([np.size(height)])
        percentile_90 = np.zeros([np.size(height)])
        percentile_99 = np.zeros([np.size(height)])
        mean = np.zeros([np.size(height)])
    #    print np.shape(winds)
    #    print np.shape(height)
        for i,ht in enumerate(height):
            #summed_z = np.ma.sum(winds[:,i,:,:])
            try:
                percentile_90[i]=np.nanpercentile(np.ravel(wdat[:,i,...]),t2) #use w/ masked arrays

                #percentile_90[i]=np.percentile(winds[:,i,:,:],90)
                #percentile_90[i]=np.percentile(winds[i],99)
            except IndexError:
                percentile_90[i]=np.nan
            try:
                percentile_99[i]=np.nanpercentile(np.ravel(wdat[:,i,...]),t1) #use w/ masked arrays
                #percentile_99[i]=np.percentile(winds[:,i,:,:],99)
                #percentile_99[i]=np.percentile(winds[i],99)
            except IndexError:
                percentile_99[i]=np.nan
            try:
                percentile_50[i]=np.nanpercentile(np.ravel(wdat[:,i,...]),t3) #use w/ masked arrays
                #percentile_50[i]=np.percentile(winds[:,i,:,:],50)
                #percentile_50[i]=np.percentile(winds[i],50)
            except IndexError:
                percentile_50[i]=np.nan
            #mean[i]=np.ma.mean(winds[:,i,:,:]) #use w/ masked arrays
            #mean[i]=np.mean(winds[:,i,:,:])

        ## This part is only necessary for the mixing ratio stuff ##
        if use_frzht==True:

            for i,hgt in enumerate(height):
    #            print np.float(hgt),frzhgt+np.float(1)
                if np.float(hgt) > frzhgt+np.float(1):
                    percentile_50[i] = 0
                    percentile_90[i] = 0
                    percentile_99[i] = 0
    #                print percentile_50[i]
    #    print percentile_50

    #         plt.xlim([-1,1])
        return percentile_99,percentile_90,percentile_50,height


#               ADD some percentile plots   
#############################################################################################################

    def percentileplots(self,ptype='W',punit='m/s',use_frzht=False,frzhgt=None):
    #    print np.ma.min(np.ma.compressed(winds))
        # Time to make the percentile lists
        height = self.data[self.z_name].data
        percentile_50 = np.zeros([np.size(height)])
        percentile_90 = np.zeros([np.size(height)])
        percentile_99 = np.zeros([np.size(height)])
        mean = np.zeros([np.size(height)])
    #    print np.shape(winds)
    #    print np.shape(height)
        for i,ht in enumerate(height):
            #summed_z = np.ma.sum(winds[:,i,:,:])
            try: 
                percentile_90[i]=np.nanpercentile(np.ravel(self.data[self.w_name].values[i,...]),90) #use w/ masked arrays

                #percentile_90[i]=np.percentile(winds[:,i,:,:],90)
                #percentile_90[i]=np.percentile(winds[i],99)
            except IndexError:
                percentile_90[i]=np.nan
            try:
                percentile_99[i]=np.nanpercentile(np.ravel(self.data[self.w_name].values[i,...]),99) #use w/ masked arrays
                #percentile_99[i]=np.percentile(winds[:,i,:,:],99)
                #percentile_99[i]=np.percentile(winds[i],99)
            except IndexError:
                percentile_99[i]=np.nan
            try:
                percentile_50[i]=np.nanpercentile(np.ravel(self.data[self.w_name].values[i,...]),50) #use w/ masked arrays
                #percentile_50[i]=np.percentile(winds[:,i,:,:],50)
                #percentile_50[i]=np.percentile(winds[i],50)
            except IndexError:
                percentile_50[i]=np.nan
            #mean[i]=np.ma.mean(winds[:,i,:,:]) #use w/ masked arrays
            #mean[i]=np.mean(winds[:,i,:,:])

        ## This part is only necessary for the mixing ratio stuff ##
        if use_frzht==True:

            for i,hgt in enumerate(height):
    #            print np.float(hgt),frzhgt+np.float(1)
                if np.float(hgt) > frzhgt+np.float(1):
                    percentile_50[i] = 0
                    percentile_90[i] = 0
                    percentile_99[i] = 0
    #                print percentile_50[i]
    #    print percentile_50 

    #         Now doing the plotting
        p2 = plt.plot(percentile_50,height,color='red',linewidth=2,label='50th')
        p1 = plt.plot(percentile_90,height,color='orange',linewidth=2,label='90th')
        p4 = plt.plot(percentile_99,height,color='black',linewidth=2,label='99th')
    #         p5 = plt.plot(mean,z,color='black',linewidth=2,label='mean')
    #         p5 = plt.plot(mean,z,color='blue',linewidth=2,label='mean')
        lns = p2+p1+p4
    #   lns = p4
        labs = [l.get_label() for l in lns]
        lgd = plt.legend(lns,labs,bbox_to_anchor=(1.6,0.85),prop={'size':14})
    #         plt.title('MC3E 1 May 2011 90th Percentile, 99th Percentile, \nand Mean of Downward Winds',size=16,y=1.08,x=0.78)
        plt.title('{c} {dd} 50th, 90th, and 99th Percentile of {p}'.format(c=self.exper,dd=self.date,p=ptype),size=16,y=1.08,x=0.78)
        #plt.xlabel('Rainrate (mm hr$^-$$^1$)',size=14)
        plt.xlabel('{p} ({u})'.format(p=ptype,u=punit))
        plt.ylabel('Height (km)',size=14)
        #plt.ylim([0,6]) # For qr plots
    #         plt.xlim([-1,1])
        return plt
#############################################################################################################

    def vertical_hid_volume(self, hid_nums, z_resolution=1.0, above=None, below=None, pick=None, cs_over = False, cscfad = None):
        # This function only returns the height and volume arrays, no plotting
        if above is None:
            above = 0
        if below is None:
            #print self.nhgts[0]
            below = self.nhgts#[0]

        #bot_index, top_index = self._get_ab_incides(above=above, below=below)
        #data = self._pick_data(self.data[self.hid_name].data, pick)
        #print above,below,np.shape(self.data[self.hid_name].data)
        holddat = deepcopy(self.data[self.hid_name].values)
        if cs_over:

            if cscfad == 'convective':
                mask= np.where(self.raintype != 2)
                self.data[self.hid_name].values[mask] = -1
    #            print 'in conv'
            elif cscfad == 'stratiform':
                mask= np.where(self.raintype != 1)
                self.data[self.hid_name].values[mask] = -1
    #            print 'in strat'
            else:
                mask = np.where(self.raintype < 100)
            
#        print(below,above)
        data = self.data[self.hid_name].sel(z=slice(above,below)).values
        
#        print('datashape',data.shape,above,below)
        msk = (np.less(self.data[self.dz_name].sel(z=slice(above,below)).values, -900))
        data[msk] = -1
        #print 'vhv dat',np.shape(data)
        # Not sure if I need this here.....
        if np.mod(z_resolution, self.dz) != 0:
            print ('Need even multiple of vertical resolution: %.1f'%self.dz)
            return

        multiple = int(z_resolution/self.dz)
        vol = np.zeros(int(self.data[self.z_name].values.shape[0]/multiple))
        hts = np.zeros(int(self.data[self.z_name].values.shape[0]/multiple))
        #print np.shape(vol)
        #print self.data[self.z_name].data.shape[1]
        if 'd' in self.data[self.z_name].dims:
#            print('in d version of looped')
            looped = np.arange(0, int(self.data[self.z_name].values.shape[1]), multiple)
#            print('looped',looped)
            vol = np.zeros(int(self.data[self.z_name].values.shape[1]/multiple)+1)
            hts = np.zeros(int(self.data[self.z_name].values.shape[1]/multiple)+1)
        else:
            looped = np.arange(0, int(self.data[self.z_name].values.shape[0]), multiple)
            vol = np.zeros(int(self.data[self.z_name].values.shape[0]/multiple))
            hts = np.zeros(int(self.data[self.z_name].values.shape[0]/multiple))
        for vi,vl in enumerate(looped):
            lev_hid = data[:,vl:vl+multiple,...] # go to vl+multiple cuz not inclusive
            #print 'lev_hid',np.shape(lev_hid)
#            print hid_nums, np.shape(lev_hid)
            where_this_hid = np.where(lev_hid == hid_nums)
#            print np.shape(where_this_hid)
            this_hid_vol = where_this_hid[0].shape[0]
            vol[vi] += this_hid_vol
            #print self.data[self.z_name].data[0][vl+multiple]
        if 'd' in self.data[self.z_name].dims:
            hts = self.data[self.z_name].sel(d=0).values[looped]
        else:
            hts = self.data[self.z_name].values[looped]

        if cscfad == 'convective' or cscfad == 'stratiform':
            self.data[self.hid_name].values = holddat

        #print np.shape(vol), np.max(vol)
        #print hts
#         print self.data[self.z_name].data[0]
#         print self.data[self.z_name].data[0][::looped]
#         hts = self.data[self.z_name].data[::multiple]
#        print 'hts in vert hid vol', np.shape(hts)
        return hts, vol


#############################################################################################################
#############################################################################################################

    def hid_vertical_fraction(self, hid_nums, z_resolution=1.0, above=None, below=None, pick=None,cs_over = False,cscfad = None):

        if np.mod(z_resolution, self.dz) != 0:
            print ('Need even multiple of vertical resolution: %.1f'%self.dz)
            return None
        else:
            #hts = self.data[self.z_name].data[0][::int(np.round(z_resolution/self.dz))]
            if 'd' in self.data[self.z_name].dims:
                hts2 = np.squeeze(self.data[self.z_name].sel(d=0,z=slice(above,below)).data)
            else:
                hts2 = np.squeeze(self.data[self.z_name].sel(z=slice(above,below)).data)
            #hts2 = (self.data[self.z_name].sel(t=slice(0,-1),z=slice(above,below))).data
            #print np.shape(hts2)
            hts = hts2[::int(np.round(z_resolution/self.dz))]
            # This gives the percent of the storm that is taken up by the combo of HID values
        hid_nums = np.asarray(hid_nums)

        hidcdf,hts = self.hid_cdf(z_resolution=z_resolution, cs_over = cs_over, cscfad = cscfad)
        hvf = np.zeros(hidcdf.shape[1])
#        print 'hvf in hvf', np.shape(hvf)
    # now loop thru each hid_num
        for hn in hid_nums:
            if hn == 1:
                hvf += hidcdf[hn-1, :]
            else:
                hvf += hidcdf[hn-1, :] - hidcdf[hn-2, :]

        return hts, hvf


#############################################################################################################

    def hid_cdf(self, z_resolution=1.0, pick=None,cs_over=False, cscfad = None):
        # vertical HID_cdf with bar plots I think
        if np.mod(z_resolution, self.dz) != 0:
            print ('Need even multiple of vertical resolution: %.1f'%self.dz)
            return None

        # loop thru the species and just call the vertical hid volume
        all_vols = []
        for sp in range(len(self.species)):
            
            hts, dat = self.vertical_hid_volume([sp+1], z_resolution=z_resolution, pick=pick,cs_over=cs_over, cscfad = cscfad)
            all_vols.append(dat) # need the +1
            

        all_vols = np.array(all_vols)
        all_cdf = np.zeros_like(all_vols)
#9        print np.shape(all_vols)
        # shape is 10,16, which is nspecies x nheights
        # need to do cdf on each level
        for iz in range(all_vols.shape[1]):
            # loop thru the vertical
            level_cum_vol = np.cumsum(all_vols[:, iz])
#            print level_cum_vol
            all_cdf[:, iz] = 100.0*level_cum_vol/level_cum_vol[-1]

        return all_cdf,hts


#############################################################################################################

    def HID_barplot_colorbar(self, figure, location = [0.9, 0.1, 0.03, 0.8], orientation='vertical', names='shortnames', lblsz=16):

        scalarMap = plt.cm.ScalarMappable(norm=self.normhid,cmap=self.hid_cmap)
        axcb = figure.add_axes(location)
        cb = mpl.colorbar.ColorbarBase(axcb, cmap=self.hid_cmap, norm=self.normhid, boundaries=self.boundshid, orientation = orientation)
        cb.set_ticks(np.arange(len(self.species))+0.5)
        
        if names.startswith('long'): 
            labs = np.array(self.species_long)
            cb.set_ticklabels(labs)
            cb.ax.tick_params(labelsize=lblsz)
        else: 
            labs = np.array(self.species)
            cb.set_ticklabels(labs)
            cb.ax.tick_params(labelsize=lblsz)
         
        return cb

    def plot_hid_cdf(self, ax=None, xlab=1, ylab=1, zmax=None, cbar=1, data=None, z_resolution=1.0, pick=None, cs_over = False, cscfad = None):
        
        ts=np.array(self.date)[0]
        if data is not None:
            pass
        else:
            pass
            data,hgt = self.hid_cdf(z_resolution=z_resolution, pick=pick,cs_over = False,cscfad = cscfad)
        if ax is None:
            fig = plt.figure(figsize=(10,8))
            ax = fig.add_subplot(111)
        else:
            fig = ax.get_figure()

        multiple = int(z_resolution/self.dz)
        
        for i in tqdm(np.arange(0,len(hgt),multiple)):
            ax.barh(hgt[i],data[0,i],left=0.,align='center',height=z_resolution,color=self.hid_colors[0],edgecolor='k')
            for spec in range(1, len(self.species)):
                ax.barh(hgt[i],data[spec,i],left=data[spec-1,i],align='center',height=z_resolution,color=self.hid_colors[spec],edgecolor='k')

        ax.set_xlim(0,100)
        if xlab == 1:
            ax.set_xlabel('Cumulative Frequency (%)',fontsize=16)
            ax.tick_params(axis='x',which='major',labelsize=16)
        else:
            ax.tick_params(axis='x',which='major',labelsize=0)
            ax.set_xticks([])
            ax.set_xticklabels([])
            
        ax.set_ylim(0,zmax)
        if ylab == 1:
            ax.set_ylabel('Height (km MSL)',fontsize=16)
            ax.tick_params(axis='y',which='major',labelsize=16)
        else:
            ax.tick_params(axis='y',which='major',labelsize=0)
            ax.set_yticks([])
            ax.set_yticklabels([])
        
        if cbar == 1:
            lur,bur,wur,hur = ax.get_position().bounds
            cbar_ax_dims = [lur+wur+0.015,bur,0.03,hur]
            self.HID_barplot_colorbar(fig,cbar_ax_dims)
        
        if cbar == 2:
            lur,bur,wur,hur = ax.get_position().bounds
            cbar_ax_dims = [lur,bur-0.125,wur,0.03]
            self.HID_barplot_colorbar(fig,cbar_ax_dims,orientation='horizontal',names='longnames')

        return fig, ax 

#############################################################################################################

    def vertical_profile(self, var, rep_func=np.average, above=None, below=None, pick=None):

        bot_index, top_index = self._get_ab_incides(above=above, below=below)
        #print bot_index,top_index
        vp = []
        heights = []

        data = self._pick_data(self.data[var], pick)

        for iz in range(bot_index, top_index):
            if isinstance(data, np.ma.core.MaskedArray):
                this_data = np.ma.compressed(data[iz,:,:])
            else:
                this_data = data[iz,:,:]
            if len(this_data):
                value = rep_func(this_data)
            else:
                value = np.nan
        vp.append(value)
        heights.append(self.data[self.z_name][iz])

        vp = np.ma.masked_invalid(np.array(vp))
        vp = np.ma.masked_where(vp < -100, vp)
        heights = np.array(heights)

        return heights, vp



#############################################################################################################
#############################################################################################################
#############################################################################################################

#############################################################################################################



    def find_updraft_downdraft(self, updraft_thresh=1.0, downdraft_thresh=-1.0):
        # This will make a mask for updraft and downdraft regions based on some thresholds
        # First, will make sure that the updraft thresh is > 0 and the downdraft is < 0
        if (updraft_thresh >= 0) and (downdraft_thresh < 0):
             self.updraft = self.data[self.w_name] >= updraft_thresh
             self.downdraft = (self.data[self.w_name] <= downdraft_thresh) & (self.data[self.w_name] >= -99.0)

        else:
            print ('Make sure your thresholds are valid')
            return



#############################################################################################################

    def updraft_percentile(self, perc):
        return np.percentile(self.data[self.w_name][self.updraft], perc)

#############################################################################################################

    def downdraft_percentile(self, perc):
        #print self.data[self.w_name][self.downdraft]
        try:
            return np.percentile(self.data[self.w_name][self.downdraft], perc)
        except IndexError:
            return None

#############################################################################################################


    def updraft_vol(self, up_val, above=None, below=None):

        bot_index, top_index = self._get_ab_incides(above=above, below=below)
    
        # if left blank, check the whole thing
        uv = 0. # initialized the updraft volume
        for vl in range(bot_index, top_index): 
            lev_w = self.data[self.w_name][vl]
            in_vol = np.where(lev_w >= up_val)      
            uv += in_vol[0].shape[0]*self.dx*self.dy*self.dz

        return uv

#############################################################################################################

    def downdraft_vol(self, down_val, above=None, below=None):

        bot_index, top_index = self._get_ab_incides(above=above, below=below)

        # if left blank, check the whole thing
        dv = 0. # initialized the updraft volume
        for vl in range(bot_index, top_index): 
            lev_w = self.data[self.w_name][vl]
            in_vol = np.where( (lev_w <= down_val) & (lev_w >= -99.0) )
            dv += in_vol[0].shape[0]*self.dx*self.dy*self.dz

        return dv

#############################################################################################################

    def w_profile(self, rep_func=np.average, masked=True):
        # returns a vertical profile of vertical velocity with a given sign (up, down, all)
        # using a specified representative function. Default is average, but could do something like a median

        # first, set up the out, which will have a spot for each z and one for up/down draft
        vert_motion = deepcopy(self.data[self.w_name])
        out = {'up': np.zeros(self.data[self.z_name].shape), 'down': np.zeros(self.data[self.z_name].shape)}

        for il in range(self.data[self.z_name].shape[0]):
            good_up = np.where(self.updraft[il])[0]
            good_dn = np.where(self.downdraft[il])[0]
            if good_up.shape[0] > 0:
                out['up'][il] = rep_func(vert_motion[il][self.updraft[il]])
            else:
                out['up'][il] = np.nan

            if good_dn.shape[0] > 0:
                out['down'][il] = rep_func(vert_motion[il][self.downdraft[il]])
            else:
                out['down'][il] = np.nan


        if masked:
            out['up'] = np.ma.masked_invalid(out['up'])
            out['down'] = np.ma.masked_invalid(out['down'])


        return out

#############################################################################################################


    def plot_w_profile(self, rep_func=np.average, masked=True, ax=None):
        # this will plot the vertical profile
        wprof = self.w_profile(rep_func=rep_func, masked=masked)

        if ax is None:
            fig, ax = plt.subplots()
        else:
        # ax has been passed in, do nothing to ax, but need to get the parent fig
            fig = ax.get_figure()

        u = ax.plot(wprof['up'], self.data[self.z_name], color='red', linewidth=2, label='updrafts')
        d = ax.plot(wprof['down'], self.data[self.z_name], color='blue', linewidth=2, label='downdrafts')

        max_abs = abs(np.array(ax.get_xlim())).max()

        ax.set_xlim(-1*max_abs, max_abs)
        #ax.set_xlim()
        ax.axvline(x=0, color='black')
        ax.set_xlabel('Vertical motion (m/s)')
        ax.set_ylabel('Altitude (km MSL)')
        ax.grid(True)
        ax.legend(loc='best')
        ax.set_title('%s %s Vertical motion profile' % (self.print_date(), self.band+'-band'))

        return fig, ax


#############################################################################################################

    def updraft_width_profile(self, thresh=5.0, temps=np.arange(20,-60,-5),thresh_dz=False):
        import scipy.interpolate as sint
        # this gets the width of the updraft as a function of temperature
        if 'd' in self.data[self.z_name].dims:
            uw = np.zeros(self.data[self.z_name].sel(d=0).shape)
        else:        
            uw = np.zeros(self.data[self.z_name].shape)
        #temps = self.data[self.z_name].data[0,:]
        # basically just loop thru the Z and get the associated temperature and area

        #data = self.data[self.w_name].data
        #data_dz = self.data[self.dz_name].data
        data = self.data[self.w_name].values
        data_dz = self.data[self.dz_name].values
        if thresh_dz == True:
            data[data_dz < self.z_thresh] = np.nan
#        print np.shape(data),'ln2075'
        for iz, z in enumerate(self.data[self.z_name].data):
            #values_above = np.where(data[iz,...] >= thresh)[0]
            values_above = np.where(data[:,iz,:,:] >= thresh)[0]
            num_above = len(values_above)
            uw[iz] = num_above*self.dx*self.dy/self.ntimes
            if self.data[self.x_name].units == "[deg]":
#                print 'changing units'
                uw[iz]=uw[iz]*110.*110.
        #print np.shape(uw),max(uw),self.data[self.x_name].units
        #print np.shape(uw)
        #print np.shape(self.T[0,:,0,0])
        # now inerpolate this to the temps listed
        self.T = xr.DataArray(data=self.T,dims=['d','z','y','x'])

        if 'd' in self.T.dims:
            print('shapes in updraft width',np.shape(uw),np.shape(self.T.sel(x=0,y=0)))
            f_temp_u = sint.interp1d(self.T.sel(d=0,x=0,y=0), uw, bounds_error=False)
        else:
            f_temp_u = sint.interp1d(self.T[:,0,0], uw, bounds_error=False)
        uwp_interp = f_temp_u(temps)
        #return temps, uw
        return temps,uwp_interp

#############################################################################################################

    def def_convstrat(self):
        #print self.conv_types,self.strat_types
        for k in self.conv_types:
            
            self.raintype[self.raintype == self.rtypes[k]] =2
        for k in self.strat_types:
            self.raintype[self.raintype == self.rtypes[k]] =1

        for k in self.mixed_types:
            self.raintype[self.raintype == self.rtypes[k]] =3

#############################################################################################################

    def calc_timeseries_stats(self,var,ht_lev = 3,thresh=-99.,cs_flag=False,make_zeros=False,areas=False):
        #First calculate the domain averaged rain rates at a given level.
        #if self.rr_name is not None:
        #    rr_timeseries_uncond = rr.mean(dim=['x','y','z'],skipna=True)
        data = deepcopy(self.data[var].sel(z=slice(int(ht_lev),int(ht_lev)+1)))
        whbad= np.where(data.values<thresh)

        if make_zeros==True:
            print('filling with zeros')
            data.values[whbad] = 0.0
            data =data.fillna(0.0)
        else:
            data.values[whbad] = np.nan
        
        #print('dat max:',data.values.max())

        if cs_flag is not False:
            datstrat = data.where(self.data[self.cs_name].sel(z=slice(1,2))==1)
            datconv = data.where(self.data[self.cs_name].sel(z=slice(1,2))==2)
            datall = data.where(self.data[self.cs_name].sel(z=slice(1,2))>0)       
            if areas==False:
                datstrat_ts= datstrat.mean(dim=['z','y','x'],skipna=True)
                datconv_ts= datconv.mean(dim=['z','y','x'],skipna=True)
                datall_ts= datall.mean(dim=['z','y','x'],skipna=True)
            else:
                cssum =(self.data[self.cs_name].max(dim='z',skipna=True))
                datstrat_ts = cssum.where(cssum==1).count(dim=['y','x'])
                datconv_ts = cssum.where(cssum==2).count(dim=['y','x'])
                datall_ts = cssum.where(cssum>0).count(dim=['y','x'])
            return datstrat_ts,datconv_ts,datall_ts
        else:

            datall = data
            if areas==True:
                datall_ts = datall.mean(dim=['z','y','x'],skipna=True)
            else:
                datall_ts = datall.count(dim=['z','y','x'])        
            return datall


#############################################################################################################

    def hid_q_compare(self):
        hidmassq = {}
        for s in np.unique(self.hid):
            whtp = np.where(self.hid == s)
            count = len(whtp)
            vol = count*self.dx*self.dy
            mass = {}
            for q in eval(self.mixr):
                #print q
                mass[q]=np.sum(self.data[q].data[whtp])*vol
            hidmassq[s] = mass
        self.hidmassq = hidmassq
        
    def score_q_corr(self):
        corr = np.zeros([len(self.scores[:,0,0,0]),len(eval(self.mixr))])
        #print np.shape(corr)
        for q,i in enumerate((self.scores[:,0,0,0])):
            #print q
            for h,j in enumerate(eval(self.mixr)):
                #print j
                #print q,i,h,j
                whgd = np.where(np.logical_and(np.isfinite(np.ravel(self.scores[q,...])),np.isfinite(np.ravel(self.data[j].data))))
                corr[q,h] = np.corrcoef(np.ravel(self.scores[q,...])[whgd],np.ravel(self.data[j].data)[whgd])[0,1]
                #print 'corr!',corr[q,h]
                #corr[q,h] = np.corrcoef(self.scores[i,...],self.data[j].data)
        self.scoreqcorr = corr

#############################################################################################################
    def get_latlon_fromxy(self):
    
        proj = Proj(init='epsg:3857')

        lon_0 = self.lon_0
        lat_0 = self.lat_0
    
####Bea little careful here. you need to make sure you are in the correct hemispheres
        print('lat / lon:',lat_0,lon_0)
        p = Proj('+proj=lcc +lon_0={n} +lon_1 = {n} +lat_1={t} +lat_2={t} +lat_0={t}'.format(t=self.lat_0,n=self.lon_0)) 
#         if self.lat_0 > 0:
#             if self.lon_0 > 0:
#                 p = Proj('+proj=lcc +a=6370000.0m +lon_0={n}e +lon_1 = {n}e +lat_1={t}n +lat_2=60n +lat_0={t}n'.format(t=self.lat_0,n=self.lon_0)) 
#             else:
#                 print('trying lat lon')
#                 lon_0 = abs(self.lon_0)
#                 p = Proj('+proj=lcc +a=6370000.0m +lon_0={n}w +lon_1 = {n}w +lat_1={t}n +lat_2=60n +lat_0={t}n'.format(t=self.lat_0,n=lon_0)) 
#         else:
#             if self.lon_0 > 0:            
#                 p = Proj('+proj=lcc +a=6370000.0m +lon_0={n}e +lon_1 = {n}e +lat_1={t}s +lat_2=60s +lat_0={t}s'.format(t=self.lat_0,n=self.lon_0)) 
#             else:
#                 p = Proj('+proj=lcc +a=6370000.0m +lon_0={n}w +lon_1 = {n}w +lat_1={t}s +lat_2=60s +lat_0={t}s'.format(t=self.lat_0,n=self.lon_0)) 
        xx, yy = np.meshgrid(self.data[self.x_name], self.data[self.y_name])
        lons, lats = p(xx*1000.,yy*1000.,inverse=True)
    
        self.data['lat'] = ((self.x_name,self.y_name),lats)
        self.data['lon'] = ((self.x_name,self.y_name),lons)
        self.lat_name='lat'
        self.lon_name='lon'
    
    #############################################################################################################
    def calc_cs_shy(self,cs_z=2.0):

        print ('Unfortunately have to run the convective stratiform per timestep. Might take a minute....{n}'.format(n=self.data.dims['d']))
        rntypetot = []
#        print(cs_z,'cs_z in 2424')
        for q in tqdm(range(self.data.dims['d'])):
            if self.lat_name in self.data.keys():
#                         lat = self.data[self.lat_name].sel(d=q).values
#                         lon = self.data[self.lon_name].sel(d=q).values
                if 'd' in self.data[self.lat_name].dims:
                    lat = self.data[self.lat_name].sel(d=q).values
                    lon = self.data[self.lon_name].sel(d=q).values
                else:
                    lat = self.data[self.lat_name].values
                    lon = self.data[self.lon_name].values
                if 'd' in self.data[self.z_name].dims:
                    zlev = np.where(self.data[self.z_name].sel(d=q).values ==cs_z)[0]
                    nlevs = np.shape(self.data[self.z_name].sel(d=q).values)[0]
                else:
                    zlev = np.where(self.data[self.z_name].values ==cs_z)[0]
                    nlevs = np.shape(self.data[self.z_name].values)[0]
            else:
                self.get_latlon_fromxy()
                lat = self.data[self.lat_name].values
                lon = self.data[self.lon_name].values
                zlev = np.where(self.data[self.z_name].values ==cs_z)[0]
                nlevs = np.shape(self.data[self.z_name].values)[0]
#            print (np.shape(self.data[self.lat_name]))
#            print 'q is '
            #print np.shape(self.data[self.z_name].sel(d=q))
#            zlev = np.where(self.data[self.z_name].sel(d=q).values ==2.25)[0]#self.cs_z)
            
#            print('cs zlev is',zlev)
#             print self.data[self.z_name].sel(d=q).values
#             print np.shape(self.data[self.dz_name])
#             print 'zlev',zlev
#            print('zlev for shy is:',zlev)
#            refl=np.squeeze(self.data[self.dz_name].sel(z=slice(zlev,zlev+1),d=q)).values
            #print('zlev',zlev)
            #refl=np.squeeze(self.data[self.dz_name].sel(z=zlev,d=q)).values
            refl = np.nanmax(np.squeeze(self.data[self.dz_name].sel(d=q).values),axis=0)
#            print('shape for shy:',np.shape(refl),np.shape(self.data[self.dz_name].values))
            if q==0:
                self.refl = refl
#            print ('refl shape',np.shape(refl))
            #print np.shape(self.data[self.dz_name].sel(d=slice(q,q+1),z=slice(zlev,zlev+1)))
            #refl = np.squeeze(self.data[self.dz_name].sel(d=q,z=zlev)).values
            #print 'refl shape',np.shape(refl),type(refl)
#             if len(np.shape(refl)) >= 3:
#                 print 'len is 3+'
#                 refl = np.squeeze(self.data[self.dz_name].sel(d=q,z=slice(zlev,zlev+1)).values)
#            refl.setflags(write=1)
#             refl[(np.isnan(refl))] = -99999.9
#             refl_missing_val = -99999.9
            #print self.data[self.dz_name].min()
 #            if self.lat_name in self.data.keys():
#                 lat = self.data[self.lat_name].values
#                 lon = self.data[self.lon_name].values
#             else:
#                 self.get_latlon_fromxy()
#                 lat = self.data[self.lat].values
#                 lon = self.data[self.lon].values
            #print np.shape(refl),'radar data 2490'
#             if np.ndim(lat)==1:
#                 lat2d = np.tile()            
#            print(np.shape(lat),np.shape(lon))
            if np.ndim(lat) == 1:
                lon2d = np.tile(lon,np.shape(lon)[0])
                lon2d = np.swapaxes(np.reshape(lon2d,(np.shape(lon)[0],np.shape(lon)[0])),1,0)

                lat2d = np.tile(lat,np.shape(lat)[0])
                lat2d = np.reshape(lat2d,(np.shape(lat)[0],np.shape(lat)[0]))
            else:
                lat2d = lat
                lon2d = lon

#            print(np.shape(lat2d),np.shape(lon2d))                
            yh_cs, yh_cc, yh_bkgnd = shy.conv_strat_latlon(refl, lat2d, lon2d, self.zconv, method='SYH', sm_rad=2,a=8, b=64)
            cs_arr = np.full(yh_cs.shape, np.nan)
            yh_conv = (yh_cs == 3) | (yh_cs == 4) | (yh_cs == 5)| (yh_cs == 2) | (yh_cs == 1)
            yh_strat = yh_cs == 0

            cs_arr[yh_conv] = 2
            cs_arr[yh_strat] = 1

#            cs_arr[np.isnan(refl)] =-1
    #        print nlevs
            rpt = np.tile(cs_arr,(nlevs,1,1))
            #print np.shape(rpt)
            #print np.nanmax(rpt)
            rntypetot.append(rpt)
            

        #self.def_convstrat()
        #np.array(rntypetot)[np.isnan(self.data[self.dz_name].values)] =-1
        #print 'shapes',np.shape(rntypetot)
        self.add_field((self.data[self.dz_name].dims,np.array(rntypetot)), 'CSS')
        mask=np.where(np.isnan(self.data[self.dz_name].values))
        hold = self.data['CSS'].values
        hold[mask] =np.nan
        self.data['CSS'].values = hold     
    
    def composite(self, var,ts=0,map_on = True,res='10m'):
        dat = np.squeeze(self.data[var].sel(d=slice(ts,ts+1)).values)
        comp = dat.fill(np.nan)
        dzcomp = np.nanmax(comp,axis=0)
        if not 'lat' in self.data.keys():
            self.get_latlon_fromxy()
        
        lats = self.data['lat']
        lons = self.data['lon']
    
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())
        from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
        ax.coastlines(resolution=res)
        if var in self.lims.keys():
            print( 'var:',var)
            range_lim = self.lims[var][1] - self.lims[var][0]
    #          print np.shape(data), np.shape(xdat),np.shape(ydat)
    #            print 'in var',var
            #print **kwargs
            vmin=self.lims[var][0]
            vmax=self.lims[var][1]
            cmap=self.cmaps[var]
        else:
            print ('unrecognized var',var)
            dat = self.data[var].data
            dat[dat<-900.0]=np.nan
            range_lim  = np.nanmax(dat) - np.nanmin(dat)
            vmin=np.nanmin(dat)
            vmax=np.nanmax(dat)
            cmap = plt.cm.gist_ncar
        
        ax.set_extent([np.min(self.lons), np.max(lons), np.min(lats), np.max(lats)])
        lon_formatter = LongitudeFormatter(number_format='.1f')
        lat_formatter = LatitudeFormatter(number_format='.1f')
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)

        cb = ax.pcolormesh(lons,lats,dzcomp, vmin = vmin, vmax = vmax, cmap = cmap,transform=ccrs.PlateCarree())
        plt.colorbar(cb)
        # ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        # ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=2, color='gray', alpha=0.5, linestyle='--')

        gl.xlabels_top = False
        gl.ylabels_right = False
        # ax.set_title('MC3E CSAPR {d:

    ######################################
    ##### plot_driver.PLOT_COMPOSITE #####
    ######################################

    # Description: plot_composite overlays a Cartopy basemap with a colormesh plot of a given variable.

    def plot_composite(self,var,time,resolution='10m',cs_over=False,statpt=False):
        
        import cartopy.crs as ccrs 
        
        dat = deepcopy(self.data[var].sel(d=time))
        if cs_over:
            whbad = np.where(self.data['CSS'].sel(d=time).values<0)
            dat.values[whbad] = np.nan
        arrdat = np.squeeze(dat.values)
        dzcomp = np.nanmax(arrdat,axis=0)

        if not self.lat_name in self.data.keys():
            print('No latitude. Calculating....')
            self.get_latlon_fromxy()
            lats = self.data['lat']
            lons = self.data['lon']
        else:
            if 'd' in self.data[self.lat_name].dims:
                lats = self.data[self.lat_name].sel(d=time).values
                lons = self.data[self.lon_name].sel(d=time).values
            else:
                lats = self.data[self.lat_name].values
                lons = self.data[self.lon_name].values
      
        minlon,maxlon,minlat,maxlat = lons[(0,0)],lons[(0,-1)],lats[(0,0)],lats[(-1,0)]

        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot(111, projection=ccrs.Mercator())
        
        cb = ax.pcolormesh(lons,lats,dzcomp,vmin=self.lims[var][0],vmax=self.lims[var][1],cmap=self.cmaps[var],transform=ccrs.PlateCarree())
       
        if cs_over: 
            cs_arr = np.squeeze(np.nanmax(np.squeeze(self.data['CSS'].sel(d=time).values),axis=0))
            ax.contour(lons,lats,cs_arr,levels=[0,1,2,3],linewidths=3,colors=['black','black'],transform=ccrs.PlateCarree())
            
        if statpt: ax.plot(self.lon_0,self.lat_0,markersize=12,marker='^',color='k',transform=ccrs.PlateCarree())

        self.co_gridlines(fig,ax,latlonyn=True,minx=minlon,maxx=maxlon,miny=minlat,maxy=maxlat)
        self.mycbar(fig,ax,cb,'Composite '+self.longnames[var]+' '+self.units[var])
        
        return fig, ax


    def co_gridlines(self,fig,ax,latlonyn=False,minx=-200,maxx=200,miny=-200,maxy=200,xlab=True,ylab=True,lnspc=2,ltspc=2,resolution='10m'):

        if latlonyn == True:

            import cartopy.crs as ccrs
            import matplotlib.ticker as ticker

            ax.coastlines(resolution=resolution)
            lur,bur,wur,hur = ax.get_position().bounds
          
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1.5, alpha=0.75, linestyle='--')
            gl.xlocator = ticker.MultipleLocator(base=lnspc)
            gl.ylocator = ticker.MultipleLocator(base=ltspc)
            gl.xlabels_top = False
            gl.ylabels_right = False 
            gl.x_inline = False
            gl.y_inline = False
            
            if xlab:
                gl.xlabel_style = {'size': 16, 'color': 'black'}
                #ax.text(0.5,bur-0.15,'Longitude',fontsize=16,ha='center',va='bottom',rotation_mode='anchor',transform=ax.transAxes)
            else:
                gl.xlabels_bottom = False
            if ylab:
                gl.ylabel_style = {'size': 16, 'color': 'black'} #, 'rotation': 90}
                #ax.text(lur-0.15,0.5,'Latitude',fontsize=16,ha='center',va='bottom',rotation_mode='anchor',rotation=90,transform=ax.transAxes)
            else:
                gl.xlabels_left = False
                gl.ylabel_style = {'size': 0}
                
            #newax = fig.add_axes(ax.get_position(), frameon=False)
            #newax.tick_params(axis='x', labelsize=0, length=0, pad=15)
            #newax.tick_params(axis='y', labelsize=0, length=0, pad=45)
            #if xlab: newax.set_xlabel('Longitude',fontsize=16)
            #if ylab: newax.set_ylabel('Latitude',fontsize=16)

            ax.set_extent([minx,maxx,miny,maxy],crs=ccrs.PlateCarree())
            ax.set_aspect('auto')

        else:

            import numpy as np
            
            ax.set_xlim([minx,maxx])
            ax.set_xticks(np.linspace(minx,maxx,5))
            ax.set_ylim([miny,maxy])
            ax.set_yticks(np.linspace(miny,maxy,5))

            ax.grid(axis='both',c='grey',linewidth=1.5, alpha=0.75, linestyle='--')

            if xlab:
                ax.set_xlabel('Distance E of radar (km)',fontsize=16)
                ax.tick_params(axis='x', which='major', labelsize=16)
            else:
                ax.set_xticklabels([])
                ax.xaxis.set_ticks_position('none')
            if ylab:
                ax.set_ylabel('Distance N of radar (km)',fontsize=16)
                ax.tick_params(axis='y', which='major', labelsize=16)
            else:
                ax.set_yticklabels([])
                ax.yaxis.set_ticks_position('none')


    def mycbar(self,fig,ax,cb,labtxt):        
        
        lur,bur,wur,hur = ax.get_position().bounds
        cbar_ax_dims = [lur+wur+0.015,bur,0.03,hur]
        cbar_ax = fig.add_axes(cbar_ax_dims)
        cbt = plt.colorbar(cb,cax=cbar_ax)
        cbt.ax.tick_params(labelsize=16)
        cbt.set_label(labtxt, fontsize=16, rotation=270, labelpad=20)

        return cbt
