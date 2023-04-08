# this houses the RadarConfig object which is basically just carrying
# along some options about names, colorbars, plotting limits, etc.

from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import datetime as dt
import pyart

class RadarConfig(object):
   
    def __init__(self, dz='DZ', zdr='DR', kdp='KD', ldr='LH', rho='RH', hid = 'HID',conv='Con',temp='T', x='x', y='y', z='z', z3d='z3d', eml='eML', u='u', v='v',rr='RR', w='w',qc='qc',qr='qr',qg='qg',qh='qh',qi='qi',qs='qs',vr='VR',exper = 'Case',rtype='obs',rsrc='nexrad',band = 'C',lat_0 = 0,lon_0=90.0,lat_r=None,lon_r=None,lat=None,lon=None,tm = None,radar_name = None,color_blind = False,dd_on=False,rr_on=False,hid_on=True,hid_cats = 'summer'):
        
        # ******** first the polarimetric stuff *************
       
        if rtype.startswith('obs'):

            if rsrc.startswith('nexrad'):
        
                self.dz_name = 'REF' # Name of the reflectivity field
                self.zdr_name = 'ZDR' # Name of the differential reflectivity field
                self.kdp_name = 'KDP' # Name of the Kdp field
                self.rho_name = 'RHO' # Name of the RhoHV field
                self.vr_name = 'VEL' # Name of radial velocity field
                self.temp_name = 'T' # Name of the temperature field
                self.lon_name = 'lon0' # File naming of the longitude variable 
                self.lat_name = 'lat0' # File naming of the latitude variable
            
            if rsrc.startswith('cwr'):
        
                self.dz_name = 'DBZH' # Name of the reflectivity field
                self.zdr_name = 'ZDR' # Name of the differential reflectivity field
                self.kdp_name = 'KDP' # Name of the Kdp field
                self.rho_name = 'RHOHV' # Name of the RhoHV field
                self.vr_name = 'VRADH' # Name of radial velocity field
                self.temp_name = 'T' # Name of the temperature field
                self.lon_name = 'lon' # File naming of the longitude variable
                self.lat_name = 'lat' # File naming of the latitude variable

        elif rtype.startswith('wrf'):
        
            self.dz_name = 'zhh01' # Name of the reflectivity field
            self.zdr_name = 'zdr01' # Name of the differential reflectivity field
            self.kdp_name = 'kdp01' # Name of the Kdp field
            self.rho_name = 'rhohv01' # Name of the RhoHV field
            self.vr_name = 'vrad01' # Name of radial velocity field
            self.temp_name = 't_air' # Name of the temperature field
            self.lon_name = 'longitude' # File naming of the longitude variable 
            self.lat_name = 'latitude' # File naming of the latitude variable

        else:    
        
            self.dz_name = dz
            self.zdr_name = zdr
            self.kdp_name = kdp
            self.rho_name = rho
            self.vr_name = vr
            self.temp_name = temp
            self.lat_name=lat
            self.lon_name=lon
        
        print('in config, lat name is ',self.lat_name)
 
        self.u_name = u # Name of the zonal wind field
        self.v_name = v # Name of the meridional wind field
        self.w_name = w # Name of the vertical wind field
        self.x_name = x # File naming of the zonal directional variable
        self.y_name = y # File naming of the meridional directional variables
        self.z_name = z # File naming of the vertical level field

        self.qc_name = qc
        self.qr_name = qr
        self.qg_name = qg
        self.qh_name = qh
        self.qi_name = qi
        self.qs_name = qs

        self.conv_name = conv
        self.cs_name = 'CSS'
        self.rr_name = rr
        self.ldr_name = 'LDR' # Name of the linear depolarization field
        self.hid_name = 'HID'
        if self.rr_name == None:
            self.rr_name = 'RR'
        self.hid_on = hid_on
        self.dd_on = dd_on
        self.rr_on = rr_on

        ########set up some experiment parameters ############
        self.exper=exper
        self.radar_lat = lat_r
        self.radar_lon = lon_r
        self.lat_0 = lat_0
        self.lon_0 = lon_0
        
        self.band = band
        self.expr =exper
        self.date = tm
        self.radar_name = radar_name
        
        if hid_cats.startswith('winter'):     
            #self.species = np.array(['DZ','RN','IC','AG','WS','PL','DR'])
            self.species = np.array(['IC','PL','DR','AG','WS','DZ','RN'])
            #self.species_long = np.array(['Drizzle','Rain','Ice\nCrystals','Snow\nAggregates','Wet\nSnow','Plates','Dendrites'])
            self.species_long = np.array(['Ice\nCrystals','Plates','Dendrites','Snow\nAggregates','Wet\nSnow','Frozen\nPrecipitation','Rain'])
            #self.hid_colors = ['LightBlue','Blue','DarkOrange','Pink','Cyan','Purple','Fuchsia']
            self.hid_colors = ['DarkOrange','Purple','Fuchsia','Pink','Cyan','LightBlue', 'Blue']
        else:
            self.species = np.array(['DZ','RN','CR','AG','WS','VI','LDG','HDG','HA','BD'])
            self.species_long = np.array(['Drizzle','Rain','Ice\nCrystals','Snow\nAggre-\ngates','Wet\nSnow','Vertical\nIce','Low-\nDensity\nGraupel','High-\nDensity\nGraupel','Hail','Big\nDrops'])
        #self.hid_colors = ['White','LightBlue','MediumBlue','Darkorange','LightPink','Cyan','DarkGray',\
        #    'Lime','Yellow','Red','Fuchsia']
            self.hid_colors = ['LightBlue','MediumBlue','Darkorange','Pink','Cyan','DarkGray','Lime','Yellow','Red','Fuchsia'] 

        self.hidwater = [1,2,10] # Group drizzle, rain and big drops
        self.hidgraup = [7,8] # Group low and high density graupel
        self.hidhail = [9] # Hail
        self.hidsnow = [3,4,5,6] # Group ice crystals, snow, wet snow and VI

        self.pol_vars = np.array([self.dz_name, self.zdr_name, self.kdp_name, self.ldr_name, self.rho_name, self.hid_name])

        self.cs_colors = ['#FFFFFF', 'DodgerBlue', 'Red', 'Khaki']
        self.cs_labels = ['', 'Strat', 'Conv', 'Mixed']

        self.cfad_levs = [0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,15.0,20.,25.]
        self.cfad_cols = ['silver','darkgray','slategrey','dimgray','blue','mediumaquamarine','yellow','orange','red','fuchsia','violet'] 

        self.set_dbz_colorbar(color_blind=color_blind)
        self.set_hid_colorbar()
        self.set_cs_colorbar()

        self.rhi_vars = [self.dz_name,self.zdr_name,self.kdp_name,self.rho_name,self.hid_name if self.hid_on else None,self.w_name if self.dd_on else None] #Names of vars for RHI plots
        self.cfad_vars = [self.dz_name,self.zdr_name,self.kdp_name,self.rho_name,self.w_name if self.dd_on else None,self.hid_name if self.hid_on else None] #Names of vars for RHI plots
        self.cfad_compare_vars = [self.hid_name if self.hid_on else None,self.dz_name,self.zdr_name,self.kdp_name,self.rho_name,self.w_name if self.dd_on else None] #Names of vars for RHI plots
        self.cappi_vars = [self.dz_name,self.zdr_name,self.kdp_name,self.rho_name,self.vr_name,self.hid_name if self.hid_on else None] #Names of vars for RHI plots
        #self.q_vars = [self.qc_name,self.qr_name,self.qg_name,self.qh_name,self.qi_name,self.rr_name if self.rr_on else None,self.hid_name if self.hid_on else None] #Mixing ratios from model to plot.
        self.q_vars = [self.qc_name,self.qr_name,self.qs_name,self.qg_name,self.qi_name,self.hid_name if self.hid_on else None] #Mixing ratios from model to plot.
        self.water_vars = [self.qr_name,self.rr_name if self.rr_on else None,self.hid_name if self.hid_on else None] #Mixing ratios from model to plot.

        # Now just set some defaults
        self.lims = {self.dz_name: [0,80], self.zdr_name: [-1, 3], self.kdp_name: [-0.5, 3], self.ldr_name: [-35, -20], self.rho_name: [0.95, 1.00], self.hid_name: [0,len(self.species)+1],self.u_name: [-20,80],self.w_name:[-2,2],self.vr_name:[-25,25],self.cs_name:[0,4],self.rr_name:[0,30],self.temp_name:[-30,30],self.qc_name:[0.0,2.0],self.qr_name:[0.0,2.0],self.qg_name:[0.0,2.0],self.qh_name:[0.0,2.0],self.qi_name:[0.0,2.0],self.qs_name:[0.0,2.0]}
        self.cfbins = {self.dz_name: np.arange(-10,60,1), self.zdr_name: np.arange(-2,6,0.05), self.kdp_name: np.arange(-2,2,0.05), self.rho_name: np.arange(0.5,1.01,0.02), self.hid_name: '' , self.u_name: np.arange(-20,80,2), self.w_name: np.arange(-2,2,0.5), self.temp_name: np.arange(20,-60,-5)} 
        self.cfticks = {self.dz_name: np.linspace(-10,60,8), self.zdr_name: np.linspace(-2,6,9), self.kdp_name: np.linspace(-2,2,5), self.rho_name: np.linspace(0.95,1,6), self.hid_name: '' , self.u_name: np.linspace(-20,80,6), self.w_name: np.linspace(-6,6,7), self.temp_name: np.linspace(20,-60,5)} 
        self.delta = {self.dz_name: 10, self.zdr_name: 1, self.kdp_name: 1, self.ldr_name: 5, self.rho_name: 0.005, self.hid_name: 1,self.w_name:5,self.vr_name:5,self.cs_name:1,self.rr_name:10,self.temp_name:5}
        self.units = {self.dz_name: '(dBZ)', self.zdr_name: '(dB)', self.kdp_name: '($^{\circ}$ km$^{-1}$)', self.ldr_name: '(dB)', self.rho_name: '', self.hid_name: '', self.u_name: '(m s$^{-1}$)', self.w_name:'(m s$^{-1}$)',self.vr_name:'(m s$^{-1}$)',self.cs_name:'',self.rr_name:'(mm hr$^{-1}$)',self.temp_name:'($^{\circ}$C)',self.qc_name:'(g kg$^{-1}$)',self.qr_name:'(g kg$^{-1}$)',self.qg_name:'(g kg$^{-1}$)',self.qh_name:'(g kg$^{-1}$)',self.qi_name:'(g kg$^{-1}$)',self.qs_name:'(g kg$^{-1}$)'}
        self.names = {self.dz_name: 'Z', self.zdr_name: 'Z$_{DR}$', self.kdp_name: 'K$_{dp}$', self.ldr_name: 'LDR', self.rho_name: r'$\rho_{hv}$', self.hid_name: '',self.w_name:'',self.vr_name:'V$_r$',self.cs_name:'',self.rr_name:'RR',self.temp_name:'T',self.qc_name:'Q$_{C}$',self.qr_name:'Q$_{R}$',self.qg_name:'Q$_{G}$',self.qh_name:'Q$_{H}$',self.qi_name:'Q$_{I}$',self.qi_name:'Q$_{S}$'}
        self.names_uc = {self.dz_name: 'REF', self.zdr_name: 'ZDR', self.kdp_name: 'KDP', self.ldr_name: 'LDR', self.rho_name: 'RHO', self.hid_name: 'HID', self.u_name:'U', self.w_name:'W',self.vr_name:'VRAD',self.cs_name:'',self.rr_name:'RR',self.temp_name:'TEMP',self.qc_name:'QCLOUD',self.qr_name:'QRAIN',self.qg_name:'QGRAUPEL',self.qh_name:'QHAIL',self.qi_name:'QICE',self.qs_name:'QSNOW'}
        self.longnames = {self.dz_name: 'Reflectivity', self.zdr_name: 'Differential Reflectivity', self.kdp_name: 'Specific Differential Phase',self.ldr_name: 'Linear Depolarization Ratio', self.rho_name: 'Correlation Coefficient', self.hid_name: 'Hydrometeor Identification',self.w_name:'Vertical Velocity',self.vr_name:'Radial Velocity',self.cs_name: 'Convective/Stratiform',self.rr_name:'Rain Rate',self.temp_name:'Temperature',self.qc_name:'Cloud Water Mixing Ratio',self.qr_name:'Rain Water Mixing Ratio',self.qg_name:'Graupel Mixing Ratio',self.qh_name:'Hail Mixing Ratio',self.qi_name:'Ice Mixing Ratio',self.qs_name:'Snow Mixing Ratio'}
        #self.cmaps = {self.dz_name: self.temp_cmap, self.zdr_name: plt.cm.Spectral_r, self.kdp_name: plt.cm.gist_heat_r, self.ldr_name: plt.cm.gist_rainbow_r, self.rho_name: plt.cm.jet, self.hid_name: self.hid_cmap,self.w_name:plt.cm.seismic,self.vr_name:plt.cm.bwr,self.cs_name: self.cs_cmap,self.rr_name:plt.cm.Spectral_r,self.temp_name:'RdYlBu_r'}
        self.cmaps = {self.dz_name: self.temp_cmap, self.zdr_name: plt.cm.Spectral_r, self.kdp_name: plt.cm.gist_heat_r, self.ldr_name: plt.cm.gist_rainbow_r, self.rho_name: plt.cm.jet, self.hid_name: self.hid_cmap,self.u_name: plt.cm.seismic, self.w_name:plt.cm.seismic,self.vr_name:plt.cm.bwr,self.cs_name: self.cs_cmap,self.rr_name:plt.get_cmap('pyart_HomeyerRainbow'),self.temp_name:'RdYlBu_r',self.qc_name:plt.cm.gist_ncar,self.qr_name:plt.cm.gist_ncar,self.qg_name:plt.cm.gist_ncar,self.qh_name:plt.cm.gist_ncar,self.qi_name:plt.cm.gist_ncar,self.qs_name:plt.cm.gist_ncar}
        self.ticklabels = {self.dz_name: np.arange(0, 90, 10), self.zdr_name: np.arange(-1, 4, 1), self.kdp_name: np.arange(-0.5, 4.5, 1),self.ldr_name: np.arange(-35, -15, 5), self.rho_name: np.arange(0.95, 1.01, 0.005), self.hid_name: np.append('', self.species),self.u_name: np.arange(-20,81,2), self.w_name:np.arange(-2,2.1,0.1),self.vr_name:np.arange(-25,30.0,5.0),self.cs_name: self.cs_labels,self.rr_name:[0.1,1,10,30,50,70,100,130,150],self.temp_name:np.arange(-30,35,5),self.qc_name:np.arange(0.0,2.1,0.25),self.qr_name:np.arange(0.0,2.1,0.25),self.qg_name:np.arange(0.0,2.1,0.25),self.qh_name:np.arange(0.0,2.1,0.25),self.qi_name:np.arange(0.0,2.1,0.25),self.qs_name:np.arange(0.0,2.1,0.25)}
#############################################################################################################

    def print_date(self,tm=None, fmt='%Y-%m-%d %H:%M:%S %Z'):
#        print tm
        if tm is not None:
            #print tm
            if len(tm) > 1:
                date1 = tm[0]
                date2 = tm[-1]
                tms = dt.datetime.strftime(date1,fmt)
                tme = dt.datetime.strftime(date2,fmt)
                tmf = '{s}-{e}'.format(s=tms,e=tme)
            else:
                tmf = dt.datetime.strftime(tm[0],fmt)
        else:
            tmf = dt.datetime.strftime(np.array(self.date)[0],fmt)
        #print tmf
        return tmf

#############################################################################################################
    def print_title(self,tm = None):
        extra = ''
        if self.exper is not None:
            extra = '{e} {x}'.format(e=extra,x = self.exper)
        if self.radar_name is not None:
            extra = '{e} {r}'.format(e=extra,r=self.radar_name)
        if self.mphys is not None:
            extra = '{e} {m}'.format(e=extra,m=self.mphys)
        
        extra = '{e}, {t}'.format(e=extra, t=self.print_date(tm))
        return extra
#############################################################################################################

    def sav_title(self,tm = None):
        extra = ''
        if self.exper is not None:
            extra = '{e}_{x}'.format(e=extra,x = self.exper)
#        if self.radar_name is not None:
#            extra = '{e}_{r}'.format(e=extra,r=self.radar_name)
        if self.mphys is not None:
            extra = '{e}_{m}'.format(e=extra,m=self.mphys)
        
#        extra = '{e}_{t}'.format(e=extra, t=self.print_date(tm,fmt = '%Y%m%d%_H%M%S'))
        return extra
#############################################################################################################

    def set_dbz_colorbar(self, color_list=None, color_blind=False):
        if color_list is None:
            # just use the default here
            if color_blind is not True:
                radarcbar = ['PeachPuff','Aqua','DodgerBlue','Blue','Lime', \
                    'LimeGreen','Green','Yellow','Orange','DarkOrange','Red', \
                    'Crimson','Fuchsia','Purple','Indigo','MidnightBlue'] 
            #radarcbar = ['PeachPuff','Aqua','DodgerBlue','MediumBlue','Lime', \
            #    'LimeGreen','Green','Yellow','Orange','OrangeRed','Red', \
            #    'Crimson','Fuchsia','Indigo','DarkCyan','White'] 
            else:
                radarcbar = ['Lavender', 'Thistle', 'Plum', 'MediumPurple', 'CornFlowerBlue', 'SkyBlue', 'PaleTurquoise', 'LightCyan', 'Yellow', 'Gold', 'Orange', 'DarkOrange', 'Chocolate', 'IndianRed', 'FireBrick', 'Maroon']
        else: 
            radarcbar = deepcopy(color_list)

        temp_cmap = colors.ListedColormap(radarcbar)
        self.temp_cmap = temp_cmap # export for use with trajectory object


#############################################################################################################


    def set_hid_colorbar(self, color_list=None):

        if color_list is None:
         hidcbar = deepcopy(self.hid_colors)

        else:
            hidcbar = deepcopy(color_list)
        self.hid_cmap = colors.ListedColormap(hidcbar)

        #self.boundshid = np.arange(0,12)
        self.boundshid = np.arange(self.hid_cmap.N+1)
        self.normhid = colors.BoundaryNorm(self.boundshid, self.hid_cmap.N)
#############################################################################################################
    def set_cs_colorbar(self, color_list=None):
        if color_list is None:
            cscbar = deepcopy(self.cs_colors)

        else:
            cscbar = deepcopy(color_list)

        self.cs_cmap = colors.ListedColormap(cscbar)
        self.cs_bounds = np.arange(0,5)
        self.cs_norm = colors.BoundaryNorm(self.cs_bounds, self.cs_cmap.N)

#############################################################################################################

    def set_lims(self, var_key, varlims):
        "Reset limits of a variable for plotting, feed it a 2 element list"
        self.lims[var_key] = varlims
#############################################################################################################

    def get_lims(self, var_key):
        return self.lims[var_key]
#############################################################################################################

    def set_delta(self, var_key, new_delta):
        self.delta[var_key] = new_delta
#############################################################################################################

    def get_delta(self, var_key):
        return self.delta[var_key]
#####################################
