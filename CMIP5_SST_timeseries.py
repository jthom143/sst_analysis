## Script to make a plot of the 50-70 average SSTs timeseries for the CMIP5 models

import iris
import iris.analysis.cartography
import iris.coord_categorisation
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import matplotlib as mpl


# Import my functions
import sys # access system routines
sys.path.append('/data1/fesd1/jthom143/python_functions')
from smooth import smooth
from sam_analysis_functions import sam_seasonal_averaging



# Load Data:
keylist = { 'ACCESS1-0': '/data1/fesd1/jthom143/SST_analysis/data/ACCESS1-0sst_5070_djf.nc',
            'ACCESS1-3': '/data1/fesd1/jthom143/SST_analysis/data/ACCESS1-3sst_5070_djf.nc',
            'bcc-csm1-1': '/data1/fesd1/jthom143/SST_analysis/data/bcc-csm1-1sst_5070_djf.nc',
            'CCSM4' : '/data1/fesd1/jthom143/SST_analysis/data/CCSM4sst_5070_djf.nc',
            'CNRM-CM5' : '/data1/fesd1/jthom143/SST_analysis/data/CNRM-CM5sst_5070_djf.nc',
            'CSIRO-Mk3-6-0':'/data1/fesd1/jthom143/SST_analysis/data/CSIRO-Mk3-6-0sst_5070_djf.nc',
            'GFDL-CM3':'/data1/fesd1/jthom143/SST_analysis/data/GFDL-CM3sst_5070_djf.nc',
            'GFDL-ESM2G':'/data1/fesd1/jthom143/SST_analysis/data/GFDL-ESM2Gsst_5070_djf.nc',
            'GFDL-ESM2M':'/data1/fesd1/jthom143/SST_analysis/data/GFDL-ESM2Msst_5070_djf.nc',
            'HadCM3':'/data1/fesd1/jthom143/SST_analysis/data/HadCM3sst_5070_djf.nc',
            'HadGEM2-AO':'/data1/fesd1/jthom143/SST_analysis/data/HadGEM2-AOsst_5070_djf.nc',
            'HadGEM2-CC':'/data1/fesd1/jthom143/SST_analysis/data/HadGEM2-CCsst_5070_djf.nc',
            'HadGEM2-ES': '/data1/fesd1/jthom143/SST_analysis/data/HadGEM2-ESsst_5070_djf.nc',
            'inmcm4': '/data1/fesd1/jthom143/SST_analysis/data/inmcm4sst_5070_djf.nc',
            #'MIROC5': '/data1/fesd1/jthom143/SST_analysis/data/MIROC5sst_5070_djf.nc',
            'MIROC-ESM': '/data1/fesd1/jthom143/SST_analysis/data/MIROC-ESMsst_5070_djf.nc',
            #'MRI-CGCM3': '/data1/fesd1/jthom143/SST_analysis/data/MRI-CGCM3sst_5070_djf.nc',
            'NorESM1-M': '/data1/fesd1/jthom143/SST_analysis/data/NorESM1-Msst_5070_djf.nc'
            }

era_int = iris.load_cube('/data1/fesd1/waugh/ERA_INTERIM/sst_197901_201406.nc')
era_djf, b, c, d = sam_seasonal_averaging(era_int)
era_djf = era_djf.collapsed('longitude', iris.analysis.MEAN)
S_5070 = iris.Constraint(latitude=lambda y: -70 < y < -50)
era_djf = era_djf.extract(S_5070)
era_djf = era_djf.collapsed('latitude', iris.analysis.MEAN)
era_djf.convert_units('celsius')
era_djf_mean = era_djf.collapsed('time', iris.analysis.MEAN)




sst = {}
length = {}
mean = {}
trend = {}
x = {}
sst_smooth = {}
detrend_sst = {}

n = 16
color = plt.cm.Paired(np.linspace(0.1,0.9,n)) # This returns RGBA; convert:
hexcolor = map(lambda rgb:'#%02x%02x%02x' % (rgb[0]*255,rgb[1]*255,rgb[2]*255),
               tuple(color[:,0:-1]))

mpl.rcParams['axes.color_cycle'] = hexcolor

fig1 = plt.figure(figsize=(10,6))
ax1 = plt.subplot(1,1,1)
ax1.set_xlabel('Time [years]')
ax1.set_ylabel('Degrees C')
ax1.set_title('CMIP 5 50-70S SST')
plt.axhline(era_djf_mean.data, color = 'k', ls = '--')



fig2 = plt.figure(figsize=(10,6))
ax2 = plt.subplot(1,1,1)
ax2.set_xlabel('Time [years]')
ax2.set_ylabel('Degrees C')
ax2.set_title('CMIP 5 50-70S SST')
plt.axhline(era_djf_mean.data, color = 'k', ls = '--')


fig3 = plt.figure(figsize=(12,5))
fs = 18
ax3 = fig3.add_axes([0.1, 0.2, 0.8, 0.7])
ax3.set_xlabel('Time [years]', fontsize = fs)
ax3.set_ylabel('SST [$^o$C]', fontsize = fs)
plt.xticks(fontsize = fs)
plt.yticks(fontsize = fs)
ax3.set_title('CMIP5 50-70S DJF SST', fontsize = fs)
#plt.axhline(era_djf_mean.data, color = 'k', ls = '--')

fig4 = plt.figure(figsize=(12,5))
ax4 = plt.subplot(1,1,1)
ax4.set_xlabel('Time [years]')
ax4.set_ylabel('Degrees C')
ax4.set_title('ERA Interim 50-70S DJF SST')
ax4.plot(era_djf.coord('season_year').points, era_djf.data, color = 'k', lw = 1.5)


for name, path in keylist.iteritems():

    sst[name] = iris.load_cube(path)
    
    sst[name].convert_units('celsius')

    length[name] = len(sst[name].coord('time').points)
    x[name] = np.arange(0, length[name])
    
    mean[name] = sst[name].collapsed('time', iris.analysis.MEAN).data
    
    # Quantify the drift in the models
    slope, b = np.polyfit(x[name], sst[name].data, 1)
    trend[name] = slope*length[name]
    
    smoothed = smooth(sst[name].data, window_len=11,window='flat')
    sst_smooth[name] = smoothed
    
    # Detrend smoothed
    detrend_sst[name] = sst_smooth[name] - slope*x[name]

    # Plot Timeseries
    ax1.plot(sst[name].data, lw = 1.5, label = name)
    ax1.legend(fontsize = 11)
    
    # Plot Smoothed Timeseries
    ax2.plot(sst_smooth[name], lw = 1.5, label = name)
    ax2.legend(fontsize = 11)


    # Plot Smoothed Timeseries
    ax3.plot(detrend_sst[name], lw = 3, label = name)
    #ax3.legend(fontsize = 11)

plt.show()


