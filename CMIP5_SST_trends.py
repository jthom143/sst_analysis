## Script to calculate the trends in SSTs

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
from sam_analysis_functions import sam_trend_analysis


# Load Data:
keylist = { 'bcc-csm1-1': '/data1/fesd1/jthom143/SST_analysis/data/bcc-csm1-1sst_5070_djf.nc',
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
#'MRI-CGCM3': '/data1/fesd1/jthom143/SST_analysis/data/MRI-CGCM3sst_5070_djf.nc'
            }

models = ['bcc-csm1-1', 'CSIRO-Mk3-6-0', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M', 'HadCM3', 'HadGEM2-AO', 'HadGEM2-CC', 'HadGEM2-ES', 'inmcm4', 'MIROC-ESM']

trend_period = 34

sst = {}
length = {}
x = {}
trends = {}
trends_mean = {}
trends_std = {}


for name, path in keylist.iteritems():
    
    sst[name] = iris.load_cube(path)
    
    sst[name].convert_units('celsius')
    
    length[name] = len(sst[name].coord('time').points)
    x[name] = np.arange(0, length[name])
    
    # Quantify the drift in the models
    slope, b = np.polyfit(x[name], sst[name].data, 1)
    
    # Detrend smoothed
    sst[name] = sst[name] - slope*x[name]
    
    trends[name], years = sam_trend_analysis(sst[name], trend_period, nonoverlapping=False)
    
    trends[name] = trends[name]/(trend_period/10)

    trends_mean[name] = np.mean(trends[name])
    trends_std[name] = np.std(trends[name])


# Create Figure
fig1 = plt.figure(figsize=(10,6))
ax1 = fig1.add_axes([0.1,0.25,0.8,0.6])
plt.xlim([0, 12])
plt.ylim([-1, 1])

x = range(1,12)

for i in range(0, 11):
    ax1.errorbar(i+1, trends_mean[models[i]], yerr=trends_std[models[i]]*2, fmt='o', color = 'b')
plt.xticks(x, models, rotation=45, ha='right')
plt.ylabel('Degrees C per decade')

''' 25-year trends ''' '''
plt.axhline(-0.096, ls = '--', color = 'g', label = 'Fan et al.')
plt.axhline(-0.11, ls = '--',color = 'r', label = 'EN3 Dataset')
plt.axhline(-0.052, ls = '--',color = '#800080' ,label = 'Met Office Reanalysis')
'''

''' 34-year trends '''
plt.axhline(-0.087, ls = '--', color = 'g', label = 'Fan et al.')
plt.axhline(-0.11, ls = '--',color = 'r', label = 'EN3 Dataset')
plt.axhline(-0.084, ls = '--',color = '#800080' ,label = 'Met Office Reanalysis')

plt.legend(fontsize = 11)
plt.title('%i year 50-70S SST trends'%trend_period)

plt.show()



