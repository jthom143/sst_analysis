# Import Functions
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from scipy import interpolate, stats
from iris import coord_categorisation
from scipy.stats import gaussian_kde
import numpy as np
import glob

# Import my functions
import sys # access system routines
sys.path.append('/data1/fesd1/jthom143/python_functions')
from cmip5_fix import latlon_fix


## Load in Data

""" Models with multiple files """

keylist = {#'inmcm4': 'sea_surface_temperature', ** Completed **
           #'MRI-CGCM3':'sea_surface_temperature', ** Completed **
           #'MIROC-ESM':'sea_surface_temperature', ** Completed **
           #'MIROC5':'sea_surface_temperature', ** Completed **
           #'HadGEM2-ES':'sea_surface_temperature', ** Completed **
           #'HadGEM2-CC':'sea_surface_temperature', ** Completed **
           #'HadCM3':'sea_surface_temperature', ** Completed **
           #'GFDL-ESM2M':'surface_temperature', ** Completed **
           #'GFDL-ESM2G':'sea_surface_temperature' ** Completed **
           #'GFDL-CM3': 'surface_temperature', ** Completed **

           #'ACCESS1-3':'sea_surface_temperature', *** Completed ***
           #'ACCESS1-0':'sea_surface_temperature'  *** Completed ***
           #'CNRM-CM5':'sea_surface_temperature'   *** Completed ***
           #'CCSM4':'sea_surface_temperature'      *** Completed ***
           #'NorESM1-M':'sea_surface_temperature', *** Completed ***
            'MPI-ESM-MR':'sea_surface_temperature'

}
# Problem with coordinates:

#'MPI-ESM-MR':'sea_surface_temperature',
#'MPI-ESM-LR':'sea_surface_temperature',

# Problem with time:
#'GISS-E2-H':'sea_surface_temperature',
#'GISS-E2-R':'sea_surface_temperature',
#'CMCC-CM':'sea_surface_temperature'





for key, var in keylist.iteritems():
    print 'working on ' +key
    

    files = sorted(glob.glob('/data1/fesd1/uhaus/tos_Omon_'+key+'_piControl_r1i1p1_*.nc'))
    sst = iris.load(files, var)
    
    for isst in sst:
        isst.attributes = {}
        #isst.remove_coord(isst.aux_coords[0])

    sst = sst.concatenate_cube()

    ## If only one file:
    #file = '/data1/fesd1/uhaus/tos_Omon_'+key+'_piControl_r1i1p1_*.nc'
    #sst = iris.load_cube(file, var)
    #sst.remove_coord(sst.aux_coords[0])
    #sst.remove_coord(sst.aux_coords[0])

    # Create seasons
    iris.coord_categorisation.add_season(sst, 'time', name = 'clim_season')
    iris.coord_categorisation.add_season_year(sst, 'time', name = 'season_year')
    
    # Average over each season
    sst_seasons = sst.aggregated_by(['clim_season', 'season_year'], iris.analysis.MEAN)
    
    # Make sure each season has 3 months
    spans_three_months = lambda t: (t.bound[1]-t.bound[0] > 3)
    three_months_bound = iris.Constraint(time = spans_three_months)
    sst_seasons_mean = sst_seasons.extract(three_months_bound)
    
    # Isolate DJF
    djf = iris.Constraint(clim_season = 'djf')
    sst_djf = sst_seasons_mean.extract(djf)

    ## If longitude and latitude are no dimensions:
    sst_djf = latlon_fix(sst_djf, key)

    stop

    # Average over 50-70S
    lats = iris.Constraint(latitude=lambda y: -70 < y < -50)
    sst_SH = sst_djf.extract(lats)

    sst_5070 = sst_SH.collapsed(['longitude','latitude'], iris.analysis.MEAN)

    iris.save(sst_5070, "/data1/fesd1/jthom143/SST_analysis/data/"+key+"sst_5070_djf.nc")
