import numpy as np
import pandas as pd
import xarray as xr
import warnings
import mpl_toolkits.basemap as basemap


################################################################################
################################################################################
#############   FUNCTIONS FOR FETCHING DATA        #############################
################################################################################
################################################################################

#extract_var will extract the variable name data stored in filename
def extract_var(filename, variablename, multiple_files=False):

    if multiple_files:
        return xr.open_mfdataset(filename)[variablename].values
    else:
        return xr.open_dataset(filename)[variablename].values

#get_data will return lat, lon, and data of file when given filename and vaible filename
def get_data(filename,variablename,multiple_files=False,center_lon=True,conversion_ratio=1):

    #extract lat and lon
    lat = extract_var(filename,'lat',multiple_files)
    lon = extract_var(filename,'lon',multiple_files)

    #extract data
    data = extract_var(filename,variablename,multiple_files)

    #center lon and data on 0degrees east
    if center_lon:
        data,lon = fix_lon(data,lon)


    return data*conversion_ratio,lat,lon



#fix_lon will be given a data array and a longitude array
#        will return arrays of longitude and data that have been reformatted
def fix_lon(array_to_fix, lon_of_array):

    #make longitude from -180 to 180 degrees east
    temp_lon = lon_of_array
    temp_lon[np.where(temp_lon>180)]= temp_lon[np.where(temp_lon>180)]-360

    #reorient data to be centered at 0 degrees east
    i_west    = np.where(temp_lon<0)
    i_east    = np.where(temp_lon>=0)
    west      = temp_lon[i_west]
    east      = temp_lon[i_east]
    fixed_lon = np.array(np.hstack((west,east)))

    #make similar adjustments so that vwnd matches new longitude
    data_west   = np.squeeze(array_to_fix[:,:,i_west])
    data_east   = np.squeeze(array_to_fix[:,:,i_east])
    fixed_array = np.concatenate((data_west,data_east), axis=2)

    return fixed_array,fixed_lon

#########################################################################################################################
# Here are some functions for loading data so that I dont' have to keep writing/running scripts. The three data sets I'm
# interested in is the Ensemble Members DJF PSL; and Ensemble Members JJA PSL
#########################################################################################################################

def get_psl_DJF_EnsembleMembers():

    return np.load('/data/ncresswell/cesm-cam5-lme/EnsembleMembers_DJF_SLP_850-2005.npy'), np.load('/data/ncresswell/cesm-cam5-lme/lat.npy'), np.load('/data/ncresswell/cesm-cam5-lme/lon.npy')

def get_psl_JJA_EnsembleMembers():

    return np.load('/data/ncresswell/cesm-cam5-lme/EnsembleMembers_JJA_SLP_850-2005.npy'), np.load('/data/ncresswell/cesm-cam5-lme/lat.npy'), np.load('/data/ncresswell/cesm-cam5-lme/lon.npy')

def get_mean_psl_JJA():

    return get_mean_psl_seasonal('JJA')

def get_mean_psl_DJF():

    return get_mean_psl_seasonal('DJF')

def get_mean_psl_seasonal(season):
    psl = np.load('/data/ncresswell/cesm-cam5-lme/EnsembleMean_SLP_850-2005.npy')/100
    lat = np.load('/data/ncresswell/cesm-cam5-lme/lat.npy')
    lon = np.load('/data/ncresswell/cesm-cam5-lme/lon.npy')

    return seasonal_averages(psl,season=season)[:,:,:],lat,lon


################################################################################
################################################################################
#############    FUNCTIONS FOR FILTERING DATA    ###############################
################################################################################
################################################################################

#filterdata below or above threshold
def filter_tresh(data_unfiltered, threshold, below=False, above=False, sentinel=np.nan):

    data = np.copy(data_unfiltered)

    if below:
        data[np.where(data<threshold)]=sentinel
        return data
    elif above:
        data[np.where(data>threshold)]=sentinel
        return data
    else:
        print('No filtering occured. Please identify filter directions ("above" or "below")')
    return data

#filter data by lat/lon box
def filter_lat_lon(data_unfiltered,lat=None,lat_range=[-90,90],lon=None,lon_range=[-180,180],sentinel=np.nan, ndims=3):

    if ndims == 4:
        return filter_lat_lon_4D(data_unfiltered,lat,lat_range,lon,lon_range,sentinel)

    data = np.copy(data_unfiltered)


    if data.shape[1]!=lat.size:
        print('BAD: lat is not the first dimension...')
        return
    if data.shape[2]!=lon.size:
        print('BAD: lon is not the first dimension...')
        return

    #filter data outside this box
    data[:,np.where(lat<lat_range[0]),:]=sentinel
    data[:,np.where(lat>lat_range[1]),:]=sentinel

    data[:,:,np.where(lon<lon_range[0])]=sentinel
    data[:,:,np.where(lon>lon_range[1])]=sentinel

    return data

#filter data by lat/lon box data is 4D
def filter_lat_lon_4D(data_unfiltered,lat,lat_range,lon,lon_range,sentinel):

    data = np.copy(data_unfiltered)

    #filter data outside this box
    data[:,:,np.where(lat<lat_range[0]),:]=sentinel
    data[:,:,np.where(lat>lat_range[1]),:]=sentinel

    data[:,:,:,np.where(lon<lon_range[0])]=sentinel
    data[:,:,:,np.where(lon>lon_range[1])]=sentinel

    return data


def filter_land(data,lat,lon):

    bm = basemap.Basemap()

    for i in range(0,lat.size):
        for j in range(0,lon.size):

            if bm.is_land(lon[j],lat[i]):

                data[:,i,j] = np.nan

    return data

def get_land_filter(data_shape, lat, lon):

    filter = np.empty(data_shape)

    if data_shape[0]!=lat.size:
        print('====== INSIDE get_land_filter ======')
        print('data_shape format must be lat by lon')
        return

    bm = basemap.Basemap()

    for i in range(0,data_shape[0]):
        for j in range(0,data_shape[1]):

            if bm.is_land(lon[j],lat[i]):
                filter[i,j] = False
            else:
                filter[i,j] = True

    return filter

#############    filtering data into seasons/months    #########################

#get_DJF will filter DJF data and either return the filtered data or calculate a seasonal mean
def get_DJF(data,dims):

    #create index arrays for dec, jan, feb
    dec_i = np.arange(11,data.shape[0],12)
    jan_i = np.arange(0,data.shape[0],12)
    feb_i = np.arange(1,data.shape[0],12)

    #initialize and populate seasonal index array
    djf_i = np.empty(dec_i.size+jan_i.size+feb_i.size,dtype=int)

    djf_i[::3]=jan_i
    djf_i[1:][::3]=feb_i
    djf_i[2::][::3]=dec_i

    if dims==3:
        #filter with djf indices
        return data[djf_i,:,:]
    elif dims==1:
        return data[djf_i]

#get_MAM will filter MAM data and return the filtered data
def get_MAM(data,dims):

    #create index arrays for dec, jan, feb
    mar_i = np.arange(2,data.shape[0],12)
    apr_i = np.arange(3,data.shape[0],12)
    may_i = np.arange(4,data.shape[0],12)

    #initialize and populate seasonal index array
    mam_i = np.empty(mar_i.size+apr_i.size+may_i.size,dtype=int)

    mam_i[::3]=mar_i
    mam_i[1:][::3]=apr_i
    mam_i[2::][::3]=may_i

    #filter with mam indices
    if dims==3:
        #filter with mam indices
        return data[mam_i,:,:]
    elif dims==1:
        return data[mam_i]

#get_JJA will filter JJA data and return filtered data
def get_JJA(data,dims):

    #create index arrays for dec, jan, feb
    jun_i = np.arange(5,data.shape[0],12)
    jul_i = np.arange(6,data.shape[0],12)
    aug_i = np.arange(7,data.shape[0],12)

    #initialize and populate seasonal index array
    jja_i = np.empty(jun_i.size+jul_i.size+aug_i.size,dtype=int)

    jja_i[::3]=jun_i
    jja_i[1:][::3]=jul_i
    jja_i[2::][::3]=aug_i

    #filter with jja indices
    if dims==3:
        #filter with jja indices
        return data[jja_i,:,:]
    elif dims==1:
        return data[jja_i]

#get_SON will filter SON data and return filtered data
def get_SON(data,dims):

    #create index arrays for dec, jan, feb
    sep_i = np.arange(8,data.shape[0],12)
    oct_i = np.arange(9,data.shape[0],12)
    nov_i = np.arange(10,data.shape[0],12)

    #initialize and populate seasonal index array
    son_i = np.empty(sep_i.size+oct_i.size+nov_i.size,dtype=int)

    son_i[::3]=sep_i
    son_i[1:][::3]=oct_i
    son_i[2::][::3]=nov_i

    if dims==3:
        #filter with son indices
        return data[son_i,:,:]
    elif dims==1:
        return data[son_i]

#break data into time intervals
def break_field_by_interval(data,interval,from_beginning=False,from_end=False):

    #initialize
    broken_data = np.empty([int(np.floor(data.shape[0]/interval)),int(interval),data.shape[1],data.shape[2]])
    num_intervals = np.floor(data.shape[0]/interval)

    if from_end:

        broken_data = data[-(int(num_intervals*interval)):,:,:]
        return broken_data.reshape([int(num_intervals),int(interval),broken_data.shape[1],broken_data.shape[2]])

    elif from_beginning:

        broken_data = data[:(int(num_intervals*interval)),:,:]
        return broken_data.reshape([int(num_intervals),interval,broken_data.shape[1],broken_data.shape[2]])

    else:
        print('Specifcy whether to filter from beginning or from end')
        return



##################    break into predetermined intervals    ####################

#break time series into intervals
def break_series_by_interval(series,interval,from_beginning=False,from_end=False):

    #initialize
    num_intervals = int(np.floor(series.shape[0]/interval))
    broken_series = np.empty([num_intervals,int(interval)])

    if from_end:

        broken_series = series[-(int(num_intervals*interval)):]
        return broken_series.reshape([num_intervals,int(interval)])

    elif from_beginning:

        broken_series = series[:(int(num_intervals*interval))]
        return broken_series.reshape([num_intervals,interval])

    else:
        print('Specifcy whether to filter from beginning or from end')
        return

################################################################################
################################################################################
#############    FUNCTIONS FOR AVERAGING DATA     ##############################
################################################################################
################################################################################

#calculate seasonal averages
def seasonal_averages(data, season=None, dims=3):

    #filter data for a indicated season
    if   season == 'DJF':
        return seasonal_averages_DJF(get_DJF(data,dims=dims),dims=dims) #DJF is treated differencely as it spans two years
    elif season == 'MAM':
        seasonal_data = get_MAM(data,dims=dims)
    elif season == 'JJA':
        seasonal_data = get_JJA(data,dims=dims)
    elif season == 'SON':
        seasonal_data = get_SON(data,dims=dims)
    else:
        print('please specify season...')
        return

    if dims == 3:
        #reshape seasonal data into years
        yearly_binned_seasonal_data = np.reshape(seasonal_data,[int(seasonal_data.shape[0]/3),3,seasonal_data.shape[1],seasonal_data.shape[2]])
    elif dims == 1:
        #reshape seasonal data into years
        yearly_binned_seasonal_data = np.reshape(seasonal_data,[int(seasonal_data.shape[0]/3),3])
    else:
        print('number of dimensions not supported')
        return

    #average for each year and return
    return np.nanmean(yearly_binned_seasonal_data,axis=1)


################################################################################
#seasonal averages for DJF
def seasonal_averages_DJF(seasonal_data,dims=3):

    if dims==3:
        #initialize array to hold seasonal average for djf
        yearly_averaged_seasonal_data = np.empty([int(seasonal_data.shape[0]/3),seasonal_data.shape[1],seasonal_data.shape[2]])

        #average first jan and feb
        yearly_averaged_seasonal_data[0,:,:]=np.nanmean(seasonal_data[0:2,:,:],axis=0)

        #compute yearly seasonal averages for other years using reshaping to bin
        seasonal_data_offset               =  seasonal_data[2:seasonal_data.shape[0]-1,:,:]
        yearly_binned_seasonal_data_offset = np.reshape(seasonal_data_offset,[int(seasonal_data_offset.shape[0]/3),3,seasonal_data_offset.shape[1],seasonal_data_offset.shape[2]])
        yearly_averaged_seasonal_data_offset = np.nanmean(yearly_binned_seasonal_data_offset,axis=1)

        #populate the rest of yearly_averaged_seasonal_data with the calculated means
        yearly_averaged_seasonal_data[1:,:,:]=yearly_averaged_seasonal_data_offset

    elif dims==1:
        #initialize array to hold seasonal average for djf
        yearly_averaged_seasonal_data = np.empty([int(seasonal_data.shape[0]/3)])

        #average first jan and feb
        yearly_averaged_seasonal_data[0]=np.nanmean(seasonal_data[0:2],axis=0)

        #compute yearly seasonal averages for other years using reshaping to bin
        seasonal_data_offset               =  seasonal_data[2:seasonal_data.shape[0]-1]
        yearly_binned_seasonal_data_offset = np.reshape(seasonal_data_offset,[int(seasonal_data_offset.shape[0]/3),3])
        yearly_averaged_seasonal_data_offset = np.nanmean(yearly_binned_seasonal_data_offset,axis=1)

        #populate the rest of yearly_averaged_seasonal_data with the calculated means
        yearly_averaged_seasonal_data[1:]=yearly_averaged_seasonal_data_offset

    return yearly_averaged_seasonal_data


def get_threshold(filename=None, var_name=None,open_as_data_array=None,fix_lon=True,calculate_seasonal_average=True, region=None, off_mean=.5,conversion=1):
        
    """
    This function will be used to calculate the threshold of the Azores High for a given file containing PSL. 
    """
    
    #load data and center longitude
    if open_as_data_array:
        psl = xr.open_dataarray(filename)
    else:
        psl = xr.open_mfdataset(filename)[var_name]
    if fix_lon:
        psl = psl.assign_coords(lon=xr.where(psl.lon<180,psl.lon,psl.lon-360))

    #filter region defined above
    psl = psl.isel(lat=((psl.lat >= region['lat_min']) & (psl.lat <= region['lat_max'])))
    psl = psl.isel(lon=((psl.lon >= region['lat_min']) & (psl.lon <= region['lat_max'])))

    #find DJF average
    if calculate_seasonal_average:
        psl_djf = seasonal_averages(psl.values,season='DJF')
        return (np.nanmean(psl_djf) + off_mean * np.nanstd(psl_djf.flatten())) * conversion
    else:
        return ((np.nanmean(psl) + off_mean * psl.std()).values * conversion)