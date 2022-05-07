import numpy as np



################################################################################
################################################################################
#############    FUNCTIONS FOR ANALYZING DATA    ###############################
################################################################################
################################################################################

# records most intense element at each time step and its location
def highest_timeseries(data,lat,lon):

    #initialize arrays to hold the highest value and its location
    highest     = np.empty(data.shape[0])
    highest_lat = np.empty(data.shape[0])
    highest_lon = np.empty(data.shape[0])


    for i in range(data.shape[0]):

        highest[i]         = np.nanmax(data[i,:,:])
        # print(i)
        highest_lat[i]     = lat[np.where(data[i,:,:]==highest[i])[0][0]] #extra index turns lat value from array to scalar
        highest_lon[i]     = lon[np.where(data[i,:,:]==highest[i])[1][0]]

    return highest,highest_lat,highest_lon

# moving calculates the running average over a window size w
def moving_average(x,w):
    return np.convolve(x,np.ones(w),'valid')/w

#returns array shape of 1st dimension of passed data with the area of non-sentinel grid points at each timestep
def get_area_of_filtered_data(data,lat,lon,sentinel=np.nan):

    area_grid = get_area_grid(data,lat,lon)
    area_grid[np.where(np.isnan(data))]=np.nan

    area_timeseries = np.empty_like(data[:,0,0])
    for i in range(0,data.shape[0]):
        area_timeseries[i] = np.nansum(area_grid[i,:,:])

    return area_timeseries

# #returns 2d map with area of each lat lon box
# def get_area_grid(data,lat,lon):
#
#     area_grid = np.empty_like(data)
#
#     for i in range(0,area_grid.shape[1]):
#
#         area_grid[:,i,:]=(np.cos(np.radians(np.abs(lat[i])))*111.321)*(111*np.abs(lat[1]-lat[0]))
#
#     return area_grid

#returns 2d map with area of each lat lon box
def get_area_grid(data,lat,lon):

    area_grid = np.empty_like(data)

    for i in range(0,area_grid.shape[1]):

        area_grid[:,i,:]=(np.cos(np.radians(np.abs(lat[i])))*abs(lon[0]-lon[1])*111.322)*(111*np.abs(lat[1]-lat[0]))

    return area_grid

#computes a running average over an indicated window
def running_average(series, window,irregular_time=False,time=None):

    if irregular_time:
        if np.all(time==None):
            print('please provide corresponding times')
            return
        else:
            smoothed = np.ones(series)*np.nan;


    smoothed = np.empty_like(series);
    for i in range(0,series.size):
        if i>window/2:
            if i<series.size-window/2:
                smoothed[i]=np.mean(series[(i-int(window/2)):(i+int(window/2))])
            else:
                smoothed[i]=np.mean(series[(i-int(window/2)):])
        else:
            if i<series.size-window/2:
                smoothed[i]=np.mean(series[:(i+int(window/2))])
            else:
                smoothed[i]=np.mean(series)

    return smoothed

#bins timeseries into window size bins. counts number of extreme events for in each bin and returns array
#     return array has the number of extreme events the are in the window centered on that date
def extremes_per_window(date_range=None,dates_of_extremes=None,window=None):

    if np.any([np.all(date_range==None),np.all(dates_of_extremes==None),window==None]):
        print("function call: extreme_per_window,start_date = extremes_per_window(date_range=None,dates_of_extremes=None,window=None)")
        return
    if np.mod(window/2,1)==0:
        print('give odd sized window..')
        return

    start_index = int(window/2)
    start_date  = date_range[start_index]
    window_range = [0,window]
    end_index   = date_range.size-int(window/2)

    extremes_per_window = np.zeros(date_range.shape)
    extremes_per_window.fill(np.nan)
    for i in range(start_index,date_range.size-start_index):

        extremes_per_window[i] = dates_of_extremes[np.logical_and(dates_of_extremes>=(window_range[0]+(i-start_index)),dates_of_extremes<=(window_range[1]+(i-start_index)))].size

    return extremes_per_window,start_date

#helper for get center of mass. turns lat vector into an array that is the same size
# as shape with lat vactors in each column
def get_lat_grid(lat,shape):

    lat_grid = np.empty(shape)

    for i in range(0,shape[1]):
        lat_grid[:,i]=lat

    return lat_grid

#helper for get center of mass. Turns lon cevotr into array that is the same size
# as shape with lon in appropriate rows
def get_lon_grid(lon,shape):

    lon_grid = np.empty(shape)

    for i in range(0,shape[0]):
        lon_grid[i,:]=lon

    return lon_grid

#given field and latitude vectr return the latitude of the center of mass
def center_of_mass_lat(data,lat,weighted):

    if weighted:
        normalized_field = data/np.nanmean(data)
        return np.nansum(normalized_field*get_lat_grid(lat,data.shape))/(data[np.logical_not(np.isnan(data))].size)
    else:
        return np.nansum((data/data)*get_lat_grid(lat,data.shape))/(data[np.logical_not(np.isnan(data))].size)


#given a longitude vector and a data field, returns the lon of the center of mass
def center_of_mass_lon(data,lon,weighted):

    if weighted:
        normalized_field = data/np.nanmean(data)
        return np.nansum(normalized_field*get_lon_grid(lon,data.shape))/(data[np.logical_not(np.isnan(data))].size)
    else:
        return np.nansum((data/data)*get_lon_grid(lon,data.shape))/(data[np.logical_not(np.isnan(data))].size)

#given a 3D field (timeXlatXlon) finds the center of mass (a weighted average) at
# each timestep
def center_of_mass(data=None,lat=None,lon=None,weighted=True):

    center_of_mass = np.empty([data.shape[0],2])

    for i in range(0,data.shape[0]):

        center_of_mass[i,0]=center_of_mass_lat(data[i,:,:],lat,weighted)
        center_of_mass[i,1]=center_of_mass_lon(data[i,:,:],lon,weighted)

    return center_of_mass


