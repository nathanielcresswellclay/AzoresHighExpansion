import numpy as np 

"""
This function will plot frequency of extreme events as reported by ERA-20C pressure data over a defined window

USAGE: to use pass an index that is defined through time, and specify what the threshold is for extreme event using the
      'percent' perameter.

"""

def running_bin_extremes(index=None, percent=20, window=11, time=np.arange(1850,2006), label_large = '',label_small = '',title = '',legend_loc=2, index_has_nan=False,ylim=None):

    running_bin_big   = np.zeros(index.shape)
    running_bin_small = np.zeros(index_shape)

    running_bin_big, running_bin_small = get_bin_timeseries(index=index, percent=percent, window=window)

    fig, ax = plt.subplots(figsize=(20,10))

    ax.plot(time,running_bin_big,color='red',linewidth=3,alpha=2,label=label_large);
    ax.plot(time,running_bin_small,color='blue',linewidth=3,alpha=10,label=label_small);

    #Alter Ticks and labels
    ax.set_xticklabels(ax.get_xticks().astype(int),fontsize=20)
    ax.set_yticklabels(ax.get_yticks().astype(int),fontsize=20)
    ax.set_ylabel('Number of Extreme Events',fontsize=25)

    ax.legend(fontsize=30,loc=legend_loc);
    ax.set_title(title,fontsize=35);

    return fig,ax

"""
This function offers the same functionality of the running_bin_extremes_era data but is altered to account for the CESM_LME
meaning it assumed the first dimension of the data is ensemble member number.

usage: pass index defined along ensemble member and time. observe the default parameters used in the calculation and
       rendering of the plot and change if desired
"""

def running_bin_extremes_lme(index=None, percent=20, window=11, label_large = '',label_small = '',title = '',legend_loc=2, index_has_nan=False,ylim=None):

    #initialize containers for running bin diagnostic
    running_bin_big   = np.zeros(index.shape)
    running_bin_small = np.zeros(index.shape)

    #iterate through ensemble members
    for e_mem in range(0,13):

        running_bin_big[e_mem,:],running_bin_small[e_mem,:] = get_bin_timeseries(index=index[e_mem,:],
                                                                                 percent=percent, window=window)

    fig, ax = plt.subplots(figsize=(20,10))

    ax.plot(time,running_bin_big.mean(axis=0),color='red',linewidth=3,alpha=2,label=label_large);
    ax.plot(time,running_bin_small.mean(axis=0),color='blue',linewidth=3,alpha=10,label=label_small);

    #Alter Ticks and labels
    ax.set_xticklabels(ax.get_xticks().astype(int),fontsize=20)
    ax.set_yticklabels(ax.get_yticks().astype(int),fontsize=20)
    ax.set_ylabel('Number of Extreme Events',fontsize=25)

    ax.legend(fontsize=30,loc=legend_loc);
    ax.set_title(title,fontsize=35);

    return fig,ax

def get_bin_timeseries(index=None, percent=20, window=11):

    running_bin_big   = np.zeros(index.shape)
    running_bin_small = np.zeros(index.shape)

    #get list of extreme years
    _,where_largest  = get_extremes(series=index,percent=percent,top=True)
    _,where_smallest = get_extremes(series=index,percent=percent,bot=True)

    running_bin_big,_  = extremes_per_window(date_range=np.arange(0,index.size), \
                                                      dates_of_extremes=where_largest, \
                                                      window=window);
    running_bin_small,_= extremes_per_window(date_range=np.arange(0,index.size), \
                                                      dates_of_extremes=where_smallest, \
                                                      window=window);

    return running_bin_big, running_bin_small

"""

This monte carlo test for frequency will help decide whether the frequency of an event in one sample is significantly
    different to other samples

USAGE: This function will use random sampling to determine the liklihood of extreme event frequency.

        if using without consideration of autocorrelation:

            specify the iterations (the more the merrier) with "n_iterations";
            identify the sample size- because the sampling draws from all model runs, multiply the timespan by the
                number of ensemble members (13);
            the bool array should be a flattened array describing at each timestep, in each ensemble member whether
                an extreme event occured;

        if considering auto correlation:

            specify the iterations (the more the merrier) with "n_iterations";
            identify the sample size- this is where consideration of autocorrelation differs from the first method-
                here you just identify the timespan that is to be sampled from each ensemble member;
            bool array is also different when considering autocorelation- calculate extremes independantly for each
                ensemble member and pass the bool array as an unflattened 13 by 1156 array.

RETURNS: Given a confidence level, this will return the lower and upper threshold for confidence. The lower threshold
        is the number of extremes in a given sample that would qualify as significantly low. THe uppwer threshold in
        is the number of extremes in a given sample that would qualify as significanly high.

"""

def monte_carlo_frequency(n_iterations, sample_size, bool_array, auto=None):

    #initialize array to hold the number of events occuring in each sample
    n_true = np.empty([n_iterations])

    if auto == False:
        for i in range(0,n_iterations):

            #sampling 'sample_size' items from bool array, we record the number of time the conditon in question
            #is satisfied
            n_true[i] = np.where(np.random.choice(bool_array,sample_size,replace=False))[0].size

        return n_true

    else:
        if bool_array.shape[0] != 13:
            print('not viable bool array. please read the "USAGE" section in the function definiton.')

        for i in range(0,n_iterations):

            n_true_temp = np.empty([13,sample_size])

            for j in range(0,13):

                start = np.random.choice(np.arange(0,bool_array.shape[1]-sample_size),1)[0]

                n_true_temp[j,:] = bool_array[j,start:start+sample_size]

            n_true[i]=np.where(n_true_temp.flatten())[0].size

        return n_true
