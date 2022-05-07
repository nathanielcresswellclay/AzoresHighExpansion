from cartopy import config
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import DataAnalyzing as da

################################################################################
################################################################################
#############    PLOTING FUNCTIONS OVERTOP OF CARTOPY PROJECTIONS    ###########
################################################################################
################################################################################

#data     is a list of arrays which contain the values to be put into the histogram
#labels   is a list of strings to label the histogram
#title    is a string- will be title
#xlabel   is a string- will be label of x axis
def multiple_hist(data,title,xlabel,labels,save_fig=False,fig_name=None):

    fig,ax = plt.subplots()
    fig.set_size_inches(8,5)

    ax.hist(data,histtype='bar',labels=labels,alpha=.5,bins=np.arange(850,2006,10))
    ax.set_xlabel(xlabel,fontsize=20)
    ax.set_title(title, fontsize=17)
    ax.legend(loc='upper left')

    if save_fig:
        plt.savefig(fig_name)

def multiple_boxplot(data,labels,title,save_fig=False,fig_name=None,fig_width=5,fig_height=10,font_size=20,x_label=None):

    fig,ax=plt.subplots()
    fig.set_size_inches(fig_width,fig_height)
    ax.set_title(title,fontsize=20)

    ax.boxplot(data);
    ax.set_xticklabels(labels,fontsize=font_size);
    ax.set_xlabel(x_label,fontsize=20)
    if save_fig:
        plt.savefig(fig_name)

    return

def north_atlantic_contour(fig=None,ax=None,data=None,lat=None,lon=None,title='',levels=None,subplot = False,cmap='bwr'):

    if fig==None:
        fig,ax = plt.subplots()

    if np.all(data==None):
        print("function call----> fig,ax = north_atlantic_contour(data=None,lat=None,lon=None,title='',levels=None)")
        return

    ax=plt.axes(projection=ccrs.Orthographic(central_longitude=-30,central_latitude=45))
    ax.set_global()

    if np.all(levels)==None:
        im = ax.contourf(lon,lat,data,levels=20,transform=ccrs.PlateCarree(),extend='both',cmap=cmap)
    else:
        im = ax.contourf(lon,lat,data,levels=levels,transform=ccrs.PlateCarree(),extend='both',cmap=cmap)


    ax.coastlines()
    ax.gridlines()

    if subplot:
        plt.show()
        return fig,ax,im
    else:
        ax.set_title(title,fontsize=17);
        ax.title.set_position([.5, 1.05])
        fig.colorbar(im,orientation='horizontal',pad=0.05)
        plt.show()
        return fig,ax


def plot_timeseries(fig=None,ax=None,series=None,color='grey',xticks=None,xtick_labels=None,yticks=None,ytick_labels=None,xlabel='',ylabel='',title='Title',extremes_percent=False,extreme_percent=2,extremes_threshold=False,threshold=None):

    if np.all(series==None):
        print("plot_timeseries function call: plot_timeseries(fig=None,ax=None,series=None,color='grey',xticks=None,xtick_labels=None,yticks=None,ytick_labels=None,xlabel='',ylabel='',title='Title',extremes=False)")
        return

    if fig==None:
        fig,ax=plt.subplots()

    ax.plot(series,color=color);

    if np.all(xticks!=None):
        ax.set_xticks(xticks);
    if np.all(xtick_labels!=None):
        ax.set_xticklabels(xtick_labels,fontsize=20);
    if np.all(yticks!=None):
        ax.set_yticks(yticks);
    if np.all(ytick_labels!=None):
        ax.set_yticklabels(ytick_labels,fontsize=20);
    if np.all(ylabel!=''):
        ax.set_ylabel(ylabel,fontsize=22);
    if np.all(xlabel!=''):
        ax.set_xlabel(xlabel,fontsize=22)
    if np.all(title!='Title'):
        ax.set_title(title,fontsize=28);

    if extremes_percent:

        highest,where_highest = get_extremes(series=series, percent=extreme_percent,top=True)
        ax.scatter(where_highest,highest,color='salmon',marker='^')
        lowest, where_lowest  = get_extremes(series=series, percent=extreme_percent,bot=True)
        ax.scatter(where_lowest,lowest,color='lightsteelblue',marker='v');

    elif extremes_threshold:

        highest,where_highest = get_extremes(series=series, percent=extreme_percent,top=True)
        ax.scatter(where_highest,highest,color='salmon',marker='^')
        lowest, where_lowest  = get_extremes(series=series, percent=extreme_percent,bot=True)
        ax.scatter(where_lowest,lowest,color='lightsteelblue',marker='v');

    return fig,ax

#helper to plot timeseries, returns extreme events and there location in time
def get_extremes(series=None, percent=5,top=False,bot=False):

    if np.any(np.isnan(series)):
        has_nan = True
    else:
        has_nan = False

    if has_nan:
        if top:
            return np.sort(series)[:-(np.where(np.isnan(series))[0].size)][int(-series.size*(percent*.01)):], np.argsort(series)[:-(np.where(np.isnan(series))[0].size)][int(-series.size*(percent*.01)):]
        if bot:
            return np.sort(series)[:-(np.where(np.isnan(series))[0].size)][:int(series.size*(percent*.01))], np.argsort(series)[:-(np.where(np.isnan(series))[0].size)][:int(series.size*(percent*.01))]
    else:
        if top:
            return np.sort(series)[int(-series.size*(percent*.01)):], np.argsort(series)[int(-series.size*(percent*.01)):]
        if bot:
            return np.sort(series)[:int(series.size*(percent*.01))], np.argsort(series)[:int(series.size*(percent*.01))]


#make north atlantic map with lat lon mapped
def northatlantic_map():

    fig = plt.figure(figsize=(20,20));

    ax=plt.subplot(1,2,1,projection=ccrs.PlateCarree());
    ax.set_global()
    ax.gridlines()
    ax.gridlines(crs=ccrs.PlateCarree(), linewidth=2, color='grey', alpha=0.5, draw_labels=True)
    ax.coastlines()
    ax.set_extent([-60,30,10,70])
    return fig,ax
