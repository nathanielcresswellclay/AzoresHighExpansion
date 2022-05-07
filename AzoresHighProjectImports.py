import warnings
warnings.filterwarnings("ignore")
import numpy as np
import scipy as sp
import netCDF4 as nc
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.io import netcdf
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import ttest_ind
from scipy.stats import kruskal
from scipy.stats import mannwhitneyu
import pandas as pd
from cartopy import config
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import seaborn as sns
from cartopy.util import add_cyclic_point
import xarray as xr
from IPython.display import clear_output
from IPython.display import display, HTML
from IPython.display import Image
import matplotlib.patches as patches
from scipy.stats import gaussian_kde
import matplotlib.ticker as mtick
import matplotlib.patches as mpatches

display(HTML("""
<style>
.output {
    display: flex;
    align-items: center;
    text-align: center;
}
</style>
"""))
