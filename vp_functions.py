# CVP functional code
# Used to extract CVP data from cf-radial files using Pyart
# Adpated by R. Neely March 2022

import numpy as np
import pyart
import h5py as h5
from netCDF4 import num2date
import kdp_functions as kdpfun
import math
import datetime as dt
from scipy import stats
from vp_io import read_file

def field_fill_to_nan(radar, mask_field):

    if mask_field in radar.fields.keys():
        try:
            if(np.ma.is_masked(radar.fields[mask_field]['data'])):
                mask = radar.fields[mask_field]['data'].mask
                radar.fields[mask_field]['data'] = radar.fields[mask_field]['data'].data
                radar.fields[mask_field]['data'][mask==True] = np.nan
        except:
            print('''field_fill_to_nan doesn't work, try again''')
            raise
    else:
        print ("{} is not in radar.fields.keys() {}".format(mask_field,radar.fields.keys()))
    return

# this will now get all the indices for the cvp
def static_index_for_csv_file(radar, el_sort_ix, lat, lon, avg_range_delta, equidistant_alt, equidistant_bound, verbose=True):
    
    lat_data = radar.gate_latitude['data'].reshape((radar.nsweeps,int(radar.nrays/radar.nsweeps),radar.ngates))
    lat_data=lat_data[el_sort_ix,:,:]
    lon_data = radar.gate_longitude['data'].reshape((radar.nsweeps,int(radar.nrays/radar.nsweeps),radar.ngates))
    lon_data=lon_data[el_sort_ix]
    #alt_data = radar.gate_altitude['data'].reshape((radar.nsweeps,int(radar.nrays/radar.nsweeps),radar.ngates))

    #in_lat = (np.abs(lat_data[0,:,:] - lat)).argmin()
    #in_lon = (np.abs(lon_data[0,:,:] - lon)).argmin()
    X = np.sqrt( np.square( lat_data[0,:,:] - lat ) +  np.square( lon_data[0,:,:] - lon ) )
    idx = np.where( X == X.min() )

    r = idx[1][0] # centre radius index
    az = idx[0][0] # centre az index

    # now find the indices in this CVP
    # get the number of elements corresponding to the 20 km ditance.
    steps_in_range_delta = int(avg_range_delta * 1000/(radar.range['data'][1]-radar.range['data'][0]))
    if (verbose):
        print ("int({} * {}) = {}".format(avg_range_delta,1000/(radar.range['data'][1]-radar.range['data'][0]),steps_in_range_delta))
        print ("steps_in_range_delta are {} and it corresponds to {} km".format(steps_in_range_delta,avg_range_delta))

    r_indxs = get_r_indexes(r,steps_in_range_delta,len(radar.range['data']),verbose=verbose)

    # get the size of the avg_az_delta based on steps_in_range_delta and r
    range_dist = (radar.range['data'][1]-radar.range['data'][0]) * r # steps_in_range_delta * r

    if(verbose):
        print ("range_dist for r = {} is {} originally it was {}".format(r,range_dist,avg_range_delta))
        print (" az resolution is {}".format((360.0/(radar.nrays/radar.nsweeps))))
        print ("1/2 beam angle is {}".format((360.0/(radar.nrays/radar.nsweeps))/2.0))
        print ("in radians {}".format(math.radians((360.0/(radar.nrays/radar.nsweeps))/2.0)))

    # Note: in the case of cvp dynamic the avg_az_delta=5 and was not modified by avg_range_delta
    az_size = 2.0 * range_dist * math.tan(math.radians((360.0/(radar.nrays/radar.nsweeps))/2.0))
    if(verbose):print ("azimuthal size for this location is 2.0 * {} / math.ctan({})= 2.0 * {} * {} = {}".format(range_dist,math.radians((360.0/(radar.nrays/radar.nsweeps))/2.0),range_dist, 1/math.tan(math.radians((360.0/(radar.nrays/radar.nsweeps))/2.0)),az_size))
    avg_az_delta = int(avg_range_delta * 1000/az_size)
    if(verbose):print ("int({} * 1000/{}) = avg_az_delta is {}".format(avg_range_delta,az_size,avg_az_delta))

    az_indxs = get_az_indexes(az,avg_az_delta, verbose=verbose)

    # get altitudes from the radar object
    altitudes = radar.fields['scan_altitude']['data']
    # reshape the field to 3D array
    altitudes = altitudes.reshape((radar.nsweeps,int(radar.nrays/radar.nsweeps),radar.ngates))

    # get altitudes for the column values
    column_altitudes = altitudes[np.ix_(range(0,radar.nsweeps),az_indxs,r_indxs)] 
    # resort the elevations
    column_altitudes=column_altitudes[el_sort_ix,:,:]

    # Note: this next part calculating the el_indxs was not done for cvp dynamic

    if verbose:
        print ("!!!! column_altitudes.shape")
        print (column_altitudes.shape)
        print ("np.unique(column_altitudes)")
        print (len(np.unique(column_altitudes)))

    # calculate mean altitudes for each range of the selected column - JCrook not sure what this is for!
    altitudes_to_smooth = np.nanmean(column_altitudes,axis=1)
    if verbose:
        print ("altitudes_to_smooth")
        print (altitudes_to_smooth.shape)
        print ("np.nanmean(altitudes_to_smooth,axis=1)")
        print (np.nanmean(altitudes_to_smooth,axis=1).shape)
    # flattened altitude values array
    flat_altitudes = altitudes_to_smooth.flatten().reshape((np.prod(altitudes_to_smooth.shape),1))
    if verbose: print ("flat_altitudes.shape is {}".format(flat_altitudes.shape))

    # get the el_indx for each height
    el_indxs=get_el_indexes(column_altitudes.flatten(), len(equidistant_alt), equidistant_bound)
    
    index_list =[r_indxs, az_indxs, el_indxs]
    print('r',r_indxs)
    print('az', az_indxs)
    return index_list

# this function converts the indices in a cvp to the voxels in a radar file
# Note that not all radar files have their elevations in the same order but
# the indices have been calculated with the elevations sorted first
def cvp_static_indexes_to_voxels(nrad,naz,nel,r_indx,az_indx,el_indxs):
    polar_data=np.zeros((nel,naz,nrad),int)+np.nan
    this_polar_data=polar_data[np.ix_(np.arange(nel), az_indx, r_indx)]
    this_polar_data=this_polar_data.flatten()
    nheights=len(el_indxs)
    for h in range(nheights):
        this_el_ixs=el_indxs[h]
        if len(this_el_ixs)>0:
            this_polar_data[this_el_ixs]=1
    polar_data[np.ix_(np.arange(nel), az_indx, r_indx)]=this_polar_data.reshape(nel,len(az_indx), len(r_indx))
    voxel_ix=np.where(polar_data==1)
    return voxel_ix
    
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except (ValueError):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')


def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.


    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    box = np.ones(window_len)/window_len
    y = np.convolve(x, box, mode='same')
    y[0] = x[0]# to remove jumps at the lowest elevation
    return y


def altitude_parameter_averaging_cvp(radar, cvp_index_data, field, azimuth_exclude = None, verbose=True):

    if verbose:
        print ("altitude_parameter_averaging_cvp")

    bin_centers = []
    bin_count = []
    bin_means = []
    bin_std = []

    altitudes = []
    observation_count = []
    mean_values = []
    std_values = []
    timeofsweep = []

    # az delta in degrees
    avg_az_delta = 5

    if(verbose):
        print ("altitude_parameter_averaging_cvp is running for the field {}".format(field))

    mask_field = field
    #expected_ele_index = int(cvp_index_data[2])
    # Note: el_indxs are not used here
    [r_indxs,az_indxs,el_indxs] = cvp_index_data

    #elevations = radar.elevation['data'].reshape((radar.nsweeps,int(radar.nrays/radar.nsweeps)))

    data = np.array(radar.fields[field]['data'].reshape((radar.nsweeps,int(radar.nrays/radar.nsweeps),radar.ngates)))

    field_fill_to_nan(radar, mask_field)
    try:data = np.where(data == radar.fields[field]['_FillValue'], np.nan, data)
    except:print ("'_FillValue' is not used in the Field {}".format(field))

    column = data[np.ix_(range(0,radar.nsweeps),az_indxs,r_indxs)]

    # get altitudes from the radar object
    altitudes = radar.fields['scan_altitude']['data']#[radar.sweep_start_ray_index['data'][sweep], :]

    # reshape the field to 3D array
    altitudes = altitudes.reshape((radar.nsweeps,int(radar.nrays/radar.nsweeps),radar.ngates))

    # get altitudes for the column values
    column_altitudes = altitudes[np.ix_(range(0,radar.nsweeps),az_indxs,r_indxs)]

    # calculate mean altitudes for each range of the selected column
    altitudes_to_smooth = np.nanmean(column_altitudes,axis=1)

    # flattened altitude values array
    flat_altitudes = altitudes_to_smooth.flatten().reshape((np.prod(altitudes_to_smooth.shape),1))

    if mask_field in ['dBuZ', 'dBZ', 'dBZ_ac', 'dBuZv', 'dBZv']:
        column = np.power(10, column / 10.0)

    if mask_field in ['uPhiDP']:
        elev = column.shape[0]
        rays = column.shape[1]
        bins = column.shape[2]
        #METEO_THRESH=0.7
        #flags = np.zeros((elev, rays, bins))
        # generate non-meteo mask: - commented out as already done in preprocessing step
        #rhohv_data = radar.fields['RhoHV']['data']
        #rhohv_3D_data = rhohv_data.reshape((radar.nsweeps,int(radar.nrays/radar.nsweeps),radar.ngates))
        #rhohv = rhohv_3D_data[np.ix_(range(0,radar.nsweeps),az_indxs,r_indxs)]
        #kdpfun.generate_meteo_mask(elev, rays, bins, flags, rhohv, METEO_THRESH)
        meteoMask = radar.fields['meteoMask']['data']
        column = kdpfun.unwrap_phidp(elev, rays, bins, meteoMask, column)

    #summed = np.zeros(column.shape)
    counted = np.zeros(column.shape)
    mask = np.ma.masked_invalid(column).mask
    inv_mask = np.where(mask, 0, 1)

    #summed = np.nansum(column, axis=1)

    counted = np.nansum(inv_mask, axis=1)

    observation_count = np.where(counted == 0, np.nan, counted)

    mean_values = np.nanmean(column, axis=1)

    std_values = np.nanstd(column, axis=1)

    if mask_field in ['dBuZ', 'dBZ', 'dBZ_ac', 'dBuZv', 'dBZv']:#mask_field == 'dBuZ' or mask_field == 'dBZ':
        mean_values = 10 * np.log10(mean_values)
        std_values = 10 * np.log10(std_values)
    #observation_count = np.nansum(counted, axis=1)

    # flattened mean array
    flat_means = mean_values.flatten().reshape((np.prod(mean_values.shape),1))

    # flattened std array
    flat_std = std_values.flatten().reshape((np.prod(std_values.shape),1))
    flat_count = observation_count.flatten().reshape((np.prod(observation_count.shape),1))

    output_2D_array = []
    output_2D_array = np.append(flat_altitudes, flat_means, axis=1)
    output_2D_array = np.append(output_2D_array, flat_std, axis=1)
    output_2D_array = np.append(output_2D_array, flat_count, axis=1)
    sorted_by_altitude_output_array = output_2D_array[output_2D_array[:,0].argsort()]

    masked_sorted_by_altitude_output_array = sorted_by_altitude_output_array[~np.isnan(sorted_by_altitude_output_array).any(axis=1)]
    masked_sorted_by_altitude_output_array = masked_sorted_by_altitude_output_array[~np.isinf(masked_sorted_by_altitude_output_array).any(axis=1)]

    if(len(masked_sorted_by_altitude_output_array)>0):
        #x = np.linspace(np.min(masked_sorted_by_altitude_output_array[:,0]),np.max(masked_sorted_by_altitude_output_array[:,0]),10000)
        bin_means, bin_edges, binnumber = stats.binned_statistic(masked_sorted_by_altitude_output_array[:,0],masked_sorted_by_altitude_output_array[:,1], statistic='mean', bins=np.linspace(0,10000,num=101))#bins=np.linspace(0,5000,num=51))
        bin_std, bin_edges, binnumber = stats.binned_statistic(masked_sorted_by_altitude_output_array[:,0],masked_sorted_by_altitude_output_array[:,1], statistic='std', bins=np.linspace(0,10000,num=101))#bins=np.linspace(0,10000,num=101))
        bin_count, bin_edges, binnumber = stats.binned_statistic(masked_sorted_by_altitude_output_array[:,0],masked_sorted_by_altitude_output_array[:,1], statistic='count', bins=np.linspace(0,10000,num=101))#bins=np.linspace(0,10000,num=101))
        bin_width = (bin_edges[1] - bin_edges[0])
        bin_centers = bin_edges[1:] - bin_width/2
        #win=3
        #smoothed_mean = smooth(bin_means,win)

    timeofsweep = num2date(np.nanmean(radar.time['data'][:]),
                                       radar.time['units'],
                                       radar.time['calendar'])

    zero_array = np.zeros((column.shape[1],))
    zero_array[:] = np.nan

    if(len(masked_sorted_by_altitude_output_array)>0):return bin_centers, bin_count, bin_means, bin_std, timeofsweep
    else:return np.nanmean(altitudes_to_smooth,axis=1), np.nanmean(observation_count,axis=1) , np.nanmean(mean_values,axis=1), np.nanmean(std_values,axis=1) , timeofsweep


def get_az_indexes(az,avg_az_delta, verbose = True):
    if (verbose):
        print ("in get_az_indexes(az = {},avg_az_delta = {})".format(az,avg_az_delta))

    if(az<avg_az_delta):
        from_az = range(359 + (az-avg_az_delta),360)
        to_az = range(0,az+avg_az_delta+1)
    elif(az+avg_az_delta>359):
        from_az = range(0,(az+avg_az_delta)-359)
        to_az = range(az-avg_az_delta,360)
    else:
        from_az = range(az-avg_az_delta,az)
        to_az = range(az,az+avg_az_delta+1)
    az_indxs = []
    az_indxs.extend(from_az)
    az_indxs.extend(to_az)
    return(az_indxs)


def get_r_indexes(r, steps_in_range_delta, max_range,verbose = True):

    if (verbose):
        print ("in get_r_indexes(r = {}, steps_in_range_delta = {})".format(r, steps_in_range_delta))

    if(r+steps_in_range_delta>max_range):
        from_r = range(r-steps_in_range_delta-1,r)
        to_r =  range(r,max_range)
    elif(r-steps_in_range_delta<0):
        from_r = range(0,r)
        to_r = range(r,r+steps_in_range_delta)
    else:
        from_r = range(r-steps_in_range_delta,r)
        to_r = range(r,r+steps_in_range_delta+1)

    r_indxs = []
    r_indxs.extend(from_r)
    r_indxs.extend(to_r)
    return(r_indxs)
    
# col_altitudes should be a flattened array of all horizontal points in cvp
def get_el_indexes(col_altitudes, nheights, height_bounds):
    el_indxs=[]
    for h in range(nheights):
        if h == nheights-1:
            this_el_indxs=np.where(col_altitudes>=height_bounds[h])
        else:
            this_el_indxs=np.where((col_altitudes>=height_bounds[h]) & (col_altitudes<height_bounds[h+1]))
        el_indxs.append(this_el_indxs[0])
        
    return el_indxs

def altitude_parameter_averaging_cvp_static(radar, field, cvp_indxs, equidistant_alt,
    azimuth_exclude = None, store_indexes=False, verbose=True):

    [r_indxs,az_indxs,el_indxs] = cvp_indxs

    if verbose:
        print ("altitude_parameter_averaging_cvp_static")

    bin_centers = []
    bin_count = []
    bin_means = []
    bin_std = []
    timeofsweep = []

    altitudes = []
    observation_count = []
    mean_values = []
    std_values = []
    timeofsweep = []

    mask_field = field
    #expected_ele_index = el

    field_fill_to_nan(radar, mask_field)

    # get the elevations in ascending order as was used for static indices
    elevations = radar.elevation['data'].reshape((radar.nsweeps,int(radar.nrays/radar.nsweeps)))
    elevations=np.unique(elevations, axis=1)
    el_sort_ix=np.argsort(elevations)

    data = np.array(radar.fields[field]['data'].reshape((radar.nsweeps,int(radar.nrays/radar.nsweeps),radar.ngates)))


    column = data[np.ix_(range(0,radar.nsweeps),az_indxs,r_indxs)]
    column=column[el_sort_ix,:,:]

    if mask_field in ['dBuZ', 'dBZ', 'dBZ_ac', 'dBuZv', 'dBZv']:column = np.power(10, column / 10.0)

    # unfold uPhiDP values
    # remove phi_dp wrap-around:
    #print 'unwrapping phidp ...'
    if mask_field in ['uPhiDP']:
        elev = column.shape[0]
        rays = column.shape[1]
        bins = column.shape[2]
        METEO_THRESH=0.7
        flags = np.zeros((elev, rays, bins))
        # generate non-meteo mask:
        rhohv_data = radar.fields['RhoHV']['data']
        rhohv_3D_data = rhohv_data.reshape((radar.nsweeps,int(radar.nrays/radar.nsweeps),radar.ngates))
        rhohv = rhohv_3D_data[np.ix_(range(0,radar.nsweeps),az_indxs,r_indxs)]
        (meteoMask) = kdpfun.generate_meteo_mask(elev, rays, bins, flags, rhohv, METEO_THRESH)
        column = kdpfun.unwrap_phidp(elev, rays, bins, meteoMask, column)

    # get mean of means for equidistant altitudes
    nheights=len(equidistant_alt)
    equdist_mean = np.zeros(nheights)*np.nan
    equdist_std = np.zeros(nheights)*np.nan
    equdist_count = np.zeros(nheights)*np.nan
    equdist_total = np.zeros(nheights)*np.nan
    equdist_indexes = [None]*nheights if store_indexes else None

    column=column.flatten()
    for h in range(nheights):
        el_ixs=el_indxs[h]
        if len(el_ixs)>0:
            temp_column = column[el_ixs]
            valid_indexes = np.where(np.isfinite(temp_column))[0]
            if len(valid_indexes)>0:
                equdist_mean[h] = np.mean(temp_column[valid_indexes])
                equdist_std[h] = np.std(temp_column[valid_indexes])
            equdist_count[h] = len(valid_indexes)
            equdist_total[h] = len(temp_column)
            if store_indexes:
                # Store indexes that are included in generating means
                index_bitmap=np.zeros(len(el_ixs),int)
                index_bitmap[valid_indexes]=1
                equdist_indexes[h] = index_bitmap

    if nheights >= 300: win = 5
    else:win = 3

    equdist_mean = smooth(equdist_mean,win)

    if mask_field in ['dBuZ', 'dBZ', 'dBZ_ac', 'dBuZv', 'dBZv']:
        equdist_mean = 10 * np.log10(equdist_mean)
        equdist_std = 10 * np.log10(equdist_std)

    timeofsweep = num2date(np.nanmean(radar.time['data'][:]),
                                       radar.time['units'],
                                       radar.time['calendar'])

    return equdist_count, equdist_mean, equdist_std, timeofsweep, equdist_indexes
    #return bin_centers, bin_count, bin_means, bin_std, timeofsweep


def add_dim(arr, arr_name, verbose=False):

    if verbose: ("Adding dimension to {} field for calculations".format(arr_name))

    if len (arr.shape) != 2:
        print ("Tried to add length-1 dummy dimension to array {}, but was not a 2d array initially.".format(arr_name))
        print ("Shape was {}".format(arr.shape))

    return arr.reshape((1,arr.shape[0],arr.shape[1]))


def altitude_parameter_averaging_qvp(radar, elevation, field, azimuth_exclude, verbose=False, meteo=False, METEO_THRESH=0.7):

    if verbose:
        print ("in altitude_parameter_averaging_qvp")

    mask_field = field
    expected_ele = elevation

    zero_array = np.zeros((radar.fields[mask_field]['data'].data.shape[1],))
    zero_array[:] = np.nan
    timeofsweep = num2date(radar.time['data'][0],
                           radar.time['units'],
                           radar.time['calendar'])

    try:
        sweep = int(np.where(abs(radar.elevation['data'] - expected_ele) <= 0.1)[0][0] / (radar.nrays / radar.nsweeps))
    except:
        if verbose: print ("leaving altitude_parameter_averaging_qvp")
        return (zero_array, zero_array, zero_array, zero_array, timeofsweep)
    else:
        if abs(radar.elevation['data'][radar.sweep_start_ray_index['data'][sweep]] - expected_ele) > 0.1:
            # Check that calculated sweep number matches metadata sweep number
            print("Elevation read was {} but expected elevation was {}".format(radar.elevation['data'][radar.sweep_start_ray_index['data'][sweep]], expected_ele))
            if verbose: print ("leaving altitude_parameter_averaging_qvp")
            return (zero_array, zero_array, zero_array, zero_array, timeofsweep)

        else:

            try:field_fill_to_nan(radar, mask_field)
            except:print('No fill value for {}'.format(mask_field))

            sweep_ind = (radar.sweep_start_ray_index['data'][sweep], radar.sweep_end_ray_index['data'][sweep] + 1)
            azimuth_data = radar.azimuth['data'][sweep_ind[0]:sweep_ind[1]]
            if not meteo:
                sweep_data_2d = radar.fields[mask_field]['data'][sweep_ind[0]:sweep_ind[1]]
                sweep_data = add_dim(sweep_data_2d, mask_field, verbose)
            else:
                this_data_2d = radar.fields[mask_field]['data'][sweep_ind[0]:sweep_ind[1]]
                this_data = add_dim(this_data_2d,mask_field, verbose)
                sweep_data = np.full_like(this_data,np.nan)

            timeofsweep = num2date(np.nanmean(radar.time['data'][sweep_ind[0]:sweep_ind[1]]),
                                   radar.time['units'],
                                   radar.time['calendar'])

            elev = sweep_data.shape[0]
            rays = sweep_data.shape[1]
            bins = sweep_data.shape[2]
            flags = np.zeros((elev,rays, bins))
            rhohv = radar.fields['RhoHV']['data'][sweep_ind[0]:sweep_ind[1]].reshape(elev,rays,bins)
            (meteoMask) = kdpfun.generate_meteo_mask(elev, rays, bins, flags, rhohv, METEO_THRESH)

            if meteo:
                snrH = radar.fields['SNR']['data'][sweep_ind[0]:sweep_ind[1]].reshape(elev,rays,bins)
                snrV = radar.fields['SNRv']['data'][sweep_ind[0]:sweep_ind[1]].reshape(elev,rays,bins)
                snrMask = np.full_like(meteoMask,0)
                snrMask[(snrH > 8) & (snrV > 8)] = 1
                sweep_data[(meteoMask == 1) & (snrMask == 1)] = this_data[(meteoMask == 1) & (snrMask == 1)]

            # unfold uPhiDP values
            # remove phi_dp wrap-around:
            if mask_field in ['uPhiDP']:
                sweep_data = kdpfun.unwrap_phidp(elev, rays, bins, meteoMask, sweep_data)
            elif mask_field in ['ZDR'] and meteo:
                sweep_data[(sweep_data < -1.5) | (sweep_data > 5)] = np.nan
            elif mask_field in ['dBuZ', 'dBZ', 'dBuZv', 'dBZv']:
                if meteo: sweep_data[(sweep_data < -10) | (sweep_data > 60)] = np.nan
                sweep_data = np.power(10, sweep_data / 10.0)

            # Apply azimuth exclude mask
            azimuth_mask = np.tile(np.isin(np.round(azimuth_data,0),azimuth_exclude).reshape(len(azimuth_data),1),(1,sweep_data.shape[2]))
            azimuth_mask = add_dim(azimuth_mask, 'azimuth_mask', verbose)
            sweep_data = np.ma.masked_where(azimuth_mask,sweep_data)

            mask = np.ma.masked_invalid(sweep_data).mask
            inv_mask = np.where(mask, 0, 1)
            summed = np.zeros(sweep_data.shape)
            counted = np.zeros(sweep_data.shape)
            summed = np.nansum([summed, sweep_data], axis=0)
            counted = np.nansum([counted, inv_mask], axis=0)
            summed = np.where(counted == 0, np.nan, summed)
            counted = np.where(counted == 0, np.nan, counted)
            mean_values = np.nanmean((summed / counted), axis=1)[0]
            std_values = np.nanstd((summed / counted), axis=1)[0]
            if mask_field in ['dBuZ', 'dBZ', 'dBuZv', 'dBZv']:
                mean_values = 10 * np.log10(mean_values)
                std_values = 10 * np.log10(std_values)
            observation_count = np.nansum(counted, axis=1)[0]
            if expected_ele == 90.0:
                altitudes = radar.range['data']
            else:
                altitudes = radar.fields['scan_altitude']['data'][radar.sweep_start_ray_index['data'][sweep], :]

            # clean up the nearest beans influenced by sidelobs
            ranges=radar.range['data']
            mean_values[ranges<400] = np.nan

            if verbose: print ("leaving altitude_parameter_averaging_qvp")

            return (altitudes, observation_count, mean_values, std_values, timeofsweep)


def time_height(list_of_files, field_list, cvp_indexes = None, equidistant_alt=None,
                elevation = None, count_threshold=0,  met_office=False,
                azimuth_exclude = [], 
                store_indexes=False,
                verbose=False, vp_mode='qvp'):

    import itertools

    if (verbose):
        print ("time_height \n function input is:{} \nwith the number of files {} \n elevation is {} \n field list is {} \n count threshold is {}".format(list_of_files,len(list_of_files),elevation,field_list,count_threshold))

    result_dict = {}
    stddev_dict = {}
    count_dict = {}
    indexes_dict = {}
    sweep_times = []
    empty_array = {}
    unit_dict = []
    long_names = []
    short_names = []
    defaulttime = [dt.datetime(1970,1,1,0,0,0)]

    if(met_office):
        testfile = h5.File(list_of_files[0], 'r')
        times = list(testfile['lp'].keys())
        list_of_files = list(itertools.chain.from_iterable(itertools.repeat(x, len(times)) for x in list_of_files))
        testfile.close()
    else:
        times = [defaulttime]


    for f, file_ in enumerate(list_of_files):
        if (verbose):
            print ("file f {} is {}".format(f,list_of_files[f]))
            print ("for loop file is {}".format(file_))


        # Edited function call so that different call not needed for metoffice=True
        (radar, unit_dict, long_names, short_names) = \
            read_file(f, times[f%len(times)], list_of_files[f], field_list, unit_dict,
                      long_names, short_names, met_office=met_office, verbose=verbose)

        sweep_time = []

        for field in field_list:
            indexes = [] if store_indexes else None
            try:
                if vp_mode == 'qvp':
                    (alts, counts, means, standard_deviations, timeofsweep) = \
                        altitude_parameter_averaging_qvp(radar, elevation, field,
                            azimuth_exclude = azimuth_exclude, verbose=verbose,
                            meteo=False, METEO_THRESH=0.7)
                elif vp_mode == 'cvp_static':
                    (counts, means, standard_deviations, timeofsweep, indexes) = \
                        altitude_parameter_averaging_cvp_static(radar, field,
                            cvp_indexes, equidistant_alt,
                            azimuth_exclude = azimuth_exclude, store_indexes=store_indexes, verbose=verbose)
                    alts=equidistant_alt
                elif vp_mode == 'cvp_dynamic':
                    (alts, counts, means, standard_deviations, timeofsweep) = \
                        altitude_parameter_averaging_cvp(radar, field,
                            cvp_indexes, azimuth_exclude = azimuth_exclude,
                            verbose=verbose)
                else:
                    print("vp_mode unrecognised")
                    raise
            except:
                if f == 0:
                    print("Starting file doesn't work, try again")
                    raise
                else:
                    print('Block1 Error running altitude parameter averaging for field {1} in file {0}'.format(file_, field))
                    empty_array = result_dict[field]
                    stddev_array = stddev_dict[field]
                    counts_array = count_dict[field]
                    indexes_list = indexes_dict[field]
                    counts = np.zeros(counts_array.shape[0])
                    means = np.zeros(counts_array.shape[0])
                    standard_deviations = np.zeros(counts_array.shape[0])
                    counts[:] = np.nan
                    means[:] = np.nan
                    standard_deviations[:] = np.nan
                    indexes = [] if store_indexes else None

                    empty_array[:, f] = np.where(counts > count_threshold, means, np.nan)
                    stddev_array[:, f] = np.where(counts > count_threshold, standard_deviations, np.nan)
                    counts_array[:, f] = counts
                    if store_indexes:
                        indexes_list[f] = indexes
                    result_dict.update({field: empty_array})
                    stddev_dict.update({field: stddev_array})
                    count_dict.update({field: counts_array})
                    indexes_dict.update({field: indexes_list})
            else:
                if f == 0:
                    empty_array = np.zeros((means.shape[0], len(list_of_files)))
                    empty_array[:] = np.nan
                    stddev_array = np.zeros((means.shape[0], len(list_of_files)))
                    counts_array = np.zeros((means.shape[0], len(list_of_files)))
                    stddev_array[:] = np.nan
                    result_dict.update({'alts': alts})
                    indexes_list = [None] * len(list_of_files) if store_indexes else None
                else:
                    empty_array = result_dict[field]
                    stddev_array = stddev_dict[field]
                    counts_array = count_dict[field]
                    indexes_list = indexes_dict[field]

                empty_array[:, f] = np.where(counts > count_threshold, means, np.nan)
                stddev_array[:, f] = np.where(counts > count_threshold, standard_deviations, np.nan)
                counts_array[:, f] = counts
                if store_indexes:
                    indexes_list[f] = indexes
                result_dict.update({field: empty_array})
                stddev_dict.update({field: stddev_array})
                count_dict.update({field: counts_array})
                indexes_dict.update({field: indexes_list})
                sweep_time.append(timeofsweep)

        try:
            sweep_times.append(max(sweep_time or defaulttime))
        except ValueError:
            print ("sweep time could not be added for file: {}".format(file_))
            print ("len(sweep_time) = {}".format(len(sweep_time)))
            print ("Skipping this record")
            raise

        radar = None

    return result_dict, stddev_dict, count_dict, sweep_times, unit_dict, long_names, short_names, indexes_dict
