"""
 Routines for reading NIMROD ODIM_H5 Aggregated files.

"""

import datetime
from netCDF4 import Dataset, date2num
import numpy as np

try:
    import h5py
    _H5PY_AVAILABLE = True
except ImportError:
    _H5PY_AVAILABLE = False

from pyart.config import FileMetadata, get_fillvalue
from pyart.io import read_sigmet, read
from pyart.io.common import make_time_unit_str, _test_arguments
from pyart.core.radar import Radar
from pyart.exceptions import MissingOptionalDependency

from preprocessing import preprocessing,shift_ppi


ODIM_H5_FIELD_NAMES = {
    'TH': 'total_power',        # uncorrected reflectivity, horizontal
    'TV': 'total_power',        # uncorrected reflectivity, vertical
    'DBZH': 'reflectivity',     # corrected reflectivity, horizontal
    'DBZV': 'reflectivity',     # corrected reflectivity, vertical
    'ZDR': 'differential_reflectivity',     # differential reflectivity
    'RHOHV': 'cross_correlation_ratio',
    'LDR': 'linear_polarization_ratio',
    'PHIDP': 'differential_phase',
    'KDP': 'specific_differential_phase',
    'SQI': 'normalized_coherent_power',
    'SNR': 'signal_to_noise_ratio',
    'VRAD': 'velocity', # radial velocity, marked for deprecation in ODIM HDF5 2.2
    'VRADH': 'velocity', # radial velocity, horizontal polarisation
    'VRADV': 'velocity', # radial velocity, vertical polarisation
    'WRAD': 'spectrum_width',
    'QIND': 'quality_index',
    'CI': 'clutter_index', #Add for NIMROD files
    'LONG_RANGE_NOISE_DBC_H': 'LONG_RANGE_NOISE_DBC_H', #Add for NIMROD files
    'LONG_RANGE_NOISE_DBC_V': 'LONG_RANGE_NOISE_DBC_V', #Add for NIMROD files
}


def read_nimrod_aggregated_odim_h5(filename,data_type, time='0000', field_names=None, additional_metadata=None,
                 file_field_names=False, exclude_fields=None,
                 include_fields=None, **kwargs):
    """
    Read a NIMROD ODIM_H5 Aggregated file.
    Editted from pyart odim read function to subset the file contents by
    pulse type and time

    Parameters
    ----------
    filename : str
        Name of the ODIM_H5 file to read.
    data_type : str
        lp (long pulse) or sp (short pulse) data to be read  #Need to expand to read vertical
    time : str
        string with the format hhmm that split the day into 144 ten minute chunks for sp data
        or 288 5 minute chunks for lp data from 0000 to 2350.
        Each period has been aggregated to make single volume.
    field_names : dict, optional
        Dictionary mapping ODIM_H5 field names to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        Py-ART configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata from during this read.
        This metadata is not used during any successive file reads unless
        explicitly included.  A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the Py-ART configuration file will be used.
    file_field_names : bool, optional
        True to use the MDV data type names for the field names. If this
        case the field_names parameter is ignored. The field dictionary will
        likely only have a 'data' key, unless the fields are defined in
        `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields specified by include_fields.
    include_fields : list or None, optional
        List of fields to include from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields not specified by exclude_fields.


    Returns
    -------
    radar : Radar
        Radar object containing data from ODIM_H5 file.

    """

    # TODO before moving to pyart.io
    # * unit test
    # * add default field mapping, etc to default config
    # * auto-detect file type with pyart.io.read function
    # * instrument parameters
    # * add additional checks for HOW attributes
    # * support for other objects (SCAN, XSEC)

    # check that h5py is available
    if not _H5PY_AVAILABLE:
        raise MissingOptionalDependency(
            "h5py is required to use read_odim_h5 but is not installed")

    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    if field_names is None:
        field_names = ODIM_H5_FIELD_NAMES
    filemetadata = FileMetadata('odim_h5', field_names, additional_metadata,
                                file_field_names, exclude_fields,
                                include_fields)

    # open the file
    with h5py.File(filename, 'r') as hfile:
        try:
            hfile=hfile[data_type][time]
        except:
            output="/gws/smf/j04/ncas_radar/rrniii/BioDAR/CVP_Extraction/Output"
            with open(output+'CVP_Cols_Done.txt', 'a') as log:
                        log.write("No Data: "+time+'for '+data_type+' in '+filename+ '\n')


        odim_object = _to_str(hfile['what'].attrs['object'])
        if odim_object not in ['PVOL', 'SCAN', 'ELEV', 'AZIM']:
            raise NotImplementedError(
                'object: %s not implemented.' % (odim_object))

        # determine the number of sweeps by the number of groups which
        # begin with dataset
        datasets = [k for k in hfile if k.startswith('dataset')]
        datasets.sort(key=lambda x: int(x[7:]))

        #Added to remove vert scans from NIMROD aggregated data as they are a different resolution
        for d in datasets:
            if hfile[d]['where'].attrs['elangle']>89.0:
                datasets.remove(d)


        nsweeps = len(datasets)

        # latitude, longitude and altitude
        latitude = filemetadata('latitude')
        longitude = filemetadata('longitude')
        altitude = filemetadata('altitude')

        h_where = hfile['where'].attrs
        latitude['data'] = np.array([h_where['lat']], dtype='float64')
        longitude['data'] = np.array([h_where['lon']], dtype='float64')
        altitude['data'] = np.array([h_where['height']], dtype='float64')

        # metadata
        metadata = filemetadata('metadata')
        metadata['instrument_name']=_get_radar_name_from_radar_number(hfile['what'].attrs['source_local_site_number'])
        metadata['source'] = _to_str(hfile['what'].attrs['source'])
        metadata['original_container'] = 'ukmo_nimrod_aggregated_odim_h5'
        metadata['odim_conventions'] = _to_str(hfile.attrs['Conventions'])

        h_what = hfile['what'].attrs
        metadata['version'] = 'nimrod_test'#_to_str(h_what['version']) ####Need to talk to josh about it
        metadata['source'] = _to_str(h_what['source'])

        try:
            ds1_how = hfile[datasets[0]]['how'].attrs
        except KeyError:
            # if no how group exists mock it with an empty dictionary
            ds1_how = {}
        if 'system' in ds1_how:
            metadata['system'] = ds1_how['system']
        if 'software' in ds1_how:
            metadata['software'] = ds1_how['software']
        if 'sw_version' in ds1_how:
            metadata['sw_version'] = ds1_how['sw_version']

        # sweep_start_ray_index, sweep_end_ray_index
        sweep_start_ray_index = filemetadata('sweep_start_ray_index')
        sweep_end_ray_index = filemetadata('sweep_end_ray_index')

        if odim_object in ['AZIM', 'SCAN', 'PVOL']:
            rays_per_sweep = [
                int(hfile[d]['where'].attrs['nrays']) for d in datasets]
        elif odim_object == 'ELEV':
            rays_per_sweep = [
                int(hfile[d]['where'].attrs['angles'].size) for d in datasets]
        total_rays = sum(rays_per_sweep)
        ssri = np.cumsum(np.append([0], rays_per_sweep[:-1])).astype('int32')
        seri = np.cumsum(rays_per_sweep).astype('int32') - 1
        sweep_start_ray_index['data'] = ssri
        sweep_end_ray_index['data'] = seri

        # sweep_number
        sweep_number = filemetadata('sweep_number')
        sweep_number['data'] = np.arange(nsweeps, dtype='int32')

        # sweep_mode
        sweep_mode = filemetadata('sweep_mode')
        sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])

        # scan_type
        if odim_object == 'ELEV':
            scan_type = 'rhi'
        else:
            scan_type = 'ppi'

        # fixed_angle
        fixed_angle = filemetadata('fixed_angle')
        if odim_object == 'ELEV':
            sweep_el = [hfile[d]['where'].attrs['az_angle'] for d in datasets]
        else:
            sweep_el = [hfile[d]['where'].attrs['elangle'] for d in datasets]
        fixed_angle['data'] = np.array(sweep_el, dtype='float32')

        # elevation
        elevation = filemetadata('elevation')
        if 'elangles' in ds1_how:
            edata = np.empty(total_rays, dtype='float32')
            for d, start, stop in zip(datasets, ssri, seri):
                edata[start:stop+1] = hfile[d]['how'].attrs['elangles'][:]
            elevation['data'] = edata
        elif odim_object == 'ELEV':
            edata = np.empty(total_rays, dtype='float32')
            for d, start, stop in zip(datasets, ssri, seri):
                edata[start:stop+1] = hfile[d]['where'].attrs['angles'][:]
            elevation['data'] = edata
        else:
            elevation['data'] = np.repeat(sweep_el, rays_per_sweep)

        # range
        _range = filemetadata('range')
        if 'rstart' in hfile['dataset1/where'].attrs:
            # derive range from rstart and rscale attributes if available

            # check that the gate spacing is constant between sweeps
            rstart = [hfile[d]['where'].attrs['rstart'] for d in datasets]
            if any(rstart != rstart[0]):
                raise ValueError('range start changes between sweeps')
            rscale = [hfile[d]['where'].attrs['rscale'] for d in datasets]
            if any(rscale != rscale[0]):
                raise ValueError('range scale changes between sweeps')
            all_sweeps_nbins = [hfile[d]['where'].attrs['nbins'] for d in datasets]
            # check for max range off all sweeps
            max_nbins = max(all_sweeps_nbins)

            if isinstance(max_nbins, np.ndarray):
                max_nbins = max_nbins[0]
            else:
                max_nbins = max(all_sweeps_nbins)

            rscenter = 1e3 * rstart[0] + rscale[0] / 2
            _range['data'] = np.arange(rscenter,
                                       rscenter + max_nbins * rscale[0],
                                       rscale[0], dtype='float32')
            _range['meters_to_center_of_first_gate'] = rstart[0] * 1000.
            _range['meters_between_gates'] = float(rscale[0])
        else:
            # if not defined use range attribute which defines the maximum range
            # in km. There is no information on the starting location of the
            # range bins so we assume this to be 0.
            # This most often occurs in RHI files, which technically do not meet
            # the ODIM 2.2 specs. Section 7.4 requires that these files include
            # the where/rstart, where/rscale and where/nbins attributes.
            max_range = [hfile[d]['where'].attrs['range'] for d in datasets]
            if any(max_range != max_range[0]):
                raise ValueError('maximum range changes between sweeps')
            # nbins is required
            nbins = hfile['dataset1/data1/data'].shape[1]
            _range['data'] = np.linspace(
                0, max_range[0] * 1000., nbins).astype('float32')
            _range['meters_to_center_of_first_gate'] = 0
            _range['meters_between_gates'] = max_range[0] * 1000. / nbins

        # azimuth
        azimuth = filemetadata('azimuth')
        az_data = np.ones((total_rays, ), dtype='float32')
        for dset, start, stop in zip(datasets, ssri, seri):
            if odim_object == 'ELEV':
                # all azimuth angles are the sweep azimuth angle
                sweep_az = hfile[dset]['where'].attrs['az_angle']
            elif odim_object == 'AZIM':
                # Sector azimuths are specified in the startaz and stopaz
                # attribute of dataset/where.
                # Assume that the azimuth angles do not pass through 0/360 deg.
                startaz = hfile[dset]['where'].attrs['startaz']
                stopaz = hfile[dset]['where'].attrs['stopaz']
                nrays = stop - start + 1
                sweep_az = np.linspace(startaz, stopaz, nrays, endpoint=True)
            elif ('startazA' in ds1_how) and ('stopazA' in ds1_how):
                # average between start and stop azimuth angles
                startaz = hfile[dset]['how'].attrs['startazA']
                stopaz = hfile[dset]['how'].attrs['stopazA']
                sweep_az = np.angle(
                    (np.exp(1.j*np.deg2rad(startaz)) +
                    np.exp(1.j*np.deg2rad(stopaz))) / 2., deg=True)
            else:
                # according to section 5.1 the first ray points to the
                # northernmost direction and proceeds clockwise for a complete
                # 360 rotation.
                try:
                    astart = hfile[dset]['how'].attrs['astart']
                except KeyError:
                    astart = 0.0
                nrays = hfile[dset]['where'].attrs['nrays']
                da = 360.0 / nrays
                sweep_az = np.arange(astart + da / 2., 360., da, dtype='float32')
            az_data[start:stop+1] = sweep_az
        azimuth['data'] = az_data

        # time
        _time = filemetadata('time')
        if ('startazT' in ds1_how) and ('stopazT' in ds1_how):
            # average between startazT and stopazT
            t_data = np.empty((total_rays, ), dtype='float32')
            for dset, start, stop in zip(datasets, ssri, seri):
                t_start = hfile[dset]['how'].attrs['startazT']
                t_stop = hfile[dset]['how'].attrs['stopazT']
                t_data[start:stop+1] = (t_start + t_stop) / 2
            start_epoch = t_data.min()
            start_time = datetime.datetime.utcfromtimestamp(start_epoch)
            _time['units'] = make_time_unit_str(start_time)
            _time['data'] = t_data - start_epoch
        else:
            t_data = np.empty((total_rays, ), dtype='int32')
            # interpolate between each sweep starting and ending time
            for dset, start, stop in zip(datasets, ssri, seri):
                dset_what = hfile[dset]['what'].attrs
                start_str = _to_str(
                    dset_what['startdate'] + dset_what['starttime'])
                end_str = _to_str(dset_what['enddate'] + dset_what['endtime'])
                start_dt = datetime.datetime.strptime(start_str, '%Y%m%d%H%M%S')
                end_dt = datetime.datetime.strptime(end_str, '%Y%m%d%H%M%S')

                time_delta = end_dt - start_dt
                delta_seconds = time_delta.seconds + time_delta.days * 3600 * 24
                rays = stop - start + 1
                sweep_start_epoch = (
                    start_dt - datetime.datetime(1970, 1, 1)).total_seconds()
                t_data[start:stop+1] = (sweep_start_epoch +
                                        np.linspace(0, delta_seconds, rays))
            start_epoch = t_data.min()
            start_time = datetime.datetime.utcfromtimestamp(start_epoch)
            _time['units'] = make_time_unit_str(start_time)
            _time['data'] = (t_data - start_epoch).astype('float32')

        # fields
        fields = {}
        h_field_keys = [k for k in hfile['dataset1'] if k.startswith('data')]
        odim_fields = [hfile['dataset1'][d]['what'].attrs['quantity'] for d in h_field_keys]
        for odim_field, h_field_key in zip(odim_fields, h_field_keys):
            field_name = filemetadata.get_field_name(_to_str(odim_field))
            if field_name is None:
                continue
            fdata = np.ma.zeros((total_rays, max_nbins), dtype='float32')
            start = 0
            # loop on the sweeps, copy data into correct location in data array
            for dset, rays_in_sweep in zip(datasets, rays_per_sweep):
                try:
                    sweep_data = _get_odim_h5_sweep_data(hfile[dset][h_field_key])
                except KeyError:
                    sweep_data = np.zeros((rays_in_sweep, max_nbins)) + np.NaN
                sweep_nbins = sweep_data.shape[1]
                fdata[start:start + rays_in_sweep, :sweep_nbins] = sweep_data[:]
                # set data to NaN if its beyond the range of this sweep
                fdata[start:start + rays_in_sweep, sweep_nbins:max_nbins] = np.nan
                start += rays_in_sweep
            # create field dictionary
            field_dic = filemetadata(field_name)
            field_dic['data'] = fdata
            field_dic['_FillValue'] = get_fillvalue()
            fields[field_name] = field_dic

    # instrument_parameters
    instrument_parameters = None

    return Radar(
        _time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)


def _to_str(text):
    """ Convert bytes to str if necessary. """
    if hasattr(text, 'decode'):
        return text.decode('utf-8')
    else:
        return text


def _get_odim_h5_sweep_data(group):
    """ Get ODIM_H5 sweet data from an HDF5 group. """

    # mask raw data
    what = group['what']
    raw_data = group['data'][:]

    if 'nodata' in what.attrs:
        nodata = what.attrs.get('nodata')
        data = np.ma.masked_equal(raw_data, nodata)
    else:
        data = np.ma.masked_array(raw_data)
    if 'undetect' in what.attrs:
        undetect = what.attrs.get('undetect')
        data[data == undetect] = np.ma.masked

    offset = 0.0
    gain = 1.0
    if 'offset' in what.attrs:
        offset = what.attrs.get('offset')
    if 'gain' in what.attrs:
        gain = what.attrs.get('gain')
    return data * gain + offset


def _get_radar_name_from_radar_number(radar_number):
    if radar_number==7:
        radar_name='Castor Bay'
    if radar_number==5:
        radar_name='Chenies'
    if radar_number==3:
        radar_name='Clee Hill'
    if radar_number==16:
        radar_name='Cobbacomebe Cross'
    if radar_number==10:
        radar_name='Crug-y-gorrllwyn'
    if radar_number==21:
        radar_name='Dean Hill'
    if radar_number==15:
        radar_name='Druim a\'Starraig'
    if radar_number==14:
        radar_name='Dudwick'
    if radar_number==4:
        radar_name='Hameldon Hill'
    if radar_number==23:
        radar_name='High Moorsley'
    if radar_number==18:
        radar_name='Holehead'
    if radar_number==9:
        radar_name='Ingham'
    if radar_number==12:
        radar_name='Jersey'
    if radar_number==19:
        radar_name='Munduff Hill'
    if radar_number==8:
        radar_name='Predannack'
    if radar_number==20:
        radar_name='Thurnham'
    if radar_number==11:
        radar_name='Wardon-hill'

    return radar_name


def named_fields(radar):

    field_dict = {'differential_phase':'PhiDP',
                  'cross_correlation_ratio':'RhoHV',
                  'normalized_coherent_power':'SQI',
                  'spectrum_width':'W',
                  'reflectivity':'dBZ',
                  'differential_reflectivity':'ZDR',
                  'specific_differential_phase':'KDP',
                  'velocity':'V'}


    for k,v in field_dict.items():
        if k in radar.fields.keys():
            radar.fields[v]=radar.fields.pop(k)

    return radar


def return_units(radar, fields):

    unit_dictionary = {}

    for field in fields:
        if field in radar.fields.keys():
            unit_dictionary.update({field:radar.fields[field]['units']})
        else:
            print ('{} not in original files, need to manually add units'.format(field))

    return unit_dictionary


def return_names(radar, fields):

    long_dictionary = {}
    short_dictionary = {}

    for field in fields:

        if field in radar.fields.keys():

            long_dictionary.update({field:radar.fields[field]['long_name']})


            try:short_dictionary.update({field: radar.fields[field]['standard_name']})
            except:
               try:short_dictionary.update({field: radar.fields[field]['proposed_standard_name']})
               except:short_dictionary.update({field: radar.fields[field]['long_name']})
            else:
               print ('{} not in original files, need to manually add names'.format(field))

    return long_dictionary, short_dictionary


def read_file(f, time, file_, fields,unit_dict=[],long_names=[], short_names=[], vp_mode='qvp', met_office=False, verbose=False):

    if ".RAW" in file_:
        radar = read_sigmet(file_)
        radar = named_fields(radar)
    elif ".h5" in file_:
        data_type='lp'
        radar = read_nimrod_aggregated_odim_h5(file_,data_type,time)
        radar = named_fields(radar)
    else:
        radar = read(file_)

    if f == 0:
        unit_dict = return_units(radar, fields)
        long_names, short_names = return_names(radar, fields)

    try:
        if (verbose): print ("preprocessing is running")
        preprocessing(radar, vp_mode)
        if not vp_mode=='qvp': shift_ppi(radar,fields)
    except:
        print('Preprocessing failed for', file_)

    return (radar, unit_dict, long_names, short_names)


# this writes the indices of the radar data that make up this CVP. 
# Within the CVP files there will be a mask of 1's and 0's for each of the el_indxs at each time
# for each variable to show which of these have valid data
# Inputs:
#     output_file is where we will store this data
#     radar_fname is the radar file we used to create the indices - we save this as a global attribute
#     poss_indxs is [r_indx, az_indx, el_indxs]
#     equidistant_alt is an array of the heights within a CVP
def write_static_cvp_indices_netcdf(output_file, radar_fname, poss_indxs, equidistant_alt):

    output = Dataset(output_file, 'w', format='NETCDF4')
    output.source='Created from '+radar_fname
    r_indx=poss_indxs[0]
    az_indx=poss_indxs[1]
    el_indxs=poss_indxs[2]
    nind=len(r_indx)
    ind_dim = output.createDimension('r_indx_dim', nind)
    nind=len(az_indx)
    ind_dim = output.createDimension('az_indx_dim', nind)
    poss_ind_r = output.createVariable('r_indx', np.int32, ('r_indx_dim'))
    poss_ind_az = output.createVariable('az_indx', np.int32, ('az_indx_dim'))
    poss_ind_r[:]=r_indx
    poss_ind_az[:]=az_indx
    nheights=len(equidistant_alt) 
    height_dim = output.createDimension('Height', nheights)
    height_var = output.createVariable('Height', np.float32, ('Height',))
    height_var[:]=equidistant_alt
    height_var.units = "metres above sea level"
    # as the el_indxs have different lengths we will write each height separately
    for h in range(nheights): 
        this_nind=len(el_indxs[h])  
        if this_nind>0:
            # netcdf doesn't like '/' in dimension names so el dimensions will have names el_indx_height_<h>_len
            # each el_indx variable will be in el_indx/height_<h>_el
            height_name='height_{h:d}'.format(h=h)
            basename='el_indx'
            dimname=basename+'_'+height_name+'_len'
            ind_dim = output.createDimension(dimname, this_nind) 
            ind_el = output.createVariable(basename+'/'+height_name+'_el', np.int32, (dimname))
            ind_el[:]=el_indxs[h]  

    
    output.close()
    
# code to read static cvp indices   
class Static_cvp_ind():
    def __init__(self, static_cvp_data):
        self.heights=[]
        self.height_units=''
        self.r_indx=[]
        self.az_indx=[]
        for var in static_cvp_data.variables:
            if var=='Height':
                heights_coord=static_cvp_data.variables[var]
                self.heights=heights_coord[:]
                self.height_units=heights_coord.units
            elif var=='r_indx':
                self.r_indx=static_cvp_data.variables[var][:]
            elif var=='az_indx':
                self.az_indx=static_cvp_data.variables[var][:]
            else:
                print('unexpected var', var)

        self.el_indxs=[[]]*len(self.heights)
        for g in static_cvp_data.groups:
            if g=='el_indx':
                for v in static_cvp_data[g].variables:
                    # get the height index from the name
                    wsplit=v.split('_')
                    if wsplit[0]=='height' and wsplit[2]=='el':
                        h=int(wsplit[1])
                        self.el_indxs[h]=static_cvp_data[g].variables[v][:]
                    else:        
                        print(wsplit)
       
    def print(self):
        print(len(self.heights),'heights:',self.heights, self.height_units)
        print('r_indx', self.r_indx)
        print('az_indx',self.az_indx)
        print('el_indxs',self.el_indxs)
        
    def get_indxs(self):
        return [self.r_indx, self.az_indx, self.el_indxs]
        
def read_static_cvp_indices(cvp_pathname):
    
    static_cvp_data=Dataset(cvp_pathname, "r", format="NETCDF4")
    this_static_cvp=Static_cvp_ind(static_cvp_data)
    return this_static_cvp
    
def output_netcdf(data_list, output_file, profile_type, field_list, elevation, azimuth_exclude, n_time, verbose=False):

    [result_dict,
     stddev_dict,
     count_dict,
     sweep_times,
     unit_dictionary,
     long_names,
     short_names,
     indexes_dict] = data_list

    # dataset of the output
    output = Dataset(output_file, 'w', format='NETCDF4')

    # Create dimensions
    output.createDimension('Height', result_dict['alts'].shape[0])
    output.createDimension('Time', len(sweep_times))
    # Create coordinate variables for 2 dimensions
    times = output.createVariable('Time', np.float64, ('Time',))
    heights = output.createVariable('Height', np.float32, ('Height',))

    time_units = 'seconds since %i-01-01T00:00:00Z' % sweep_times[0].year
    number_times = [date2num(sweeptime, time_units, calendar='gregorian') for sweeptime in sweep_times]

    ## Add values to the dimension variables
    heights[:] = result_dict['alts']
    Nh=len(result_dict['alts'])
    times[:] = np.array(number_times)
    Nt=len(number_times)
    times.units = time_units
    times.calendar = "gregorian"
    heights.units = "metres above sea level"

    # create dimensions for el_indices - all fields should have the same length of el_indices so don't need to define them in each field
    field=field_list[0]
    if indexes_dict[field] is not None:
        el_ind_for_height_dimnames=[None]*Nh 
        el_ind_for_height_len=np.zeros(Nh, int) 
        for h in range(Nh): 
            if isinstance(indexes_dict[field][0][h], np.ndarray): 
                this_len=len(indexes_dict[field][0][h]) 
                height_name='height_{h:d}'.format(h=h)
                basename='el_indx'
                dimname=basename+'_'+height_name+'_len'
                ind_dim = output.createDimension(dimname, this_len) 
                el_ind_for_height_dimnames[h]=dimname
                el_ind_for_height_len[h]=this_len
 
    if profile_type == "QVP":
        name_str = 'Quasi-vertical {} {} at an elevation angle of %.1f degrees' % elevation
    elif profile_type == "CVP":
        name_str = 'Column-vertical {} {}'

    # Create the field variables
    # for each variable in the list we create several fields
    for field in field_list:
        group_data_construct_means = '/'+field+'/Means'
        group_data_construct_deviations = '/' + field + '/StdDevs'
        group_data_construct_counts = '/' + field + '/Counts'
        group_data_construct_indexes= '/' + field + '/VP_Indexes'
        temp_means = output.createVariable(group_data_construct_means, np.float32, ('Height','Time'))
        temp_stds = output.createVariable(group_data_construct_deviations, np.float32, ('Height', 'Time'))
        temp_counts = output.createVariable(group_data_construct_counts, np.float32, ('Height', 'Time'))
        temp_means[:] = result_dict[field]
        temp_stds[:] = stddev_dict[field]
        temp_counts[:] = count_dict[field]
        temp_means.long_name = name_str.format('mean', field)
        stdev_name = name_str.format('standard deviation of', field)
        count_name = name_str.format('number of observations of', field)
        temp_stds.long_name = stdev_name
        temp_counts.long_name = count_name
        if field in unit_dictionary.keys():
            output[field].units = unit_dictionary[field]
            output[field].long_name = long_names[field]
            output[field].standard_name = short_names[field]
        else:
            output[field].units = 'Default'
            output[field].long_name = 'Default'
            output[field].standard_name = 'Default'
        if indexes_dict[field] is not None:
            # Create variables to store indexes used to calculate vertical profile for this field
            # We have a bitmap of valid indices with a different length for each height
            for h in range(Nh):
                # all times for this height should have same length indexes
                if isinstance(indexes_dict[field][0][h], np.ndarray):
                    height_name='height_{h:d}'.format(h=h)
                    varname=group_data_construct_indexes+'/'+height_name+'_el'
                    vp_ind = output.createVariable(varname, np.byte, (el_ind_for_height_dimnames[h],'Time'))
                    these_indices=np.zeros((el_ind_for_height_len[h],Nt), int)
                    for t in range(Nt):
                        if len(indexes_dict[field][t][h])!=el_ind_for_height_len[h]:
                            print('ERROR: height', h, 'mismatching lenind', len(indexes_dict[field][t][h]))
                        these_indices[:,t]=indexes_dict[field][t][h]
                    vp_ind[:] = these_indices

    if profile_type == "QVP":
        # add elevation attribute
        output.elevation = elevation
        output.excluded_azimuths = azimuth_exclude
        output.excluded_azimuth_number = len(azimuth_exclude)

    output.close()

    if verbose:
        print('New netcdf file created at {}'.format(output_file))
