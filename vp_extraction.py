#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Vertical profile extraction

Extract QVP/CVP from radar scan data and store the profile in a
specified format.

This script was developed by CEMAC as part of the Drivers and
Repercussions of UK Insect Declines (DRUID) project (NE/V006916/1).

Authors:
 * Written by David Dufton, Nov. 2016
 * Adapted by M. Lukach, May 2018
 * Adapted by R. Neely, March 2022
 * Adapted by T.D. James <t.d.james1@leeds.ac.uk>, June 2022

:copyright: © 2022 University of Leeds.
:license: BSD3

"""
from __future__ import (absolute_import, division, print_function)
#from six.moves import (filter, input, map, range, zip)

import sys
import os
import re
import glob
import datetime as dt
import argparse

import numpy as np
from dateutil.parser import parse as dateparse

import vp_functions
import vp_io

DEFAULT_TSTART = '20180214T000000'
DEFAULT_TSTOP = '20180214T235959'
DEFAULT_FIELD_LIST = {
    'CVP': ['PhiDP', 'RhoHV', 'SQI', 'W', 'dBZ', 'ZDR', 'V'],
    'QVP': ['dBZ', 'ZDR', 'RhoHV', 'PhiDP', 'V', 'W', 'SQI']
}
DEFAULT_COUNT_THRESHOLD = 0
DEFAULT_AVG_RANGE_DELTA = 2.5
DEFAULT_COLUMN_LAT = 51.78928
DEFAULT_COLUMN_LON = -0.38672
DEFAULT_STATIC_POINT = "Rothamsted"
DEFAULT_MIN_H = 0
DEFAULT_MAX_H = 2000
DEFAULT_H_STEP = 200

def parse_args():
    formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter)

    parser.add_argument('profile_type', type=str,
                        choices=['QVP', 'CVP'], default="QVP",
                        help='''Profile type (QVP or CVP)''')

    parser.add_argument('input_dir', type=str,
                        help="Path to input data")

    parser.add_argument('output_dir', type=str,
                        help="Path to output directory")

    parser.add_argument("-d", "--debug",
                        action="store_true",
                        help=argparse.SUPPRESS)

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help='''Print messages about the program execution
                        to the console (stdout)''')

    parser.add_argument("-s", "--start-time",
                        dest="tstart", default=DEFAULT_TSTART,
                        help='''Start date and time (included) in format:
                        "YYYYMMDDThhmmss"''')

    parser.add_argument("-e", "--end-time",
                        dest="tstop", default=DEFAULT_TSTOP,
                        help='''End date and time (included) in format:
                        "YYYYMMDDThhmmss"''')

    parser.add_argument("-f", "--fields",
                        dest="field_list",
                        default=None,
                        help='''Space-separated list of variables to be
                        used in profile extraction''')

    parser.add_argument("-z", "--zoom-interval",
                        default=None,
                        help='''Time interval to zoom, specified as
                        comma-separated values:
                        \"YYYYMMDDThhmmss,YYYYMMDDThhmmss\"''')

    parser.add_argument("-m", "--met-office",
                        action="store_true",
                        help='''Indicate that the source of the data is
                        the Met office radars''')

    # Options for CVP extraction
    cvp_group = parser.add_argument_group('CVP', 'Options for CVP extraction')

    cvp_group.add_argument("-r", "--column-radius",
                           dest="avg_range_delta", default=DEFAULT_AVG_RANGE_DELTA,
                           help='''The radius of the column for CVP extraction (km)''')

    cvp_group.add_argument("-a", "--column-latitude",
                           dest="lat", default=DEFAULT_COLUMN_LAT,
                           help='''The latitude of the column for CVP extraction''')

    cvp_group.add_argument("-o", "--column-longitude",
                           dest="lon", default=DEFAULT_COLUMN_LON,
                           help='''The longitude of the column for CVP extraction''')

    cvp_group.add_argument("-p", "--column-position",
                           dest="static_point", type=str,
                           default=DEFAULT_STATIC_POINT,
                           help='''The name of the static position of the column
                           (station or lighting trap) for CVP extraction''')

    cvp_group.add_argument("-k", "--column-minimum-altitude",
                            dest="min_h", default = DEFAULT_MIN_H,
                            help="Column lowest altitude in m")

    cvp_group.add_argument("-j", "--column-maximum-altitude",
                            dest="max_h", default = DEFAULT_MAX_H,
                            help="Column highest altitude in m")

    cvp_group.add_argument("-u", "--column-profile-resolution",
                            dest="h_step", default = DEFAULT_H_STEP,
                            help="Column resulting profile resolution in m")

    cvp_group.add_argument("-i", "--store-column-indexes",
                           dest="store_indexes", action="store_true",
                           help="Store indexes of points contributing to CVP values")

    # Options for QVP extraction
    qvp_group = parser.add_argument_group('QVP', 'Options for QVP extraction')

    qvp_group.add_argument("-l", "--elevation", type=float,
                           help='''Elevation for QVP extraction''')

    qvp_group.add_argument("-c", "--count-threshold",
                           dest="count_threshold", default=DEFAULT_COUNT_THRESHOLD,
                           help='''Minimal number of points for the mean value
                           calculation at each range for QVP extraction''')

    qvp_group.add_argument("-b", "--azimuth-bounds-to-exclude",
                           dest="azimuth_bounds_to_exclude",
                           help='''Azimuths that should be excluded in the mean
                           value calculation for QVP extraction, specified as:
                           \"start1,end1;start2,end2\"''')

    args = parser.parse_args()

    # Check if input directory exists
    if not os.path.exists(args.input_dir):
        err_msg = "Input dir {0} does not exist\n"
        err_msg = err_msg.format(args.input_dir)
        raise ValueError(err_msg)

    if args.profile_type == "QVP" and args.elevation == None:
        raise argparse.ArgumentError(None,
                                     "argument --elevation is required for QVP extraction")

    return args


def get_field_list(field_list,
                   profile_type,
                   verbose=False):
    '''
    Get list of fields to include in profile.

    '''

    if field_list is None:
        field_list = DEFAULT_FIELD_LIST[profile_type]

    if not type(field_list) is list:
        field_list = field_list.split(',')

    return field_list


def get_event_date(tstart,
                   tstop,
                   verbose=False):
    '''
    Get event date string based on the specified start time.

    Also returns start and end times as datetime objects.
    '''

    start_datetime = dateparse(tstart)
    stop_datetime = dateparse(tstop)

    if verbose:
        print("start datetime is ", start_datetime)
        print("stop datetime is ", stop_datetime)

    # Define event date
    date_fmt = '%Y%m%d'
    event_date = start_datetime.strftime(date_fmt)

    # Check if start and end datetimes belong to the same day
    start_date = dt.datetime.strptime(event_date, date_fmt)
    stop_date = dt.datetime.strptime(stop_datetime.strftime(date_fmt),
                                     date_fmt)

    if not stop_date - start_date < dt.timedelta(1):
        # TODO define and handle a list of event dates
        if verbose:
            print("Start and end dates span more than one day")

    return (start_datetime, stop_datetime, event_date)


def get_zoom_interval(zoom_interval,
                      start_datetime,
                      stop_datetime,
                      verbose=False):
    '''
    Parse zoom interval and return start and end datetimes for the
    interval.

    '''
    if type(zoom_interval) is str:
        zoom_interval = tuple(zoom_interval.split(","))

    if zoom_interval and len(zoom_interval) == 2:
        zoom_start = dateparse(zoom_interval[0])
        zoom_end = dateparse(zoom_interval[1])
    else:
        if verbose:
            print('''Zoom interval should be specified as a list of two dates.
            Using start and end dates instead.''')
        zoom_start = start_datetime
        zoom_end = stop_datetime

    return zoom_start, zoom_end


def get_qvp_options(elevation,
                    count_threshold,
                    azimuth_bounds_to_exclude,
                    verbose=False):
    '''
    Get options specific to QVP extraction: elevation, count
    threshold and azimuth bounds to exclude.

    '''
    # elevation
    elevation = float(elevation)
    if verbose:
        print("elevation is ", elevation)

    # minimal number of valuable pixels in each range
    count_threshold = int(count_threshold)
    if verbose:
        print("count threshold is ", count_threshold)

    # exclude parts of domain that are partially or completely blocked
    azimuth_exclude = []
    if azimuth_bounds_to_exclude != None:
        azimuth_bounds = [x.split(',') for x in azimuth_bounds_to_exclude.split(';')]
        azimuth_chunks = []
        for chunk in azimuth_bounds:
            if len(chunk) != 2:
                if verbose:
                    print ('''Azimuth exclusion configuration incorrect.
                    Should be input as "start1,end1;start2,end2;..."\n
                    Skipping azimuth exclusion chunk: %s''' % chunk)
            else:
                index = np.sort(np.array(chunk, dtype=int))
                if index[0] == index[1]:
                    azimuth_chunks.append(np.array([index[0]]))
                else:
                    azimuth_chunks.append(np.arange(index[0], index[1]))
        if len(azimuth_chunks):
            azimuth_exclude = np.concatenate(azimuth_chunks)
    if verbose:
        print("azimuth values to exclude are ", azimuth_exclude)

    return (elevation, count_threshold, azimuth_exclude)


def get_cvp_options(avg_range_delta,
                    lat,
                    lon,
                    static_point,
                    min_h,
                    max_h,
                    h_step,
                    store_indexes,
                    verbose=False):
    '''
    Get options specific to CVP extraction: average range delta,
    latitude, longitude and static point label.

    '''

    avg_range_delta = float(avg_range_delta)
    lat = float(lat)
    lon = float(lon)
    static_point = static_point
    min_h = float(min_h)
    max_h = float(max_h)
    h_step = float(h_step)
    store_indexes = store_indexes

    return (avg_range_delta, lat, lon, static_point, min_h, max_h, h_step, store_indexes)


def get_input_folder_glob_spec(input_dir,
                               profile_type,
                               elevation,
                               event_date,
                               met_office,
                               verbose=False):
    '''
    Get pattern for matching files in input folder.

    The matching pattern depends on profile type and source of radar scan data.
    '''

    # Path to the directory with the radar scans
    folder_with_files = os.path.normpath(input_dir)
    if not os.path.isdir(folder_with_files):
        folder_with_files = os.path.dirname(folder_with_files)
    if verbose:
        print("Path to input directory is", folder_with_files)

    # Handle different options for input directory structure
    if met_office:
        folder_glob_spec = '{}/*.h5'.format(folder_with_files)
    else:
        if profile_type == 'QVP' and elevation == 90.0:
            folder_glob_spec = '{}/{}/ver/*.nc'.format(folder_with_files, event_date)
        else:
            folder_glob_spec = '{}/{}/*.nc'.format(folder_with_files, event_date)
    if verbose:
        print(("Input folder glob spec is {}").format(folder_glob_spec))

    return folder_glob_spec


def get_file_list(input_dir,
                  profile_type,
                  elevation,
                  event_date,
                  start_datetime,
                  stop_datetime,
                  met_office,
                  verbose=False):
    '''
    Get list of files for profile extraction.

    '''

    # Get file matching pattern
    folder_glob_spec = get_input_folder_glob_spec(input_dir,
                                                  profile_type,
                                                  elevation,
                                                  event_date,
                                                  met_office,
                                                  verbose)

    # select all files available in the directory
    file_list = glob.glob(folder_glob_spec)
    file_list.sort()

    if (met_office):
        # match files in specified date range
        match_file_list = []
        for f in file_list:
            file = os.path.basename(f)
            regex = r'(\d{8})(_)|(\d{6})(\d{6})|(\d{8})(-*)(\d{6})'
            timestamp_search = re.search(regex, file, re.IGNORECASE)
            if timestamp_search:
                # TODO check selection of date string
                if timestamp_search.group(1) != None:
                    timestr = timestamp_search.group(1)
                elif timestamp_search.group(3) != None:
                    timestr = timestamp_search.group(3)
                elif timestamp_search.group(5) != None:
                    timestr = timestamp_search.group(5)
                timestamp = dateparse(timestr)
                if timestamp >= start_datetime and timestamp <= stop_datetime:
                    match_file_list.append(f)
        file_list = match_file_list

    if verbose:
        print("file_list")
        print(file_list)
        if file_list:
            print('First file is:', file_list[0])
            print('Last file is:', file_list[-1])
            print("There are ", len(file_list), " files")

    return file_list


def get_output_filepath(output_dir,
                        output_filename,
                        profile_type,
                        event_date,
                        verbose=False):
    '''
    Get output file path.

    '''

    # Path to output directory
    output_dir = os.path.normpath(output_dir)

    if profile_type == 'QVP':
        output_dir = '{}/{}_QVP'.format(output_dir, event_date)

    if verbose:
        print("Path to output directory is", output_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # path and filename for the output file
    output_file = '{}/{}'.format(output_dir, output_filename)

    if verbose:
        print("Output will be placed in {}".format(output_file))

    return output_file


def main():

    #### Input management

    args = parse_args()

    debug = args.debug
    verbose = args.verbose

    # True if the data are from the met office radar
    met_office = args.met_office

    # Profile type (QVP or CVP)
    profile_type = args.profile_type

    # List of variables in original scan that will be used in profile
    # extraction
    field_list = get_field_list(args.field_list,
                                profile_type,
                                verbose)

    # Event date(s)
    start_datetime, stop_datetime, event_date = get_event_date(args.tstart,
                                                               args.tstop,
                                                               verbose)

    # Zoom interval
    zoom_start, zoom_end = get_zoom_interval(args.zoom_interval,
                                             start_datetime,
                                             stop_datetime,
                                             verbose)

    if profile_type == 'QVP':
        # QVP specific options
        elevation, count_threshold, azimuth_exclude = \
            get_qvp_options(args.elevation,
                            args.count_threshold,
                            args.azimuth_bounds_to_exclude,
                            verbose)

        # Filename for the output QVP file
        # e.g. 20170517_QVP_20.0deg.nc
        output_filename = '{}_QVP_{:.1f}deg.nc'.format(event_date, elevation)

    elif profile_type == 'CVP':
        # CVP specific options
        avg_range_delta, lat, lon, static_point, min_h, max_h, h_step, store_indexes = \
            get_cvp_options(args.avg_range_delta,
                            args.lat,
                            args.lon,
                            args.static_point,
                            args.min_h,
                            args.max_h,
                            args.h_step,
                            args.store_indexes,
                            verbose)

        # Filename for the output CVP file
        # e.g. Rothamsted_10km_20170517.nc
        output_filename = '{}_{}km_{}.nc'.format(static_point,
                                                 avg_range_delta,
                                                 event_date)

        # Dummy values for elevation and azimuth_exclude
        elevation = 0
        azimuth_exclude = []

        # For static CVP, set alt to zero
        alt = 0

    # File list
    file_list = get_file_list(args.input_dir,
                              profile_type,
                              elevation,
                              event_date,
                              start_datetime,
                              stop_datetime,
                              met_office,
                              verbose)

    # Output file path
    output_filepath = get_output_filepath(args.output_dir,
                                          output_filename,
                                          profile_type,
                                          event_date,
                                          verbose)

    #### Processing

    if profile_type == "QVP":
        # Extract QVP
        data_list = vp_functions.time_height(
            file_list,
            field_list,
            elevation=elevation,
            count_threshold=count_threshold,
            azimuth_exclude=azimuth_exclude,
            met_office=met_office,
            verbose=verbose,
            vp_mode='qvp'
        )
    elif profile_type == "CVP":

        # Get CVP indexes
        static_cvp_indices_filename = '{}_{}km_{}_{}_{}_static_cvp_indices.nc'.format(static_point,
                                                                         avg_range_delta,
                                                                         min_h, max_h, h_step)
        # if output dir has a year at the end then remove that to get to directory where cvp
        # static indices are stored, otherwise use the same directory
        wsplit=args.output_dir.split('/')
        if wsplit[-2]==start_datetime.strftime('%Y'):
            # directory ends with year so remove it
            s='/'
            static_dir=s.join(wsplit[:-2])+'/'
        else:
            static_dir=args.output_dir
        static_cvp_indx_pathname=static_dir+'/'+static_cvp_indices_filename
        static_cvp=vp_io.read_static_cvp_indices(static_cvp_indx_pathname)

        # Extract CVP
        data_list = vp_functions.time_height(
            file_list,
            field_list,
            cvp_indexes=static_cvp.get_indxs(),
            equidistant_alt=static_cvp.heights,
            met_office=met_office,
            store_indexes = store_indexes,
            verbose=verbose,
            vp_mode='cvp_static'
        )
    else:
        raise ValueError("Unrecognised profile type")

    #### Output management

    # Output to NetCDF
    vp_io.output_netcdf(data_list,
                            output_filepath,
                            profile_type,
                            field_list,
                            elevation,
                            azimuth_exclude,
                            len(file_list),
                            verbose=verbose)

    # end main()


if __name__ == "__main__":
    main()
