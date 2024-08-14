#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 
Column Vertical profile extraction of static indices for given radar site

This creates the static cvp_indexes for each cvp within the radar site and should be run as follows:
python vp_create_static_cvp_indices <input_dir> <output_dir> -p <params_file> -m -v -k <min_h> -j <max_h> -u <h_step>
where:
input_dir is where to find a radar file for this radar site, eg /gws/nopw/j04/ncas_radar_vol3/ukmo-nimrod/raw_h5_data/single-site/chenies/2020
output_dir is where to store the static_cvp_indices files eg /gws/smf/j07/ncas_radar/data/ukmo-nimrod/raw_cvp_data/chenies
params_file is an xlsx file containing the definition of the CVPs
-m=met_office
-v=verbose
min_h, max_h and h_step define the altitude levels for the CVP

This script was developed by CEMAC as part of the Drivers and
Repercussions of UK Insect Declines (DRUID) project (NE/V006916/1).

Authors:
 * Written by Julia crook, May 2024

:copyright: © 2024 University of Leeds.
:license: BSD3

"""

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
import vp_params

DEFAULT_MIN_H = 0
DEFAULT_MAX_H = 2000
DEFAULT_H_STEP = 200

def parse_args():
    formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=formatter)

    parser.add_argument('input_dir', type=str,
                        help="Path to input data")

    parser.add_argument('output_dir', type=str,
                        help="Path to output directory")

    parser.add_argument("-p", "--params_file",
                        dest="params_file", default=vp_params.cvp_params_file,
                        help=".xlsx file defining where CVPs are")

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help='''Print messages about the program execution
                        to the console (stdout)''')

    parser.add_argument("-m", "--met-office",
                        action="store_true",
                        help='''Indicate that the source of the data is
                        the Met office radars''')

    parser.add_argument("-k", "--column-minimum-altitude",
                            dest="min_h", default = DEFAULT_MIN_H,
                            help="Column lowest altitude in m")

    parser.add_argument("-j", "--column-maximum-altitude",
                            dest="max_h", default = DEFAULT_MAX_H,
                            help="Column highest altitude in m")

    parser.add_argument("-u", "--column-profile-resolution",
                            dest="h_step", default = DEFAULT_H_STEP,
                            help="Column resulting profile resolution in m")

    args = parser.parse_args()

    # Check if input directory exists
    if not os.path.exists(args.input_dir):
        err_msg = "Input dir {0} does not exist\n"
        err_msg = err_msg.format(args.input_dir)
        raise ValueError(err_msg)


    return args


def get_cvp_options(min_h,
                    max_h,
                    h_step,
                    verbose=False):
    '''
    Get options specific to CVP indices: min_h, max_h, h_step

    '''

    min_h = float(min_h)
    max_h = float(max_h)
    h_step = float(h_step)

    return (min_h, max_h, h_step)


def get_input_folder_glob_spec(input_dir,
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
        folder_glob_spec = '{}/*.nc'.format(folder_with_files)
    if verbose:
        print(("Input folder glob spec is {}").format(folder_glob_spec))

    return folder_glob_spec


def get_radar_filename(input_dir,
                       met_office,
                       verbose=False):
    '''
    Get file from which we can work out cvp static indexes.

    '''

    # Get file matching pattern
    folder_glob_spec = get_input_folder_glob_spec(input_dir,
                                                  met_office,
                                                  verbose)

    # select all files available in the directory
    file_list = glob.glob(folder_glob_spec)
    file_list.sort()

    match_file=file_list[0] # just choose the first file
        
    if verbose:
        print("match_file")
        print(match_file)

    return match_file

def main():

    #### Input management

    args = parse_args()
    
    # CVP specific options
    min_h, max_h, h_step = \
            get_cvp_options(args.min_h,
                            args.max_h,
                            args.h_step,
                            args.verbose)

    # radar file to read
    radar_fname = get_radar_filename(args.input_dir,
                                  args.met_office,
                                  args.verbose)
        
    # Output file path
    # Path to output directory
    output_dir = os.path.normpath(args.output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    

    #### Processing

    # calculate the heights in the cvp
    if args.verbose:
        print('min_h, max_h and h_step:', min_h, max_h, h_step)
    equidistant_alt = np.linspace((min_h + h_step/2), (max_h - h_step/2), num=int(max_h/h_step))
    equidistant_bound = np.linspace((min_h), (max_h), num=int(max_h/h_step)+1)
        
    # read the params file defining the cvps
    if args.verbose:
        print('reading',args.params_file)
    sites = vp_params.read_params('CVP', args.params_file)
    avg_range_delta=sites[0]['col_radius'][0] # not sure why this is read as an array
    
    # read a radar file
    (radar, unit_dict, long_names, short_names) = vp_io.read_file(0, '0000', radar_fname, [],met_office=args.met_office, verbose=False)
    # unfortunately the elevations are not always in the same order - they need to be for working out
    # static indices that are consistent across all times for this radar site cvp, so we need to sort them
    elevations = radar.elevation['data'].reshape((radar.nsweeps,int(radar.nrays/radar.nsweeps)))
    elevations, counts=np.unique(elevations, return_counts=True)
    ix=np.where(counts!=int(radar.nrays/radar.nsweeps))
    if len(ix[0])>0:
        print('WARNING: azimuths have different elevations!')
    # now sort them
    el_sort_ix=np.argsort(elevations)
    
    for s in range(len(sites)):
        this_site=sites[s]
        # Get CVP indexes for this CVP
        if args.verbose:
            print('getting static cvp indices for', this_site['col_pos_name'], this_site['col_lat'], this_site['col_long'], avg_range_delta)
        cvp_indexes = vp_functions.static_index_for_csv_file(
            radar,
            el_sort_ix,
            this_site['col_lat'],
            this_site['col_long'],
            avg_range_delta,
            equidistant_alt,
            equidistant_bound,
            verbose=args.verbose
        )
        cvp_name=this_site['col_pos_name']
        output_filename = '{}_{}km_{}_{}_{}_static_cvp_indices.nc'.format(cvp_name,
                                                                         avg_range_delta,
                                                                         min_h, max_h, h_step)
        vp_io.write_static_cvp_indices_netcdf(output_dir+'/'+output_filename, radar_fname, cvp_indexes, equidistant_alt)

    # end main()


if __name__ == "__main__":
    main()
