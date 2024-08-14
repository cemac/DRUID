#! /usr/bin/env python3
from sys import argv
import pandas as pd

vp_script="/gws/smf/j04/ncas_radar/rrniii/BioDAR/VP_Extraction/vp_extraction.py"
#log output location
log_output="/gws/smf/j04/ncas_radar/rrniii/BioDAR/VP_Extraction/Output/"

data_input = "/gws/nopw/j04/ncas_radar_vol3/ukmo-nimrod/raw_h5_data/single-site/"
data_output = "/gws/smf/j07/ncas_radar/data/ukmo-nimrod/"

cvp_params_file = "./Parameters_CVP.xlsx"

# cvp_sites = {
             # 0: {'col_pos_name': "rothamsted",            'col_lat' : 51.806908, 'col_long' : -0.3609735, 'radar' : "chenies",       'col_radius' : [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]},
             # 1: {'col_pos_name': "rothamsted_fixed_15km", 'col_lat' : 51.78928,  'col_long' : -0.38672,   'radar' : "chenies",       'col_radius' : [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]},
             # 2: {'col_pos_name': "rothamsted_fixed_25km", 'col_lat' : 51.85451,  'col_long' : -0.28863,   'radar' : "chenies",       'col_radius' : [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]},
             # 3: {'col_pos_name': "rothamsted_fixed_35km", 'col_lat' : 51.92049,  'col_long' : -0.18684,   'radar' : "chenies",       'col_radius' : [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]},
# }
 #            4: {'col_pos_name': "silwood",               'col_lat' : 51.409459, 'col_long' : -0.6426696, 'radar' : "chenies",       'col_radius' : [2.5, 5, 7.5, 10]},
 #            5: {'col_pos_name': "silwood_fixed_15km",    'col_lat' : 51.55378,  'col_long' : -0.57675,   'radar' : "chenies",       'col_radius' : [2.5, 5, 7.5, 10]},
 #            6: {'col_pos_name': "silwood_fixed_25km",    'col_lat' : 51.46861,  'col_long' : -0.60742,   'radar' : "chenies",       'col_radius' : [2.5, 5, 7.5, 10]},
 #            7: {'col_pos_name': "silwood_fixed_35km",    'col_lat' : 51.37699,  'col_long' : -0.65068,   'radar' : "chenies",       'col_radius' : [2.5, 5, 7.5, 10]},
 #            8: {'col_pos_name': "wye",                   'col_lat' : 51.185132, 'col_long' : 0.9448683,  'radar' : "thurnham",      'col_radius' : [2.5, 5, 7.5, 10]},
 #            9: {'col_pos_name': "wye_fixed_15km",        'col_lat' : 51.22918,  'col_long' : 0.79166,    'radar' : "thurnham",      'col_radius' : [2.5, 5, 7.5, 10]},
 #           10: {'col_pos_name': "wye_fixed_25km",        'col_lat' : 51.18637,  'col_long' : 0.9314,     'radar' : "thurnham",      'col_radius' : [2.5, 5, 7.5, 10]},
 #           11: {'col_pos_name': "wye_fixed_35km",        'col_lat' : 51.14978,  'col_long' : 1.05019,    'radar' : "thurnham",      'col_radius' : [2.5, 5, 7.5, 10]},
 #           12: {'col_pos_name': "hereford",              'col_lat' : 52.124081, 'col_long' : -2.6382625, 'radar' : "clee-hill",     'col_radius' : [2.5, 5, 7.5, 10]},
 #           13: {'col_pos_name': "hereford_fixed_15km",   'col_lat' : 52.25514,  'col_long' : -2.6118,    'radar' : "clee-hill",     'col_radius' : [2.5, 5, 7.5, 10]},
 #           14: {'col_pos_name': "hereford_fixed_25km",   'col_lat' : 52.17268,  'col_long' : -2.62416,   'radar' : "clee-hill",     'col_radius' : [2.5, 5, 7.5, 10]},
 #           15: {'col_pos_name': "hereford_fixed_35km",   'col_lat' : 52.08205,  'col_long' : -2.64367,   'radar' : "clee-hill",     'col_radius' : [2.5, 5, 7.5, 10]},
 #           16: {'col_pos_name': "preston",               'col_lat' : 53.854688, 'col_long' : -2.7646454, 'radar' : "hameldon-hill", 'col_radius' : [2.5, 5, 7.5, 10]},
 #           17: {'col_pos_name': "preston_fixed_15km",    'col_lat' : 53.80162,  'col_long' : -2.51104,   'radar' : "hameldon-hill", 'col_radius' : [2.5, 5, 7.5, 10]},
 #           18: {'col_pos_name': "preston_fixed_25km",    'col_lat' : 53.82787,  'col_long' : -2.65404,   'radar' : "hameldon-hill", 'col_radius' : [2.5, 5, 7.5, 10]},
 #           19: {'col_pos_name': "preston_fixed_35km",    'col_lat' : 53.86671,  'col_long' : -2.79909,   'radar' : "hameldon-hill", 'col_radius' : [2.5, 5, 7.5, 10]},
 #           20: {'col_pos_name': "starcross",             'col_lat' : 50.629506, 'col_long' : -3.4548264, 'radar' : "cobbacombe",    'col_radius' : [2.5, 5, 7.5, 10]},
 #           21: {'col_pos_name': "starcross_fixed_15km",  'col_lat' : 50.82359,  'col_long' : -3.45459,   'radar' : "cobbacombe",    'col_radius' : [2.5, 5, 7.5, 10]},
 #           22: {'col_pos_name': "starcross_fixed_25km",  'col_lat' : 50.73279,  'col_long' : -3.45605,   'radar' : "cobbacombe",    'col_radius' : [2.5, 5, 7.5, 10]},
 #           23: {'col_pos_name': "starcross_fixed_55km",  'col_lat' : 50.63798,  'col_long' : -3.45654,   'radar' : "cobbacombe",    'col_radius' : [2.5, 5, 7.5, 10]},
 #       }
qvp_sites = {
             0: {'radar':'chenies',      'azimuth_exclude':None, 'count_threshold':1, 'elevation':[0.5, 1, 2, 3, 4]},
             1: {'radar':'clee-hill',    'azimuth_exclude':None, 'count_threshold':1, 'elevation':[0.5, 1, 2, 3, 4]},
             2: {'radar':'cobbacombe',   'azimuth_exclude':None, 'count_threshold':1, 'elevation':[0.5, 1, 2, 3, 4]},
             3: {'radar':'deamhill',     'azimuth_exclude':None, 'count_threshold':1, 'elevation':[0.5, 1, 2, 3, 4]},
             4: {'radar':'hameldon-hill','azimuth_exclude':None, 'count_threshold':1, 'elevation':[0.5, 1, 2, 3, 4]},
             5: {'radar':'high-moorsley','azimuth_exclude':None, 'count_threshold':1, 'elevation':[0.5, 1, 2, 3, 4]},
             6: {'radar':'thurnham',     'azimuth_exclude':None, 'count_threshold':1, 'elevation':[0.5, 1, 2, 3, 4]}
        }
# Azimuths to exclude should be put in in the format start1,end1;start2,end2;start3,end3....


def read_params(mode, path=None):
    if mode == "CVP":
        if path is None:
            path = cvp_params_file
        params = pd.read_excel(path, index_col=0)
        params['col_radius'] = params['col_radius'].apply(lambda x: x if isinstance(x, list) else [x])
        return params.to_dict(orient='index')
    elif mode == "QVP":
        return qvp_sites
    else:
        raise ValueError('VP mode not recognised')


def main():
    try:
        args = argv[1:]
    except IndexError:
        raise SystemExit(f"Usage: {argv[0]} <mode> [<param_file>]")

    mode = args[0]
    fp = args[1] if len(args) > 1 else None

    if mode == 'CVP':
        cvp_sites = read_params(mode, fp)
        print (len(cvp_sites))
    elif mode == 'QVP':
        print (len(qvp_sites))

if __name__ == '__main__':
    main()
