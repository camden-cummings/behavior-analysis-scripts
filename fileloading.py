#!/usr/bin/python

import os
import argparse
import numpy as np
import re
from Fish import Fish  # fish object
import matplotlib.pyplot as plt
from fileloading_labview import well_conversion
from fileloading_helpers import load_sections_file, load_highspeed_data, convert_to_polar, load_genotypes
from fileloading_python import load_python_data
from fileloading_labview import load_labview_data

# Input arguments
parser = argparse.ArgumentParser(description='loading for fish behavior files')
parser.add_argument('-python', action='store_true', dest='python', default=False)
parser.add_argument('-r', type=str, action="store", dest="rois_file")
parser.add_argument('-e', type=str, action="store", dest="events_file")
parser.add_argument('-c', type=str, action="store", dest="centroid_file")
parser.add_argument('-m', type=str, action="store", dest="movie_prefix", default="")
parser.add_argument('-g', type=str, action="store", dest="genotype_file")
parser.add_argument('-s', type=str, action="store", dest="sections_file", default="sectionsfile")

# labview
parser.add_argument('-t', type=str, action="store", dest="timestamp_file")
parser.add_argument('-d', type=str, action="store", dest="dpix_file")


parser.add_argument('-longmovie', type=str, action="store",
                    dest="long_movie_name", default="nomovie")
parser.add_argument('-outputmovies', action="store_true",
                    dest="output_movies", default=False)
# Tracked data from before code was updated to have output ROIs, irrelevant for new users, only compatible with 96-well plates
parser.add_argument('-oldtracking', action="store_true",
                    dest="old_tracking", default=False)
parser.add_argument('-xyhm', action="store_true", dest="xyhm", default=False)
parser.add_argument('-social', action="store_true",
                    dest="social", default=False)
parser.add_argument('-graphonly', action="store_true",
                    dest="graph_only", default=False)
# CURRENTLY NOT COMPATIBLE WITH STIMULI THAT NEED FILTERING
parser.add_argument('-graphmulti', type=str, action="store", dest="graph_multi")
parser.add_argument('-j', type=str, action="store",
                    dest="graph_parameters", default="PlotParameters")

parser.add_argument('-i', type=float, action="store",
                    dest="msec_per_frame", default=3.508772)

# not used in fileloading but used in processmotiondata
# The classic Schier (and now Prober) sleep plots are inactive min / hour (sleep) and active sec / hour (with sleep bouts not counted)
parser.add_argument('-a', type=str, action="store", dest="activity_times",
                    default="1/60,60/600,60/3600,1/3600")  # list of comparisons for activity data in seconds
# list of times bins for the bout data (ie, ave bout speed / minute) in seconds
parser.add_argument('-b', type=str, action="store",
                    dest="bout_bins", default="60,600,3600")
# List of thresholds for the distance and dpix calculations, which typically should not be touched. First is for distance, which is less robust than dpix, and after is the dpix value. This is the threshold value to be counted for a bout if greater than or equal to.
parser.add_argument('-v', type=str, action="store",
                    dest="threshold_values", default="0.5,3.0")
# List of high-speed thresholds for the distance and dpix calculations, which typically should not be touched. First is for distance, which is less robust than dpix, and after is the dpix value. This is the threshold value to be counted for a bout if greater than or equal to.
parser.add_argument('-w', type=str, action="store",
                    dest="hs_threshold_values", default="0.9,3.0")
# Same as above but for number of frames.
parser.add_argument('-f', type=str, action="store",
                    dest="threshold_frames", default="1,3")
# Same as above but for number of frames and high-speed data.
parser.add_argument('-x', type=str, action="store",
                    dest="hs_threshold_frames", default="2,3")
# thresholds for activity data, first distance and second dpix (differs from bout thresholds because we don't have frame considerations)
parser.add_argument('-y', type=str, action="store",
                    dest="activity_times_thresholds", default="1,10")
# ((boutrev > 4.0) and (0.3 < (boutspeed) < 1.3) and (boutdistance > 70)):
parser.add_argument('-z', type=str, action="store",
                    dest="seizure_filters", default="4.0,0.3,1.3,70")
# baseline level of light, used to determine what is a dark flash for filtering O-bend
parser.add_argument('-l', type=int, action="store",
                    dest="light_baseline", default=200)
# two measures that are intersected (both must be true), and responses with greater magnitude than both are considered true O-bends. The two measures are "responsetime" and "sumabsha" (sum of absolute value of heading angles)
parser.add_argument('-o', type=str, action="store", dest="obend_filter",
                    default="60,>:responsetime,10,>:responsesumabsheadingangle")
# two measures that are intersected (both must be true), and responses with greater magnitude than both are considered true C-bends. The two measures are "responsevelocity" and "responsecumdpix"
parser.add_argument('-p', type=str, action="store", dest="cbend_filter",
                    default="0.2,>:responsevelocity,1500,>:responsecumulativedpix")
parser.add_argument('-k', type=str, action="store",
                    dest="movie_filter", default="1,=:boutseizurecount")

args = parser.parse_args()
python = args.python
long_movie_name = args.long_movie_name
longmovie = False
if long_movie_name != "nomovie":
    longmovie = True
output_movies = args.output_movies

rois_file = args.rois_file
timestamp_file = args.timestamp_file

old_tracking = args.old_tracking
graph_only = args.graph_only
xyhm = args.xyhm
social = args.social
graph_parameters = args.graph_parameters
graph_multi = args.graph_multi
if not graph_only:
    events_file = args.events_file
    centroid_file = args.centroid_file
    dpix_file = args.dpix_file
    genotype_file = args.genotype_file
    sections_file = args.sections_file

movie_prefix = args.movie_prefix
msec_per_frame = args.msec_per_frame
# not used in fileloading but used in processmotiondata
activity_times = list(map(int, re.split(',|, |/|/', args.activity_times)))
# these bins and activity bins are not going to be less than a second (that doesn't work in code well), so it's fine to use int instead of float
bout_bins = list(map(int, args.bout_bins.split(',')))
threshold_values = list(map(float, args.threshold_values.split(',')))
hs_threshold_values = list(map(float, args.hs_threshold_values.split(',')))
threshold_frames = list(map(int, args.threshold_frames.split(',')))
hs_threshold_frames = list(map(int, args.hs_threshold_frames.split(',')))
activity_times_thresholds = list(
    map(float, args.activity_times_thresholds.split(',')))

seizure_filters = list(map(float, args.seizure_filters.split(',')))
light_baseline = args.light_baseline
obend_filter = list(map(str, args.obend_filter.split(',')))
cbend_filter = list(map(str, args.cbend_filter.split(',')))
movie_filter = list(map(str, args.movie_filter.split(',')))

# The Fish object carries all the data around for each fish, including their genotype and ID number
# Later (after processmotiondata.py) the data inside this Fish object is analyzed (bouts counted, binned, responses counted) and the AnalyzedFish object carries that data
# This analysis code only compares two groups: a control group and a test group
# Or it can analyze a single group, but no statistics will be done
def generate_fish_objects(dp_data_array, rho_array, theta_array, x_array, y_array, hs_dpix, hs_pos, rois_dict, num_of_wells):
    f = open(genotype_file, 'r')
    lines = f.readlines()
    f.close()
    genotype_list = {}
    for line in lines:
        # The file must use "controlgroup_geno: #,#,#" and "testgroup_geno: #,#,#,#" to emphasize to users that only two are allowed
        genogroup = line.split(':')[0].split('_')[0]
        if len(line.split(':')[0].split('_')[0]) > 1:
            realgenotype = line.split(':')[0].split('_')[1]
        else:
            if (line.split(':')[0].split('_')[0]) == "controlgroup":
                realgenotype = "control"
            elif (line.split(':')[0].split('_')[0]) == "testgroup":
                realgenotype = "test"
            else:
                print(
                    "Not using correct nomenclature of controlgroup_genotype: and testgroup_genotype:")
                realgenotype = line.split(':')[0].split('_')[0]
        fishidslist = line.split(':')[1].strip().split(',')
#        print(fishidslist)
        inputids = []
        for fish_id in fishidslist:
            print(fish_id)
            # Have to subtract 1 so that we can start from 0
            # Keeping it like this, but then adding 1 back to the ID that is saved
            inputids.append(int(fish_id) - 1)
        genotype_list[genogroup + "_" + realgenotype] = inputids
    fish_list = []
    for n in range(0, num_of_wells):
        split_hs_dpix = {}
        split_hs_pos_x = {}
        split_hs_pos_y = {}
        for d in hs_dpix.keys():
            if old_tracking:
                split_hs_dpix[d] = hs_dpix[d][:, well_conversion[n]]
                split_hs_pos_x[d] = hs_pos[d][:, 2 * well_conversion[n]]
                split_hs_pos_y[d] = hs_pos[d][:, 2 * well_conversion[n] + 1]
            else:
                split_hs_dpix[d] = hs_dpix[d][:, n]
                split_hs_pos_x[d] = hs_pos[d][:, 2 * n]
                split_hs_pos_y[d] = hs_pos[d][:, 2 * n + 1]
        for x in genotype_list.keys():
            if n in genotype_list[x]:
                # Adding 1 back onto the fish.idnumber, because all we use it for later is to connect to original input and we want it to match
                if social:
                    lower_bound = rois_dict[n + 1][0]
                    upper_bound = rois_dict[n + 1][2]
                    roimask = (x_array[:, n] < lower_bound) | (
                            x_array[:, n] > upper_bound)
                    x_array[:, n][roimask] = np.nan

                    if rois_dict[n + 1][0] < 400:
                        x_array[:, n] = x_array[:, n] - \
                                        np.nanmax(x_array[:, n])
                        x_array[:, n] = abs(-1 * x_array[:, n])
                    else:
                        x_array[:, n] = x_array[:, n] - \
                                        abs(np.nanmin(x_array[:, n]))

                newfish = Fish(n + 1, x.split('_')[0], x.split('_')[1], dp_data_array[:, n], rho_array[:, n],
                               theta_array[:, n], x_array[:, n], y_array[:, n], split_hs_dpix, split_hs_pos_x,
                               split_hs_pos_y)
                if xyhm:
                    nan_mask = np.logical_not(np.isnan(x_array[:, n]))
                    heatmap, xedges, yedges = np.histogram2d(
                        x_array[:, n][nan_mask], y_array[:, n][nan_mask])  #, bins=(35, 5))
                    fig, ax1 = plt.subplots()
                    ax1.imshow(
                        heatmap.T, interpolation='lanczos', cmap='jet')
                    plt.savefig("fish_" + str(n + 1) + "_XYheatmap.png",
                                transparent=True, format="png")
                    plt.close()
                if rois_file or longmovie:
                    if social:
                        rois_dict[n + 1][2] = rois_dict[n + 1][2] - \
                                              rois_dict[n + 1][0]
                        rois_dict[n + 1][0] = 0
                        newfish.add_rois(rois_dict[n + 1])
                    else:
                        newfish.add_rois(rois_dict[n + 1])
                fish_list.append(newfish)
    return fish_list

# Start here
def loading_procedures():
#    startdate = "6/3/2025"
#    start_time = datetime.datetime(2025, 6, 3, 18, 2, 5)
#    end_time = datetime.datetime(2025, 6, 5, 13, 36, 15)
#    events = load_sections_file(startdate, end_time, start_time, sections_file)

    if python:
        rois_dict, num_of_wells, cen_data_array, dp_data_array, tuple_timestamps = load_python_data(centroid_file, rois_file)
    else:
        rois_dict, num_of_wells, cen_data_array, dp_data_array, tuple_timestamps = load_labview_data(timestamp_file, rois_file, dpix_file, centroid_file, social, longmovie)

    print('num of wells', num_of_wells)

    dp_data_array.resize(dp_data_array.size //
                         num_of_wells, num_of_wells)

    cen_data_array = cen_data_array.reshape(
        cen_data_array.size // (num_of_wells * 2), (num_of_wells * 2))

    print("Done loading timestamp file")

    # just setting to zero to make it easier to ignore
    cen_data_array[cen_data_array == 65535] = 0
    # This is because they are all values of -16 and the initial type of the array is unsigned int, but it should be clear that it means the fish hasn't moved yet
    # converting them to zero for now so that it makes it easier to deal with the downstream max/min tests
    tuple_rho_theta = convert_to_polar(cen_data_array, num_of_wells, social)
    print("Done converting to polar coordinates")

    # This is when the run started, calculated from the timestamp file and it is useful for putting event data on correct day
    # Day "0" in the sections file corresponds to this start date from the very beginning of the timestamp file
    startdate = str(tuple_timestamps[0][0].month) + "/" + str(
        tuple_timestamps[0][0].day) + "/" + str(tuple_timestamps[0][0].year)
    end_time = tuple_timestamps[3]
    print('startdate', startdate, 'endDT', end_time)
    start_time = tuple_timestamps[4]
    print('startDT', start_time)

    events = load_sections_file(startdate, end_time, start_time, sections_file)
    hs_dpix, hs_pos = load_highspeed_data(startdate, events, msec_per_frame, movie_prefix, events_file, num_of_wells, python)

    #    global_tuple_events = load_event_data(startdate, endDT, startDT)
    print("Done loading events")

    fish_list = generate_fish_objects(dp_data_array, tuple_rho_theta[0], tuple_rho_theta[1], tuple_rho_theta[2],
                                      tuple_rho_theta[3], hs_dpix, hs_pos, rois_dict, num_of_wells)
    #    print(fish_list)
    #    print(global_tuple_events[2])
    return fish_list, tuple_timestamps[0], tuple_timestamps[1], tuple_timestamps[2], events, movie_prefix

def get_movie_prefix():
    return movie_prefix

def initialize_args():
    print("Initializing arguments")


