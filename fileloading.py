#!/usr/bin/python

import os
import glob
import argparse
import numpy as np
import re
import datetime
import os.path
from os import path
from Fish import Fish  # fish object
from EventSection import EventSection  # event section object
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import math
from fileloading_labview import well_conversion


# Input arguments
parser = argparse.ArgumentParser(description='loading for fish behavior files')
parser.add_argument('-longmovie', type=str, action="store",
                    dest="longmoviename", default="nomovie")
parser.add_argument('-outputmovies', action="store_true",
                    dest="outputmovies", default=False)
parser.add_argument('-r', type=str, action="store", dest="roisfile")
# Tracked data from before code was updated to have output ROIs, irrelevant for new users, only compatible with 96-well plates
parser.add_argument('-oldtracking', action="store_true",
                    dest="oldtracking", default=False)
parser.add_argument('-xyhm', action="store_true", dest="xyhm", default=False)
parser.add_argument('-social', action="store_true",
                    dest="social", default=False)
parser.add_argument('-graphonly', action="store_true",
                    dest="graphonly", default=False)
# CURRENTLY NOT COMPATIBLE WITH STIMULI THAT NEED FILTERING
parser.add_argument('-graphmulti', type=str, action="store", dest="graphmulti")
parser.add_argument('-j', type=str, action="store",
                    dest="graphparameters", default="PlotParameters")

parser.add_argument('-e', type=str, action="store", dest="eventsfile")
parser.add_argument('-c', type=str, action="store", dest="centroidfile")
parser.add_argument('-d', type=str, action="store", dest="dpixfile")
parser.add_argument('-m', type=str, action="store",
                    dest="movieprefix", default="")
parser.add_argument('-g', type=str, action="store", dest="genotype_file")
parser.add_argument('-s', type=str, action="store",
                    dest="sectionsfile", default="sectionsfile")
parser.add_argument('-n', type=int, action="store",
                    dest="num_of_wells", default=96)
parser.add_argument('-i', type=float, action="store",
                    dest="msecperframe", default=3.508772)

# not used in fileloading but used in processmotiondata
# The classic Schier (and now Prober) sleep plots are inactive min / hour (sleep) and active sec / hour (with sleep bouts not counted)
parser.add_argument('-a', type=str, action="store", dest="activitytimes",
                    default="1/60,60/600,60/3600,1/3600")  # list of comparisons for activity data in seconds
# list of times bins for the bout data (ie, ave bout speed / minute) in seconds
parser.add_argument('-b', type=str, action="store",
                    dest="boutbins", default="60,600,3600")
# List of thresholds for the distance and dpix calculations, which typically should not be touched. First is for distance, which is less robust than dpix, and after is the dpix value. This is the threshold value to be counted for a bout if greater than or equal to.
parser.add_argument('-v', type=str, action="store",
                    dest="thresholdvalues", default="0.5,3.0")
# List of high-speed thresholds for the distance and dpix calculations, which typically should not be touched. First is for distance, which is less robust than dpix, and after is the dpix value. This is the threshold value to be counted for a bout if greater than or equal to.
parser.add_argument('-w', type=str, action="store",
                    dest="hsthresholdvalues", default="0.9,3.0")
# Same as above but for number of frames.
parser.add_argument('-f', type=str, action="store",
                    dest="thresholdframes", default="1,3")
# Same as above but for number of frames and high-speed data.
parser.add_argument('-x', type=str, action="store",
                    dest="hsthresholdframes", default="2,3")
# thresholds for activity data, first distance and second dpix (differs from bout thresholds because we don't have frame considerations)
parser.add_argument('-y', type=str, action="store",
                    dest="activitytimesthresholds", default="1,10")
# ((boutrev > 4.0) and (0.3 < (boutspeed) < 1.3) and (boutdistance > 70)):
parser.add_argument('-z', type=str, action="store",
                    dest="seizurefilters", default="4.0,0.3,1.3,70")
# baseline level of light, used to determine what is a dark flash for filtering O-bend
parser.add_argument('-l', type=int, action="store",
                    dest="lightbaseline", default=200)
# two measures that are intersected (both must be true), and responses with greater magnitude than both are considered true O-bends. The two measures are "responsetime" and "sumabsha" (sum of absolute value of heading angles)
parser.add_argument('-o', type=str, action="store", dest="obendfilter",
                    default="60,>:responsetime,10,>:responsesumabsheadingangle")
# two measures that are intersected (both must be true), and responses with greater magnitude than both are considered true C-bends. The two measures are "responsevelocity" and "responsecumdpix"
parser.add_argument('-p', type=str, action="store", dest="cbendfilter",
                    default="0.2,>:responsevelocity,1500,>:responsecumulativedpix")
parser.add_argument('-k', type=str, action="store",
                    dest="moviefilter", default="1,=:boutseizurecount")


args = parser.parse_args()
longmoviename = args.longmoviename
longmovie = False
if longmoviename != "nomovie":
    longmovie = True
outputmovies = args.outputmovies
roisfile = args.roisfile
oldtracking = args.oldtracking
graphonly = args.graphonly
xyhm = args.xyhm
social = args.social
graphparameters = args.graphparameters
graphmulti = args.graphmulti
if not graphonly:
    eventsfile = args.eventsfile
    centroidfile = args.centroidfile
    dpixfile = args.dpixfile
    movieprefix = args.movieprefix
    genotype_file = args.genotype_file
    sectionsfile = args.sectionsfile

num_of_wells = args.num_of_wells
msecperframe = args.msecperframe
# not used in fileloading but used in processmotiondata
activitytimes = list(map(int, re.split(',|, |/|/', args.activitytimes)))
# these bins and activity bins are not going to be less than a second (that doesn't work in code well), so it's fine to use int instead of float
boutbins = list(map(int, args.boutbins.split(',')))
thresholdvalues = list(map(float, args.thresholdvalues.split(',')))
hsthresholdvalues = list(map(float, args.hsthresholdvalues.split(',')))
thresholdframes = list(map(int, args.thresholdframes.split(',')))
hsthresholdframes = list(map(int, args.hsthresholdframes.split(',')))
activitytimesthresholds = list(
    map(float, args.activitytimesthresholds.split(',')))

seizurefilters = list(map(float, args.seizurefilters.split(',')))
lightbaseline = args.lightbaseline
obendfilter = list(map(str, args.obendfilter.split(',')))
cbendfilter = list(map(str, args.cbendfilter.split(',')))
moviefilter = list(map(str, args.moviefilter.split(',')))

# Convert cartesian coordinates to polar
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return rho, theta

# Loading of the centroid position data for the high-speed movies
def load_movie_motion_pos(hs_pos, thistime, fn):
    # At this point I actually don't care if the data is centered
    # All I want it for is distances
    moviename = fn[:-12] + ".avi.centroid2"
    cenfile = open(moviename, 'r')
    hscen_data_list = []
    lines = cenfile.readlines()
    for line in lines:
        if "\n" != line:
           hscen_data_list.append(int(line))
    hscen_data_array = np.array(hscen_data_list)
    hscen_data_array.resize(
        hscen_data_array.size // (num_of_wells*2), (num_of_wells*2))

    hs_pos[thistime] = hscen_data_array

# Loading of the delta pixel motion data for the high-speed movies
def load_movie_motion(hs_dpix, thistime, moviename):
    if not path.exists(moviename):
        print("Not all high-speed movies are tracked. Please check. Exiting.")
        exit()
    firstdpix = open(moviename, 'r')
    hsdp_data_list = []
    lines = firstdpix.readlines()
    for line in lines:
        hsdp_data_list.append(int(line))

    hsdp_data_array = np.array(hsdp_data_list)

    hsdp_data_array = hsdp_data_array.reshape(
        hsdp_data_array.size // num_of_wells, num_of_wells)
    hs_dpix[thistime] = hsdp_data_array

# Loading in the high-speed movies and the sections we will analyze (sectionsfile)
def load_sections_file(startdate, end_time, start_time):
    # Load in the original events file
    # Tab separation vs spacing is KEY in the events file (LabView will fail if not right, and so will this code)
    events = []
    # hs_dpix is a dictionary of numpy array for each high-speed movie event
    s = open(sectionsfile, 'r')
    slines = s.readlines()
    for sline in slines:
        # The EventSection object contains all of the high-speed events that are linked to each other
        # For example, the prepulse test has both strong tap, weak tap, and prepulse tap events in that section
        # Sections are designated by the user in the sections file and have a specific naming convention
        try:
            eventsection = EventSection(sline.strip(), startdate)
        except:
            print("The format of the input of your sections files is wrong. Every line should be as follows: habituation_daytaps=1_15:25:30-1_16:59:00")
            print("Exiting")
            exit()
        # Makes sure that the last section is not after the true end of the data
        eventsection.check_endtime(end_time)
        # Makes sure that the last section is not after the true end of the data
        eventsection.check_starttime(start_time)
        events.append(eventsection)

    return events

def load_highspeed_data(startdate, events):
    f = open(eventsfile, 'r')
    lines = f.readlines()

    hscounter = 0
    # 1:06:24\tPM\t0\t2\n'
    # Getting list of all the ID numbers associated with just the high-speed short movies
    highspeedlist = []
    for file in glob.glob(movieprefix + '*motion2'):
        highspeedlist.append(file)
    
    highspeedlist.sort(key=lambda name: int(name.split('/')[-1].split('-')[1].split('.')[0]))

    lastAMorPM0 = None
    lastAMorPMcounter0 = 0
    hs_dpix = {}
    hs_pos = {}

    for line in lines:
        TAPES = line.split('	')
        print(TAPES)
        # Need to be splitting on tab only for the new format of the input file with the spaces in last column
        # Converting date to the format that we were using for the
        # 6/18/2016_12:59:34_PM
        dateplustime = startdate + "_" + \
                       TAPES[0][0:len(TAPES[0])] + "_" + TAPES[1]
        thistime = datetime.datetime.strptime(
            dateplustime, '%m/%d/%Y_%I:%M:%S_%p')
        thisAMorPM0 = TAPES[1]
        if lastAMorPM0 is not None:
            # If we transition either AM->PM or PM->AM
            # The only time we will ever have issues with this is if the events file doesn't correspond well with beginning timestamp input
            # Or there could be issues with start after midnight
            # These situations should not happen
            # Basically things become a mess if the program does not read the events file in as it should have and it skips starting at beginning
            if thisAMorPM0 != lastAMorPM0:
                if thisAMorPM0 == "AM":
                    lastAMorPMcounter0 = lastAMorPMcounter0 + 1
        lastAMorPM0 = thisAMorPM0
        thistime = thistime + datetime.timedelta(days=lastAMorPMcounter0)
#        print(thistime, events)
        for eventsection in events:
#            print(thistime, eventsection.starttime, eventsection.endtime)

            if eventsection.starttime <= thistime <= eventsection.endtime:
                eventsection.add_event(TAPES[2], TAPES[3].strip(';'), thistime, msecperframe)
#                print("!")

        if int(TAPES[2]) != 0 and hscounter < len(highspeedlist):
            # Load in only the high-speed movie data marked '1' and with the high-speed movie prefix
            # This would need to be modified if a labview code other than '1' was used to make high-speed movies that should analyzed in this way
            # Would also need to modify code in processmotiondata.py, specifically the CalculateEventProperties function relies on this "1" to trigger analysis
            if int(TAPES[2]) == 1:
#                print("-----", highspeedlist[hscounter])
                load_movie_motion(hs_dpix, thistime, highspeedlist[hscounter])
                load_movie_motion_pos(hs_pos, thistime, highspeedlist[hscounter])
                hscounter = hscounter + 1

    if hscounter != len(highspeedlist):
        print(
            "ERROR, the number of high-speed movies does not correspond to the number expected based on the event file")
        print(
            "Please make sure that the event file is accurate and run was fully completed. Modify events file to reflect reality of data output if needed")
        print("Exiting")

    return hs_dpix, hs_pos

# Find max and min value for each fish in order to identify well edges
# This is actually not really the well edges, but the edges of the fish's radius of movement
# I prefer this approach to reading in the originally designated ROI, just in case ROI isn't accurate (ie, includes extra plastic of well edge)
def max_min(cen_data_array):
    maxxys = []
    minxys = []
    for n in range(0, num_of_wells*2, 2):
        maxtest = np.amax(cen_data_array[:, n])
        mintest = np.amin(cen_data_array[:, n])
        # This is the check for if nothing ever moves (the well is empty)
        if maxtest == mintest and maxtest == 0:
            maxxys.append(0)
            maxxys.append(0)
            minxys.append(0)
            minxys.append(0)
        else:
            maxrealx = maxtest
            minrealx = np.amin(
                cen_data_array[:, n][np.nonzero(cen_data_array[:, n])])
            maxrealy = np.amax(cen_data_array[:, n+1])
            minrealy = np.amin(
                cen_data_array[:, n+1][np.nonzero(cen_data_array[:, n+1])])
            maxxys.append(maxrealx)
            maxxys.append(maxrealy)
            minxys.append(minrealx)
            minxys.append(minrealy)
    maxxysnp = np.array(maxxys)
    minxysnp = np.array(minxys)
    return maxxysnp, minxysnp

# Polar coordinates are essential for easy calculation of well edge/center preferences
def convert_to_polar(cen_data_array):
    (maxxysnp, minxysnp) = max_min(cen_data_array)
    midcoords = (maxxysnp + minxysnp) / 2
    midcoords = midcoords.astype(np.int16)
    cen_data_array = cen_data_array.astype(np.int16)
    # just setting to very low value to make it easier to skip later
    cen_data_array[cen_data_array == 0] = -10000
    # subtract middle coordinate to get everything centered about 0
    zerodcoords = np.zeros(np.shape(cen_data_array))
    for i in range(0, num_of_wells*2):
        zerodcoords[:, i] = cen_data_array[:, i] - midcoords[i]
    zerodcoords[zerodcoords < -5000] = 0
    # zerodcoords currently contains negative numbers, which I think mean that the fish hasn't moved yet
    thetadata = np.zeros((len(cen_data_array), num_of_wells))
    rhodata = np.zeros((len(cen_data_array), num_of_wells))
    xzerod = np.zeros((len(cen_data_array), num_of_wells))
    yzerod = np.zeros((len(cen_data_array), num_of_wells))
    for i in range(0, num_of_wells):
        (rhodata[:, i], thetadata[:, i]) = cart2pol(
            zerodcoords[:, 2*i], zerodcoords[:, 2*i+1])
        xzerod[:, i] = zerodcoords[:, 2*i]
        yzerod[:, i] = zerodcoords[:, 2*i+1]
        if social:
            # RETURNING X/Y NOT ZEROED!!
            xzerod[:, i] = cen_data_array[:, 2*i]
            yzerod[:, i] = cen_data_array[:, 2*i+1]
    return rhodata, thetadata, xzerod, yzerod

# The Fish object carries all the data around for each fish, including their genotype and ID number
# Later (after processmotiondata.py) the data inside this Fish object is analyzed (bouts counted, binned, responses counted) and the AnalyzedFish object carries that data
# This analysis code only compares two groups: a control group and a test group
# Or it can analyze a single group, but no statistics will be done
def generate_fish_objects(dp_data_array, rho_array, theta_array, x_array, y_array, hs_dpix, hs_pos, rois_dict):
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
        print(fishidslist)
        inputids = []
        for fish_id in fishidslist:
            print(fish_id)
            # Have to subtract 1 so that we can start from 0
            # Keeping it like this, but then adding 1 back to the ID that is saved
            inputids.append(int(fish_id)-1)
        genotype_list[genogroup + "_" + realgenotype] = inputids
    fish_list = []
    for n in range(0, num_of_wells):
        split_hs_dpix = {}
        split_hs_pos_x = {}
        split_hs_pos_y = {}
        for d in hs_dpix.keys():
            if oldtracking:
                split_hs_dpix[d] = hs_dpix[d][:, well_conversion[n]]
                split_hs_pos_x[d] = hs_pos[d][:, 2*well_conversion[n]]
                split_hs_pos_y[d] = hs_pos[d][:, 2*well_conversion[n]+1]
            else:
                split_hs_dpix[d] = hs_dpix[d][:, n]
                split_hs_pos_x[d] = hs_pos[d][:, 2*n]
                split_hs_pos_y[d] = hs_pos[d][:, 2*n+1]
        for x in genotype_list.keys():
            if n in genotype_list[x]:
                # Adding 1 back onto the fish.idnumber, because all we use it for later is to connect to original input and we want it to match
                if social:
                    lower_bound = rois_dict[n+1][0]
                    upper_bound = rois_dict[n+1][2]
                    roimask = (x_array[:, n] < lower_bound) | (
                        x_array[:, n] > upper_bound)
                    x_array[:, n][roimask] = np.nan

                    if rois_dict[n+1][0] < 400:
                        x_array[:, n] = x_array[:, n] - \
                            np.nanmax(x_array[:, n])
                        x_array[:, n] = abs(-1 * x_array[:, n])
                    else:
                        x_array[:, n] = x_array[:, n] - \
                            abs(np.nanmin(x_array[:, n]))

                newfish = Fish(n + 1, x.split('_')[0], x.split('_')[1], dp_data_array[:, n], rho_array[:, n],
                               theta_array[:, n], x_array[:, n], y_array[:, n], split_hs_dpix, split_hs_pos_x, split_hs_pos_y)
                if xyhm:
                    nan_mask = np.logical_not(np.isnan(x_array[:, n]))
                    heatmap, xedges, yedges = np.histogram2d(
                        x_array[:, n][nan_mask], y_array[:, n][nan_mask])#, bins=(35, 5))
                    fig, ax1 = plt.subplots()
                    ax1.imshow(
                        heatmap.T, interpolation='lanczos', cmap='jet')
                    plt.savefig("fish_" + str(n+1) + "_XYheatmap.png",
                                transparent=True, format="png")
                    plt.close()
                if roisfile or longmovie:
                    if social:
                        rois_dict[n+1][2] = rois_dict[n+1][2] - \
                            rois_dict[n+1][0]
                        rois_dict[n+1][0] = 0
                        newfish.add_rois(rois_dict[n+1])
                    else:
                        newfish.add_rois(rois_dict[n+1])
                fish_list.append(newfish)
    return fish_list

def get_dists_and_posns(no_dups):
    posns_by_frame = []
    dists_by_frame = []

    indic = [i for i in range(0, int(len(no_dups)), num_of_wells)]

    for count in range(len(indic) - 2):
        frame_arr = no_dups.iloc[indic[count]:indic[count + 1]]
        next_frame_arr = no_dups.iloc[indic[count + 1]:indic[count + 2]]

        xs = list(frame_arr['pos_x'])
        ys = list(frame_arr['pos_y'])
        next_xs = list(next_frame_arr['pos_x'])
        next_ys = list(next_frame_arr['pos_y'])

        frame_xys = np.array(list(zip(xs, ys))).flatten()

        dists = []
        for i in range(len(xs)):
            dists.append(math.dist((xs[i], ys[i]), (next_xs[i], next_ys[i])))

        posns_by_frame.append(frame_xys)
        dists_by_frame.append(dists)

    return np.array(posns_by_frame), np.array(dists_by_frame)


def get_timestamps_from_csv(no_dups):
    frame_nums = list(no_dups['frame'].dropna())
    prev_time = "-1_-1_-1"  # times[0]
    #    print(prev_time)
    prev_hour = 0
    year = 2025
    month = 6
    day_count = 3

    count = 0
    frame_start = int(frame_nums[0])
    final_frame = int(frame_nums[-1])
    print("first frame", frame_start)
    print("final frame", final_frame)
    print("len frames", len(frame_nums))
    hour, minute, second = prev_time.split("_")
    overalltime = f"{year}_{month}_{day_count}_{hour}_{minute}_{second}"

    prev_second = 0
    timestamp_data_dict = {}
    timestamp_data_array = []
    ms_timestamp_data_array = []

    dropped_seconds = []

    for i in range(frame_start, final_frame, 1000):
        curr_first_frame = i
        if i + 999 < final_frame:
            curr_final_frame = i + 999
        else:
            curr_final_frame = final_frame

        start_ind = no_dups[no_dups['frame'] == curr_first_frame].index[0]
        end_ind = no_dups[no_dups['frame'] == curr_final_frame].index[-1]

        trunc_no_dups = no_dups.truncate(before=start_ind, after=end_ind)
        gen_range = range(curr_first_frame, curr_final_frame)
        for frame_num in gen_range:
            curr_time = trunc_no_dups.loc[(trunc_no_dups['frame'] == frame_num) & (trunc_no_dups['row'] == 0) & (
                        trunc_no_dups['col'] == 0), 'time'].iloc[0]

            if prev_time != curr_time:
                hour, minute, second = prev_time.split("_")

                if hour != prev_hour:
                    if int(hour) == 0:
                        day_count += 1

                if int(hour) == -1:
                    hour, minute, second = curr_time.split("_")

                overalltime = f"{year}_{month}_{day_count}_{hour}_{minute}_{second}"

                second = int(second)
                if second - prev_second > 1:
                    for drop in range(prev_second + 1, second):
                        dropped_time = datetime.datetime(int(year), int(month), int(day_count), int(hour), int(minute),
                                                         int(drop), 0)
                        dropped_seconds.append(dropped_time)

                count = 0
                prev_hour = hour
                prev_second = second

            year, month, day, hour, minute, second = overalltime.split("_")
            date = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))
            ms_date = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second),
                                        int(frame_num - frame_start) % 100000)
            timestamp_data_array.append(date)
            ms_timestamp_data_array.append(ms_date)
            timestamp_data_dict[date] = int(frame_num - frame_start)
            count += 1
            prev_time = curr_time

            if frame_num == int(frame_start):
                start_time = date

            end_time = date

    timests = (ms_timestamp_data_array, timestamp_data_dict, dropped_seconds, end_time, start_time)
    with open(centroidfile[:-4] + ".p", "wb") as fp:
        pickle.dump(timests, fp)

    return timests


def get_rois_from_csv(cell_filename):
    roi_bboxs = {}
    with open(cell_filename, 'rb') as f:
        rois = pickle.load(f)
        i = 1
        for roi in rois:
            xs = [point[0] for point in roi]
            ys = [point[1] for point in roi]
            roi_bboxs[i] = [min(xs), min(ys), max(xs), max(ys)]
            i += 1
    return roi_bboxs

# Start here
def loading_procedures():
    if roisfile or longmovie or social:
        rois_dict = get_rois_from_csv(roisfile)

    startdate = "6/3/2025"
    start_time =  datetime.datetime(2025, 6, 3, 18, 2, 5)
    end_time = datetime.datetime(2025, 6, 6, 13, 36, 15)
    section_events = load_sections_file(startdate, end_time, start_time)
    
    for e in section_events:
        print('e', e.events, e.name, e.starttime, e.endtime)
        for ev in e.events.keys():
            print(ev, e[ev])
    hs_dpix, hs_pos = load_highspeed_data(startdate, section_events)

    df = pd.read_csv(centroidfile)
    df.head()
    no_dups = df.drop_duplicates()

    if longmovie:
        firstdpix = open(dpixfile, 'r')
        dp_data_list = []
        dlines = firstdpix.readlines()
        for dline in dlines:
            dp_data_list.append(int(dline))
        dp_data_array = np.array(dp_data_list)
        dp_data_array = dp_data_array.reshape(
            dp_data_array.size // num_of_wells, num_of_wells)
        cenfile = open(centroidfile, 'r')
        cen_data_list = []
        clines = cenfile.readlines()
        for cline in clines:
            cen_data_list.append(int(cline))
        cen_data_array = np.array(cen_data_list)
        cen_data_array = cen_data_array.reshape(
            cen_data_array.size // (num_of_wells*2), (num_of_wells*2))
    else:
        cen_data_array, dp_data_array = get_dists_and_posns(no_dups)

        dp_data_array.resize(dp_data_array.size //
                             num_of_wells, num_of_wells)

        cen_data_array = cen_data_array.reshape(
            cen_data_array.size // (num_of_wells*2), (num_of_wells*2))

    centroid_pickle = centroidfile[:-4]+".p"
    if os.path.exists(centroid_pickle):
       with open(centroid_pickle, "rb") as fp:
           tuple_timestamps = pickle.load(fp)
    else:
        tuple_timestamps = get_timestamps_from_csv(no_dups)

    print("Done loading timestamp file")

    # just setting to zero to make it easier to ignore
    cen_data_array[cen_data_array == 65535] = 0
    # This is because they are all values of -16 and the initial type of the array is unsigned int, but it should be clear that it means the fish hasn't moved yet
    # converting them to zero for now so that it makes it easier to deal with the downstream max/min tests
    tuple_rho_theta = convert_to_polar(cen_data_array)
    print("Done converting to polar coordinates")

    # This is when the run started, calculated from the timestamp file and it is useful for putting event data on correct day
    # Day "0" in the sections file corresponds to this start date from the very beginning of the timestamp file
    startdate = str(tuple_timestamps[0][0].month) + "/" + str(
        tuple_timestamps[0][0].day) + "/" + str(tuple_timestamps[0][0].year)
    endDT = tuple_timestamps[3]
    print('startdate', startdate, 'endDT', endDT)
    startDT = tuple_timestamps[4]
    print('startDT', startDT)

#    global_tuple_events = load_event_data(startdate, endDT, startDT)
    print("Done loading events")

    fish_list = generate_fish_objects(dp_data_array, tuple_rho_theta[0], tuple_rho_theta[1], tuple_rho_theta[2],
                                      tuple_rho_theta[3], hs_dpix, hs_pos, rois_dict)
#    print(fish_list)
#    print(global_tuple_events[2])

    return fish_list, tuple_timestamps[0], tuple_timestamps[1], tuple_timestamps[2], section_events

def initialize_args():
    print("Initializing arguments")

if __name__ == "__main__":
    rois_dict = get_rois_from_csv(roisfile)
    print(rois_dict)


