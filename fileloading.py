#!/usr/bin/python

import os
import glob
import argparse
import numpy as np
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
parser.add_argument('-graphmulti', type=str, action="store", dest="dirlist")
parser.add_argument('-j', type=str, action="store",
                    dest="graphparameters", default="PlotParameters")

parser.add_argument('-e', type=str, action="store", dest="eventsfile")
parser.add_argument('-c', type=str, action="store", dest="centroidfile")
parser.add_argument('-d', type=str, action="store", dest="dpixfile")
parser.add_argument('-m', type=str, action="store",
                    dest="movieprefix", default="")
parser.add_argument('-g', type=str, action="store", dest="genotypefile")
parser.add_argument('-s', type=str, action="store",
                    dest="sectionsfile", default="sectionsfile")
parser.add_argument('-n', type=int, action="store",
                    dest="numberofwells", default=96)

args = parser.parse_args()
longmoviename = args.longmoviename
longmovie = False
if (longmoviename != "nomovie"):
    longmovie = True
roisfile = args.roisfile
oldtracking = args.oldtracking
graphonly = args.graphonly
xyhm = args.xyhm
social = args.social
graphmulti = args.dirlist
if ((not graphonly) and (not graphmulti)):
    eventsfile = args.eventsfile
    centroidfile = args.centroidfile
    dpixfile = args.dpixfile
    movieprefix = args.movieprefix
    genotypefile = args.genotypefile
    sectionsfile = args.sectionsfile
elif (graphmulti):
    graphmulti = list(map(str, args.dirlist.split(',')))

numberofwells = args.numberofwells

# Convert cartesian coordinates to polar
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return (rho, theta)

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
        hscen_data_array.size // (numberofwells*2), (numberofwells*2))
#    print(hscen_data_array)
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
#        print(line)
    hsdp_data_array = np.array(hsdp_data_list)
#    print(hsdp_data_array)
    hsdp_data_array = hsdp_data_array.reshape(
        hsdp_data_array.size // numberofwells, numberofwells)
    hs_dpix[thistime] = hsdp_data_array

# Loading in the high-speed movies and the sections we will analyze (sectionsfile)
def load_event_data(startdate, endDT, startDT):
    # Load in the original events file
    # Tab separation vs spacing is KEY in the events file (LabView will fail if not right, and so will this code)
    lastAMorPM0 = None
    lastAMorPMcounter0 = 0
    hs_dpix = {}
    hs_pos = {}
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
        eventsection.check_endtime(endDT)
        # Makes sure that the last section is not after the true end of the data
        eventsection.check_starttime(startDT)
        events.append(eventsection)
    #print(eventsection)
    f = open(eventsfile, 'r')
    lines = f.readlines()
    #print(lines)
    hscounter = 0
    # 1:06:24\tPM\t0\t2\n'
    # Getting list of all the ID numbers associated with just the high-speed short movies
    highspeedlist = []
    for file in glob.glob(movieprefix + '*motion2'):
        highspeedlist.append(file)
#    highspeedlist.sort()
#    print("hslist", highspeedlist
    highspeedlist.sort(key=lambda name : int(name.split('/')[3].split('-')[1].split('.')[0]))
    print(highspeedlist)
    for line in lines:
        TAPES = line.split('	')
        # Need to be splitting on tab only for the new format of the input file with the spaces in last column
        # Converting date to the format that we were using for the
        # 6/18/2016_12:59:34_PM
        dateplustime = startdate + "_" + \
            TAPES[0][0:len(TAPES[0])] + "_" + TAPES[1]
        thistime = datetime.datetime.strptime(
            dateplustime, '%m/%d/%Y_%I:%M:%S_%p')
        # thistime = faststrptime(dateplustime)
        thisAMorPM0 = TAPES[1]
        if lastAMorPM0 != None:
            # If we transition either AM->PM or PM->AM
            # The only time we will ever have issues with this is if the events file doesn't correspond well with beginning timestamp input
            # Or there could be issues with start after midnight
            # These situations should not happen
            # Basically things become a mess if the program does not read the events file in as it should have and it skips starting at beginning
            if thisAMorPM0 != lastAMorPM0:
                if thisAMorPM0 == "AM":
                    lastAMorPMcounter0 = lastAMorPMcounter0 + 1
        lastAMorPM0 = thisAMorPM0
        thistime = thistime + datetime.timedelta(days=(lastAMorPMcounter0))
        for eventsection in events:
            if eventsection.starttime <= thistime <= eventsection.endtime:
                # Debug
                # print(TAPES[2],TAPES[3].strip(';'),  thistime, e.starttime, e.endtime)
                eventsection.add_event(TAPES[2], TAPES[3].strip(';'), thistime)
        lasttime = thistime
#        highspeedlist.sort()
#        highspeedlist.sort(key=lambda name : int(name.split('/')[3].split('-')[1].split('.')[0]))
        #print(TAPES)
#        print(hscounter, highspeedlist[hscounter])
        if int(TAPES[2]) != 0 and hscounter < len(highspeedlist):
            # Load in only the high-speed movie data marked '1' and with the high-speed movie prefix
            # This would need to be modified if a labview code other than '1' was used to make high-speed movies that should analyzed in this way
            # Would also need to modify code in processmotiondata.py, specifically the CalculateEventProperties function relies on this "1" to trigger analysis
            if int(TAPES[2]) == 1:
                print(hscounter, int(highspeedlist[hscounter].split("-")[1].split(".")[0]), highspeedlist[hscounter], TAPES)
#                print("1-------------------")
                load_movie_motion(hs_dpix, thistime, highspeedlist[hscounter])
#                print(hs_dpix)
#                print("2--------------------")
                load_movie_motion_pos(
                    hs_pos, thistime, highspeedlist[hscounter])
                hscounter = hscounter + 1

#    print(highspeedlist)
#    print(hscounter, len(highspeedlist))
    if hscounter != len(highspeedlist):
        print("ERROR, the number of high-speed movies does not correspond to the number expected based on the event file")
        print("Please make sure that the event file is accurate and run was fully completed. Modify events file to reflect reality of data output if needed")
        print("Exiting")
#        exit()

    # Debug to make sure eveything is added correctly
#    for e in events:
#        print("e: ", e.name, e.type, e.events)
#    print(hs_dpix)
#    print(hs_pos)
#    print(events)
    return (hs_dpix, hs_pos, events)

# Find max and min value for each fish in order to identify well edges
# This is actually not really the well edges, but the edges of the fish's radius of movement
# I prefer this approach to reading in the originally designated ROI, just in case ROI isn't accurate (ie, includes extra plastic of well edge)
def max_min(cen_data_array):
    maxxys = []
    minxys = []
    for n in range(0, numberofwells*2, 2):
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
    return (maxxysnp, minxysnp)

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
    for i in range(0, numberofwells*2):
        zerodcoords[:, i] = cen_data_array[:, i] - midcoords[i]
    zerodcoords[zerodcoords < -5000] = 0
    # zerodcoords currently contains negative numbers, which I think mean that the fish hasn't moved yet
    thetadata = np.zeros((len(cen_data_array), numberofwells))
    rhodata = np.zeros((len(cen_data_array), numberofwells))
    xzerod = np.zeros((len(cen_data_array), numberofwells))
    yzerod = np.zeros((len(cen_data_array), numberofwells))
    for i in range(0, numberofwells):
        (rhodata[:, i], thetadata[:, i]) = cart2pol(
            zerodcoords[:, 2*i], zerodcoords[:, 2*i+1])
        xzerod[:, i] = zerodcoords[:, 2*i]
        yzerod[:, i] = zerodcoords[:, 2*i+1]
        if social:
            # RETURNING X/Y NOT ZEROED!!
            xzerod[:, i] = cen_data_array[:, 2*i]
            yzerod[:, i] = cen_data_array[:, 2*i+1]
    return (rhodata, thetadata, xzerod, yzerod)

# The Fish object carries all the data around for each fish, including their genotype and ID number
# Later (after processmotiondata.py) the data inside this Fish object is analyzed (bouts counted, binned, responses counted) and the AnalyzedFish object carries that data
# This analysis code only compares two groups: a control group and a test group
# Or it can analyze a single group, but no statistics will be done
def generate_fish_objects(dp_data_array, rho_array, theta_array, x_array, y_array, hs_dpix, hs_pos, genotypefile, rois_dict):
    f = open(genotypefile, 'r')
    lines = f.readlines()
    f.close()
    genotype_list = {}
    for line in lines:
        # The file must use "controlgroup_geno: #,#,#" and "testgroup_geno: #,#,#,#" to emphasize to users that only two are allowed
        genogroup = line.split(':')[0].split('_')[0]
        if (len(line.split(':')[0].split('_')[0]) > 1):
            realgenotype = line.split(':')[0].split('_')[1]
        else:
            if ((line.split(':')[0].split('_')[0]) == "controlgroup"):
                realgenotype = "control"
            elif ((line.split(':')[0].split('_')[0]) == "testgroup"):
                realgenotype = "test"
            else:
                print(
                    "Not using correct nomenclature of controlgroup_genotype: and testgroup_genotype:")
                realgenotype = line.split(':')[0].split('_')[0]
        fishidslist = line.split(':')[1].strip().split(',')
        print(fishidslist)
        inputids = []
        for id in fishidslist:
            print(id)
            # Have to subtract 1 so that we can start from 0
            # Keeping it like this, but then adding 1 back to the ID that is saved
            inputids.append(int(id)-1)
        genotype_list[genogroup + "_" + realgenotype] = inputids
    fish_list = []
    for n in range(0, numberofwells):
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
                    # print("BEFORE ",np.amax(x_array[:,n]), " ", np.amin(x_array[:,n]))
                    lower_bound = rois_dict[n+1][0]
                    upper_bound = rois_dict[n+1][2]
                    roimask = (x_array[:, n] < lower_bound) | (
                        x_array[:, n] > upper_bound)
                    # roimask = ((x_array[:,n] < rois_dict[n+1][0]) | (x_array[:,n] > rois_dict[n+1][2]))
                    # print(x_array[:,n])
                    x_array[:, n][roimask] = np.nan
                    # print(roimask)
                    # print(x_array[:,n])
                    # print("AFTER ",np.nanmax(x_array[:,n]), " ", np.nanmin(x_array[:,n]))
                    # print(x_array[:,n])
                    if rois_dict[n+1][0] < 400:
                        x_array[:, n] = x_array[:, n] - \
                            np.nanmax(x_array[:, n])
                        x_array[:, n] = abs(-1 * x_array[:, n])
                    else:
                        x_array[:, n] = x_array[:, n] - \
                            abs(np.nanmin(x_array[:, n]))
                    # print(x_array[:,n])
                    # print("EARLYMAXMIN",np.nanmin(x_array[:,n]), np.nanmax(x_array[:,n]))
                    # print("EARLY ROIS ", rois_dict[n+1])
                newfish = Fish(n + 1, x.split('_')[0], x.split('_')[1], dp_data_array[:, n], rho_array[:, n],
                               theta_array[:, n], x_array[:, n], y_array[:, n], split_hs_dpix, split_hs_pos_x, split_hs_pos_y)
                if xyhm:
                    nan_mask = np.logical_not(np.isnan(x_array[:, n]))
                    heatmap, xedges, yedges = np.histogram2d(
                        x_array[:, n][nan_mask], y_array[:, n][nan_mask])#, bins=(35, 5))
                    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                    fig, ax1 = plt.subplots()
                    im = ax1.imshow(
                        heatmap.T, interpolation='lanczos', cmap='jet')
                    # im = ax1.pcolormesh(xedges,yedges,heatmap.T,cmap='jet',shading='auto')
                    plt.savefig("fish_" + str(n+1) + "_XYheatmap.png",
                                transparent=True, format="png")
                    plt.close()
                if (roisfile or longmovie):
                    # print(rois_dict[n+1])
                    # print("EARLY ROIS b", rois_dict[n+1])
                    if social:
                        rois_dict[n+1][2] = rois_dict[n+1][2] - \
                            rois_dict[n+1][0]
                        rois_dict[n+1][0] = 0
                        # print("EARLY ROIS c", rois_dict[n+1])
                        newfish.add_rois(rois_dict[n+1])
                    else:
 #                       print(n+1)
                        newfish.add_rois(rois_dict[n+1])
                fish_list.append(newfish)
    # if xyhm:
    # print("Done with making xy heatmaps. Exiting without further analysis.")
    # exit()
    return fish_list


def get_dists_and_posns(numberofwells, preprocessedfile):
    df = pd.read_csv(preprocessedfile)
    df.head()
    no_dups = df.drop_duplicates()

    posns_by_frame = []
    dists_by_frame = []
    indic = [i for i in range(0, int(len(no_dups)), numberofwells)]
    for count in range(len(indic)-2):
        frame = no_dups.iloc[indic[count]:indic[count+1]]
        next_frame = no_dups.iloc[indic[count+1]:indic[count+2]]
#        print(frame, next_frame)

        xs = list(frame['pos_x'])
        ys = list(frame['pos_y'])
        next_xs = list(next_frame['pos_x'])
        next_ys = list(next_frame['pos_y'])

        frame_xys = []
        dists = []
        for i in range(len(xs)):
            frame_xys.append(xs[i])
            frame_xys.append(ys[i])

            dists.append(math.dist((xs[i], ys[i]), (next_xs[i], next_ys[i])))
#        print(dists)
        posns_by_frame.append(frame_xys)
        dists_by_frame.append(dists)
    #print(posns_by_frame)
    return np.array(posns_by_frame), np.array(dists_by_frame)


def get_timestamps_from_csv(filename):
    df = pd.read_csv(filename)
    df.head()
    no_dups = df.drop_duplicates()
#    print(df)
#    print(no_dups)
    times = list(df['time'])
    frame_nums = list(no_dups['frame'].dropna())
    prev_time = "-1_-1_-1" #times[0]
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

    frame_nums = set(frame_nums)
    prev_second = 0
    timestamp_data_dict = {}
    timestamp_data_array = []
    ms_timestamp_data_array = []

    dropped_seconds = []
#    print([i for i in range(frame_start, final_frame, 1000)][:-1])
    for i in range(frame_start, final_frame, 1000):
        curr_first_frame = i
        if i+999 < final_frame:
            curr_final_frame = i+999
        else:
            curr_final_frame = final_frame
#        curr_final_frame = i+999
        start_ind = no_dups[no_dups['frame'] == curr_first_frame].index[0]
        end_ind = no_dups[no_dups['frame'] == curr_final_frame].index[-1]

        trunc_no_dups = no_dups.truncate(before=start_ind, after=end_ind)
#        print(trunc_no_dups)
        gen_range = range(curr_first_frame, curr_final_frame)
        for frame_num in gen_range:
            #print(frame_num)
#            print(trunc_no_dups.loc[(trunc_no_dups['frame'] == frame_num) & (trunc_no_dups['row'] == 0) & (trunc_no_dups['col'] == 0), 'time'])
            curr_time = trunc_no_dups.loc[(trunc_no_dups['frame'] == frame_num) & (trunc_no_dups['row'] == 0) & (trunc_no_dups['col'] == 0), 'time'].iloc[0]
#            print(prev_time, curr_time)
            if prev_time != curr_time:
                hour, minute, second = prev_time.split("_")
#                print(hour, minute, second)
                if hour != prev_hour:
                    #print(day_count, hour, minute, second, frame_num)
                    if int(hour) == 0:
                        day_count += 1
                    print(day_count, hour, minute, second, frame_num)
                if int(hour) == -1:
                    print("begin here", curr_time)
                    #curr_time = prev_time
                    hour, minute, second = curr_time.split("_")

                    # indices.append(reduced_count)

                # reduced_times.append(f"{year}_{month}_{day}_{hour}_{minute}_{second}_{count // numberofwells}")
                overalltime = f"{year}_{month}_{day_count}_{hour}_{minute}_{second}"

                second = int(second)
                #print(second, prev_second)
                if second - prev_second > 1:
                    for drop in range(prev_second+1, second):
                        dropped_time = datetime.datetime(int(year), int(month), int(day_count), int(hour), int(minute), int(drop), 0)
                        dropped_seconds.append(dropped_time)

                count = 0
                prev_hour = hour
                prev_second = second

            year, month, day, hour, minute, second = overalltime.split("_")
            date = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))
#            print(date)
            ms_date = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), int(frame_num-frame_start)%100000)
#            print(ms_date)
            timestamp_data_array.append(date)
            ms_timestamp_data_array.append(ms_date)
            timestamp_data_dict[date] = int(frame_num-frame_start)
            count += 1
            prev_time = curr_time

            if frame_num == int(frame_start):
                start_time = date
                print("start time set")
                print(frame_num, frame_start, start_time) 
            end_time = date
#    print(timestamp_data_dict[datetime.datetime(2025, 6, 3, 9, 29, 49, 3848)])
    print("start", start_time)
    print("end", end_time)
    print("len", len(ms_timestamp_data_array))
    print("lendict", len(timestamp_data_dict))
    
    timests = (ms_timestamp_data_array, timestamp_data_dict, dropped_seconds, end_time, start_time)
    with open(centroidfile[:-4]+".p", "wb") as fp:
        pickle.dump(timests, fp)

    return timests

def get_rois_from_csv(cell_filename):
    print(cell_filename)
    dict = {}
    with open(cell_filename, 'rb') as f:
        rois = pickle.load(f)
        i = 1
        for roi in rois:
            xs = [point[0] for point in roi]
            ys = [point[1] for point in roi]
            dict[i] = [min(xs), min(ys), max(xs), max(ys)]
            i += 1
    return dict

# Start here
def loading_procedures():
    rois_dict = {}
#    if (roisfile or longmovie or social):
#        load_rois(rois_dict)

    rois_dict = get_rois_from_csv(roisfile)
    print(rois_dict)

#    tuple_timestamps = load_timestamp_file()
    startdate = "6/3/2025"
    startDT =  datetime.datetime(2025, 6, 3, 18, 2, 5)
    endDT = datetime.datetime(2025, 6, 5, 13, 36, 15)
    global_tuple_events = load_event_data(startdate, endDT, startDT)
    print("0", global_tuple_events[0])
    print("1", global_tuple_events[1])
    print("2", global_tuple_events[2])

    if (longmovie):
        firstdpix = open(dpixfile, 'r')
        dp_data_list = []
        dlines = firstdpix.readlines()
        for dline in dlines:
            dp_data_list.append(int(dline))
        dp_data_array = np.array(dp_data_list)
        dp_data_array = dp_data_array.reshape(
            dp_data_array.size // numberofwells, numberofwells)
        cenfile = open(centroidfile, 'r')
        cen_data_list = []
        clines = cenfile.readlines()
        for cline in clines:
            cen_data_list.append(int(cline))
        cen_data_array = np.array(cen_data_list)
        cen_data_array = cen_data_array.reshape(
            cen_data_array.size // (numberofwells*2), (numberofwells*2))
    else:
        cen_data_array, dp_data_array = get_dists_and_posns(numberofwells, centroidfile)

        dp_data_array.resize(dp_data_array.size //
                             numberofwells, numberofwells)

        cen_data_array = cen_data_array.reshape(
            cen_data_array.size // (numberofwells*2), (numberofwells*2))
    
    print(cen_data_array.shape)
    print(dp_data_array.shape)
    
    centroid_pickle = centroidfile[:-4]+".p"
    if os.path.exists(centroid_pickle):
       with open(centroid_pickle, "rb") as fp:
           tuple_timestamps = pickle.load(fp)
    else:
       tuple_timestamps = get_timestamps_from_csv(centroidfile)

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
#    endDT = tuple_timestamps[3]
#    print(startdate, endDT)
#    startDT = tuple_timestamps[4]
#    print(startDT)

#    global_tuple_events = load_event_data(startdate, endDT, startDT)
    print("Done loading events")

    fish_list = generate_fish_objects(dp_data_array, tuple_rho_theta[0], tuple_rho_theta[1], tuple_rho_theta[2],
                                      tuple_rho_theta[3], global_tuple_events[0], global_tuple_events[1], genotypefile, rois_dict)
#    print(fish_list)
#    print(global_tuple_events[2])
    return (fish_list, tuple_timestamps[0], tuple_timestamps[1], tuple_timestamps[2], global_tuple_events[2])



def initialize_args():
    print("Initializing arguments")
