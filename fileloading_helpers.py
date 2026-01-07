#!/usr/bin/python

import glob
import numpy as np
import datetime
from os import path
from EventSection import EventSection  # event section object

# Convert cartesian coordinates to polar
def cart2pol(x, y):
    rho = np.sqrt(x ** 2 + y ** 2)
    theta = np.arctan2(y, x)
    return rho, theta


# Loading of the centroid position data for the high-speed movies
def load_movie_motion_pos(hs_pos, thistime, fn, num_of_wells):
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
        hscen_data_array.size // (num_of_wells * 2), (num_of_wells * 2))

    hs_pos[thistime] = hscen_data_array


# Loading of the delta pixel motion data for the high-speed movies
def load_movie_motion(hs_dpix, thistime, moviename, num_of_wells):
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
def load_sections_file(startdate, end_time, start_time, sectionsfile):
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
            print(
                "The format of the input of your sections files is wrong. Every line should be as follows: habituation_daytaps=1_15:25:30-1_16:59:00")
            print("Exiting")
            exit()
        # Makes sure that the last section is not after the true end of the data
        eventsection.check_endtime(end_time)
        # Makes sure that the last section is not after the true end of the data
        eventsection.check_starttime(start_time)
        events.append(eventsection)

    return events


def load_highspeed_data(startdate, events, msecperframe, movieprefix, events_fn, num_of_wells):
    f = open(events_fn, 'r')
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
        for eventsection in events:
            if eventsection.starttime <= thistime <= eventsection.endtime:
                eventsection.add_event(TAPES[2], TAPES[3].strip(';'), thistime, msecperframe)

        if int(TAPES[2]) != 0 and hscounter < len(highspeedlist):
            # Load in only the high-speed movie data marked '1' and with the high-speed movie prefix
            # This would need to be modified if a labview code other than '1' was used to make high-speed movies that should analyzed in this way
            # Would also need to modify code in processmotiondata.py, specifically the CalculateEventProperties function relies on this "1" to trigger analysis
            if int(TAPES[2]) == 1:
                load_movie_motion(hs_dpix, thistime, highspeedlist[hscounter], num_of_wells)
                load_movie_motion_pos(hs_pos, thistime, highspeedlist[hscounter], num_of_wells)
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
def max_min(cen_data_array, num_of_wells):
    maxxys = []
    minxys = []
    for n in range(0, num_of_wells * 2, 2):
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
            maxrealy = np.amax(cen_data_array[:, n + 1])
            minrealy = np.amin(
                cen_data_array[:, n + 1][np.nonzero(cen_data_array[:, n + 1])])
            maxxys.append(maxrealx)
            maxxys.append(maxrealy)
            minxys.append(minrealx)
            minxys.append(minrealy)
    maxxysnp = np.array(maxxys)
    minxysnp = np.array(minxys)
    return maxxysnp, minxysnp


# Polar coordinates are essential for easy calculation of well edge/center preferences
def convert_to_polar(cen_data_array, num_of_wells, social):
    (maxxysnp, minxysnp) = max_min(cen_data_array, num_of_wells)
    midcoords = (maxxysnp + minxysnp) / 2
    midcoords = midcoords.astype(np.int16)
    cen_data_array = cen_data_array.astype(np.int16)
    # just setting to very low value to make it easier to skip later
    cen_data_array[cen_data_array == 0] = -10000
    # subtract middle coordinate to get everything centered about 0
    zerodcoords = np.zeros(np.shape(cen_data_array))
    for i in range(0, num_of_wells * 2):
        zerodcoords[:, i] = cen_data_array[:, i] - midcoords[i]
    zerodcoords[zerodcoords < -5000] = 0
    # zerodcoords currently contains negative numbers, which I think mean that the fish hasn't moved yet
    thetadata = np.zeros((len(cen_data_array), num_of_wells))
    rhodata = np.zeros((len(cen_data_array), num_of_wells))
    xzerod = np.zeros((len(cen_data_array), num_of_wells))
    yzerod = np.zeros((len(cen_data_array), num_of_wells))
    for i in range(0, num_of_wells):
        (rhodata[:, i], thetadata[:, i]) = cart2pol(
            zerodcoords[:, 2 * i], zerodcoords[:, 2 * i + 1])
        xzerod[:, i] = zerodcoords[:, 2 * i]
        yzerod[:, i] = zerodcoords[:, 2 * i + 1]
        if social:
            # RETURNING X/Y NOT ZEROED!!
            xzerod[:, i] = cen_data_array[:, 2 * i]
            yzerod[:, i] = cen_data_array[:, 2 * i + 1]
    return rhodata, thetadata, xzerod, yzerod
