import datetime
import argparse
import numpy as np

# This dictionary is only used if -oldtracking is on, for old datasets before ROIs were saved by LabVIEW. Irrelevant for new users.
well_conversion = {0: 0, 8: 1, 16: 2, 24: 3, 32: 4, 40: 5, 48: 6, 56: 7, 64: 8, 72: 9, 80: 10, 88: 11, 1: 12, 9: 13, 17: 14, 25: 15, 33: 16, 41: 17, 49: 18, 57: 19, 65: 20, 73: 21, 81: 22, 89: 23, 2: 24, 10: 25, 18: 26, 26: 27, 34: 28, 42: 29, 50: 30, 58: 31, 66: 32, 74: 33, 82: 34, 90: 35, 3: 36, 11: 37, 19: 38, 27: 39, 35: 40, 43: 41, 51: 42, 59: 43, 67: 44, 75: 45, 83: 46, 91: 47,
                   4: 48, 12: 49, 20: 50, 28: 51, 36: 52, 44: 53, 52: 54, 60: 55, 68: 56, 76: 57, 84: 58, 92: 59, 5: 60, 13: 61, 21: 62, 29: 63, 37: 64, 45: 65, 53: 66, 61: 67, 69: 68, 77: 69, 85: 70, 93: 71, 6: 72, 14: 73, 22: 74, 30: 75, 38: 76, 46: 77, 54: 78, 62: 79, 70: 80, 78: 81, 86: 82, 94: 83, 7: 84, 15: 85, 23: 86, 31: 87, 39: 88, 47: 89, 55: 90, 63: 91, 71: 92, 79: 93, 87: 94, 95: 95}


# A fast function for reading time data
# This depends on have the below format as output, which it is from our LabView program
def faststrptime(val):
    # Example
    # 6/18/201612:59:34 PM
    splits1 = val.split("/")
    splits2 = splits1[2].split(":")
    return datetime.datetime(
        int(splits1[2][0:4]),  # %Year
        int(splits1[0]),  # %Month
        int(splits1[1]),  # %Day
        int(splits2[0][4:len(splits2[0])]),  # %Hour
        int(splits2[1]),  # %Minute
        int(splits2[2][0:2]),  # %Second
    )

# The milliseconds are estimated based on the number of frames per that second
# It is clearly not accurate if a chunk of frames were lost in just the middle, for example
# However it is useful for downstream analyses to have milliseconds, in particular if one has to analyze responses in the slow speed data
def convert_to_ms_time(timestamp_data_array, timestamp_data_dict):
    mstimestamp_data_array = timestamp_data_array
    mseccounter = 0
    for position in range(0, len(timestamp_data_array)):
        if (position + 1) != len(timestamp_data_array):
            if timestamp_data_array[position] == timestamp_data_array[position + 1]:
                mseccounter = mseccounter + 1
            else:
                startpos = position - mseccounter
                msec = 1000.0 / ((position - startpos) + 1)
                for i in range(startpos, position + 1):
                    newsec = timestamp_data_array[i] + \
                        datetime.timedelta(milliseconds=msec*(i-startpos))
                    mstimestamp_data_array[i] = newsec
                mseccounter = 0
        else:
            mseccounter = mseccounter + 1
            startpos = position - mseccounter
            msec = 1000.0 / ((position - startpos) + 1)
            for i2 in range(startpos, position + 1):
                newsec = timestamp_data_array[i2] + \
                    datetime.timedelta(milliseconds=msec*(i2-startpos))
                mstimestamp_data_array[i2] = newsec
            mseccounter = 0
            break
    return mstimestamp_data_array

def load_dpix(dpixfile, numberofwells):
    with open(dpixfile, 'rb') as fid:
        dp_data_array = np.fromfile(fid, dtype='>u2')
    dp_data_array = dp_data_array.reshape(dp_data_array.size // numberofwells, numberofwells)
    return dp_data_array

def load_centroids(centroidfile, numberofwells):
    with open(centroidfile, 'rb') as fid:
        cen_data_array = np.fromfile(fid, '>u2')
    cen_data_array = cen_data_array.reshape(cen_data_array.size // (numberofwells * 2), (numberofwells * 2))
    return cen_data_array

# Load the roi data, for analyzing long movies
def load_rois(roisfile, social, roi_dict):
    f = open(roisfile, 'r')
    lines = f.readlines()
    i = 1
    for line in lines:
        try:
            print(int(line.split(' ')[0]))
        except ValueError:
            continue
        minx = int(line.split(' ')[0])
        miny = int(line.split(' ')[1])
        maxx = int(line.split(' ')[2])
        maxy = int(line.split(' ')[3])
        if social:
            # print("ROIS BEFORE ", minx, miny, maxx, maxy)
            if minx < 400:  # could make this the halfway point of the image or something, but just going with a value I'm sure will work is ok until the code needs to be more flexible.
                # VALUE USED FOR ALMOST ALL ANALYSIS PREVIOUSLY
                maxx = maxx - (0.13 * (maxx - minx))
                # minx = minx + (0.01 * (maxx - minx)) # Also trimming on this side to avoid calls of the fish in the acrylic window in some cases
            else:
                minx = minx + (0.13 * (maxx - minx))
                # maxx = maxx - (0.01 * (maxx - minx))
            # print("ROIS AFTER ", minx, miny, maxx, maxy)
        roi_dict[i] = [minx, miny, maxx, maxy]
        i += 1

# The sections file is in 24 hour time, as is Python code
def load_timestamp_file(tstampfile):
    # Loading in the timestamp file data
    # Used to have to get rid of the ^M character that is between the times, but not with updated Anaconda
    timestamp_data_array = []
    dropped_seconds = []
    f = open(tstampfile, 'r')
    lines = f.readlines()
    f.close()
    for line in lines:
#        print(line)
        timestamp_data_array.append(line.strip())
    n = 0
    timestamp_data_dict = {}
    lasttime = None
    starttime = None
    for t in timestamp_data_array:
        thistime = faststrptime(t)
        # I can't include the AM/PM in the faststrptime loading for speed reasons
        # So input times are assumed to just be AM (other than 12, assumed PM), and this adding/subtracting converts to 24 hour
        thisAMorPM0 = t.split()[len(t.split())-1]
        if thistime.hour == 12:
            if thisAMorPM0 == "AM":
                thistime = thistime + datetime.timedelta(hours=-12)
        elif thisAMorPM0 == "PM":
            thistime = thistime + datetime.timedelta(hours=12)
        timestamp_data_array[n] = thistime
        timestamp_data_dict[thistime] = n
        # Account for situation at beginning, so we don't enter code below and try to subtract
        if n == 0:
            n = n + 1
            lasttime = thistime
            starttime = thistime
            continue
        # This step is important later for the fast slicing of the data
        # Missing seconds need to be accounted for so we know not to expect it
        # Ideally we do not lose that amount of frames, but it can happen
        tdelta1 = thistime - lasttime
        testtime = thistime - \
            datetime.timedelta(seconds=tdelta1.total_seconds())
        if thistime != lasttime:
            timestamp_data_dict[thistime] = n
            if tdelta1.total_seconds() > 1:
                for x in range(0, int(tdelta1.total_seconds()-1)):
                    print("DROPPED A SECOND: ", thistime, lasttime, testtime, testtime + datetime.timedelta(
                        seconds=1), timestamp_data_array[n], n-1, timestamp_data_array[n-1])
                    dropped_seconds.append(
                        testtime + datetime.timedelta(seconds=1))
                    testtime = testtime + datetime.timedelta(seconds=1)
            lasttime = thistime
        n = n + 1
    mstimestamp_data_array = convert_to_ms_time(
        timestamp_data_array, timestamp_data_dict)
#    print(len(timestamp_data_array), len(mstimestamp_data_array), len(timestamp_data_dict))
    return mstimestamp_data_array, timestamp_data_dict, dropped_seconds, thistime, starttime

def load_labview_data(timestamp_file, rois_file, dpix_file, centroid_file, num_of_wells, social, longmovie):
    tuple_timestamps = load_timestamp_file(timestamp_file)
    rois_dict = {}
    load_rois(rois_file, social, rois_dict)

    if longmovie:
        firstdpix = open(dpix_file, 'r')
        dp_data_list = []
        dlines = firstdpix.readlines()
        for dline in dlines:
            dp_data_list.append(int(dline))
        dp_data_array = np.array(dp_data_list)
        dp_data_array = dp_data_array.reshape(
            dp_data_array.size // num_of_wells, num_of_wells)
        cenfile = open(centroid_file, 'r')
        cen_data_list = []
        clines = cenfile.readlines()
        for cline in clines:
            cen_data_list.append(int(cline))
        cen_data_array = np.array(cen_data_list)
        cen_data_array = cen_data_array.reshape(
            cen_data_array.size // (num_of_wells * 2), (num_of_wells * 2))
    else:
        cen_data_array = load_centroids(centroid_file, num_of_wells)
        dp_data_array = load_dpix(dpix_file, num_of_wells)

    print(rois_dict, cen_data_array, dp_data_array)
    return rois_dict, cen_data_array, dp_data_array, tuple_timestamps
