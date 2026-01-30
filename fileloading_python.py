import pickle
import math
import os
import datetime

import pandas as pd
import numpy as np

from collections import Counter


def distwrapper(a, b, c, d):
    return math.dist((a, b), (c, d))

def get_dpix_and_posns(no_dups, num_of_wells):
    posx = np.array(no_dups['pos_x'])
    posy = np.array(no_dups['pos_y'])
    dpix = np.array(no_dups['dpix'])

    """
    vfunc = np.vectorize(distwrapper)

    curr_x = posx[:-num_of_wells]
    curr_y = posy[:-num_of_wells]

    next_x = posx[num_of_wells:]
    next_y = posy[num_of_wells:]

    final_frame = len(posx)

    dists = []
    for i in range(0, final_frame, 96*100):
        curr_first_frame = i
        if i + 96*100 < final_frame:
            curr_final_frame = i + 96*100
        else:
            curr_final_frame = final_frame

        dists.extend(vfunc(curr_x[curr_first_frame:curr_final_frame], curr_y[curr_first_frame:curr_final_frame], next_x[curr_first_frame:curr_final_frame], next_y[curr_first_frame:curr_final_frame]))#np.concatenate([a[:-1], a[1:]], axis=1))

    print(np.array(dists).shape)
    dists = np.resize(np.array(dists), [int(len(dists)/96), 96])
    """

    a = np.zeros(len(posx)*2)
    a[0::2] = posx
    a[1::2] = posy
    posns_by_frame = a.reshape([int((a.shape[0])/(num_of_wells*2)), (num_of_wells*2)])
    return np.array(posns_by_frame), np.array(dpix)

def sanity_check(df,
                 num_of_cells):  # there should never be any rows with same frame number, and if we know the expected number of cells per frame, it's trivial to check that that case isn't the case
    c = Counter(df['frame'].tolist())

    if len(set(c.values())) > 1:
        print(
            f'There are multiple runs in the file given. Trimming dataframe to the the most recent of all of the runs.')

        for frame, count in c.items():
            if count == num_of_cells:
                return df[df['frame'] == frame].index[0]

    return None


def get_timestamps_from_csv(no_dups, savefn=""):
    frame_nums = list(no_dups['frame'].dropna())
    prev_time = "-1_-1_-1"  # times[0]
    #    print(prev_time)
    prev_hour = 0
    year = 2025
    month = 6
    day = 3

    count = 0
    frame_start = int(frame_nums[0])
    final_frame = int(frame_nums[-1])

    prev_second = 0
    timestamp_data_dict = {}
    timestamp_data_array = []
    ms_timestamp_data_array = []

    dropped_seconds = []

    div = 10000
    for i in range(frame_start, final_frame, div):
        curr_first_frame = i
        if i + div < final_frame:
            curr_final_frame = i + div
            gen_range = range(0, div)
        else:
            curr_final_frame = final_frame
            gen_range = range(0, final_frame - i + 1)

        start_ind = no_dups[no_dups['frame'] == curr_first_frame].index[0]
        end_ind = no_dups[no_dups['frame'] == curr_final_frame].index[-1]

        trunc_no_dups = no_dups.truncate(before=start_ind, after=end_ind)

        times = list(trunc_no_dups.loc[(trunc_no_dups['row'] == 0) & (
                trunc_no_dups['col'] == 0), 'time'])

        frames = list(trunc_no_dups.loc[(trunc_no_dups['row'] == 0) & (
                trunc_no_dups['col'] == 0), 'frame'])

        for r in gen_range:
#            print(r)
            curr_time = times[r]
            frame_num = int(frames[r])
#            print(frame_num, frame_num-frame_start)
            if prev_time != curr_time:
                hour, minute, second = prev_time.split("_")

                if hour != prev_hour:
                    if int(hour) == 0:
                        day += 1

                if int(hour) == -1:
                    hour, minute, second = curr_time.split("_")

                date = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))
                
                second = int(second)
                if second - prev_second > 1:
                    for drop in range(prev_second + 1, second):
                        dropped_time = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute),
                                                         int(drop), 0)
                        dropped_seconds.append(dropped_time)

                allframes = trunc_no_dups.loc[(trunc_no_dups['time'] == curr_time), 'frame'].tolist()
                totalframes = int(allframes[-1]) - int(allframes[0])

                mstimer = 0
                count = 0
                prev_hour = hour
                prev_second = second


            if mstimer > 0:
                if mstimer >= totalframes:
                    totalframes += 1
                mstime = int((mstimer/(totalframes+1))*1000000)
                   
            else:
                mstime = 0

            ms_date = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), mstime)

            timestamp_data_array.append(date)
            ms_timestamp_data_array.append(ms_date)

            timestamp_data_dict[date] = int(frame_num - frame_start)
            mstimer += 1
            prev_time = curr_time

            if frame_num == int(frame_start):
                start_time = date

            end_time = date

    timests = (ms_timestamp_data_array, timestamp_data_dict, dropped_seconds, end_time, start_time)

    if len(savefn) > 0:
        with open(savefn, "wb") as fp:
            pickle.dump(timests, fp)

    return timests


def get_rois_from_csv(cell_filename):
    roi_bboxs = {}
    with open(cell_filename, 'rb') as f:
        rois = pickle.load(f)
        i = 1
        for roi in rois:
            print(roi)
            xs = [point[0] for point in roi]
            ys = [point[1] for point in roi]
            roi_bboxs[i] = [min(xs), min(ys), max(xs), max(ys)]
            print(roi_bboxs[i])
            i += 1
    i -= 1
    return roi_bboxs, i


def load_python_data(centroid_file, rois_file):
    df = pd.read_csv(centroid_file) #, header=5183042)
#    df.columns = ['time', 'frame', 'row', 'col', 'pos_x', 'pos_y', 'dpix']

#    df = df.drop_duplicates()
    print(df)
    print('df loaded')

    rois_dict, num_of_wells = get_rois_from_csv(rois_file)
    print('rois loaded')

    val = sanity_check(df, num_of_wells)
    if val:
        df = df.truncate(before=val)

    print('sanity check completed')

    timestamp_pickle = centroid_file[:-4] + ".p"
    if os.path.exists(timestamp_pickle):
        print("opening existing timestamp file")
        with open(timestamp_pickle, "rb") as fp:
            tuple_timestamps = pickle.load(fp)
    else:
        tuple_timestamps = get_timestamps_from_csv(df, timestamp_pickle)

    print('timestamps loaded')
    
    cen_data_array, dp_data_array = get_dpix_and_posns(df, num_of_wells)
    print('posns loaded')

    return rois_dict, num_of_wells, cen_data_array, dp_data_array, tuple_timestamps

