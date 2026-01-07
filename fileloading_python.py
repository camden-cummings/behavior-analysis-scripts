import pickle
import math
import os
import datetime

import pandas as pd
import numpy as np


def get_dists_and_posns(no_dups, num_of_wells):
    posns_by_frame = []
    dists_by_frame = []

    indic = [i for i in range(0, int(len(no_dups)), num_of_wells)]
    frame_arr = no_dups.iloc[indic[0]:indic[1]]
    xs = list(frame_arr['pos_x'])
    ys = list(frame_arr['pos_y'])

    for count in range(len(indic) - 2):
        next_frame_arr = no_dups.iloc[indic[count + 1]:indic[count + 2]]

        next_xs = list(next_frame_arr['pos_x'])
        next_ys = list(next_frame_arr['pos_y'])

        frame_xys = np.array(list(zip(xs, ys))).flatten()
        dists = [math.dist((xs[i], ys[i]), (next_xs[i], next_ys[i])) for i in range(len(xs))]

        posns_by_frame.append(frame_xys)
        dists_by_frame.append(dists)

        xs = next_xs
        ys = next_ys

    return np.array(posns_by_frame), np.array(dists_by_frame)

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

    for i in range(frame_start, final_frame, 1000):
        curr_first_frame = i
        if i + 1000 < final_frame:
            curr_final_frame = i + 1000
            gen_range = range(0,1000)
        else:
            curr_final_frame = final_frame
            gen_range = range(0,final_frame-i+1)


        start_ind = no_dups[no_dups['frame'] == curr_first_frame].index[0]
        end_ind = no_dups[no_dups['frame'] == curr_final_frame].index[-1]

        trunc_no_dups = no_dups.truncate(before=start_ind, after=end_ind)

        times = list(trunc_no_dups.loc[(trunc_no_dups['row'] == 0) & (
                trunc_no_dups['col'] == 0), 'time'])

        frames = list(trunc_no_dups.loc[(trunc_no_dups['row'] == 0) & (
                trunc_no_dups['col'] == 0), 'frame'])

        for r in gen_range:
            curr_time = times[r]
            frame_num = int(frames[r])
            if prev_time != curr_time:
                hour, minute, second = prev_time.split("_")

                if hour != prev_hour:
                    if int(hour) == 0:
                        day += 1

                if int(hour) == -1:
                    hour, minute, second = curr_time.split("_")

                date = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))
                ms_date = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second),
                                            int(frame_num - frame_start) % 100000)

                second = int(second)
                if second - prev_second > 1:
                    for drop in range(prev_second + 1, second):
                        dropped_time = datetime.datetime(int(year), int(month), int(day), int(hour), int(minute),
                                                         int(drop), 0)
                        dropped_seconds.append(dropped_time)

                count = 0
                prev_hour = hour
                prev_second = second

            timestamp_data_array.append(date)
            ms_timestamp_data_array.append(ms_date)
            timestamp_data_dict[date] = int(frame_num - frame_start)
            count += 1
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
            xs = [point[0] for point in roi]
            ys = [point[1] for point in roi]
            roi_bboxs[i] = [min(xs), min(ys), max(xs), max(ys)]
            i += 1
    return roi_bboxs

def load_python_data(centroid_file, rois_file, num_of_wells):
    df = pd.read_csv(centroid_file)
    df.head()
    no_dups = df.drop_duplicates()

    rois_dict = get_rois_from_csv(rois_file)

    cen_data_array, dp_data_array = get_dists_and_posns(no_dups, num_of_wells)

    timestamp_pickle = centroid_file[:-4] + ".p"
    if os.path.exists(timestamp_pickle):
        with open(timestamp_pickle, "rb") as fp:
            tuple_timestamps = pickle.load(fp)
    else:
        tuple_timestamps = get_timestamps_from_csv(no_dups, timestamp_pickle)

    return rois_dict, cen_data_array, dp_data_array, tuple_timestamps