import numpy as np
from ProcessedData import ProcessedData  # ProcessedData object

def socialactivityb(x_array, rois, intervalnumindices, intervalindices, timebintup, prefix=""):
    activitydict = {
        prefix + "socialpreferenceframesb": []
    }
    binnedlist = []
    qpt = rois[2]*0.75  # maxx
    for m in range(0, len(intervalindices)-1, 2):
        startindex = np.searchsorted(intervalnumindices, intervalindices[m])
        # Add one because range below is non-inclusive
        endindex = np.searchsorted(
            intervalnumindices, intervalindices[m+1]) + 1
        shortarray = intervalnumindices[startindex:endindex]

        sframes = 0
        asframes = 0
        tframes = 0
        for m3 in range(0, len(shortarray)-1, 2):
            for m4 in range(shortarray[m3], shortarray[m3+1]+1):
                tframes = tframes + 1
                if x_array[m4] > qpt:
                    sframes = sframes + 1
                if x_array[m4] < qpt:
                    asframes = asframes + 1
        socialpreference = (sframes - asframes) / tframes
        activitydict[prefix +
                     "socialpreferenceframesb"].append(socialpreference)
    for k2, v2 in activitydict.items():
        binnedlist.append(ProcessedData(k2, v2, timebintup))
    return binnedlist

def socialactivity(x_array, rois, intervalnumindices, intervalindices, timebintup, prefix=""):
    activitydict = {
        prefix + "socialpreferenceframes": []
    }
    binnedlist = []
    qpt = rois[2]*0.75  # maxx
    cqpt = rois[2]*0.25
    for m in range(0, len(intervalindices)-1, 2):
        startindex = np.searchsorted(intervalnumindices, intervalindices[m])
        # Add one because range below is non-inclusive
        endindex = np.searchsorted(
            intervalnumindices, intervalindices[m+1]) + 1
        shortarray = intervalnumindices[startindex:endindex]

        sframes = 0
        asframes = 0
        tframes = 0
        for m3 in range(0, len(shortarray)-1, 2):
            for m4 in range(shortarray[m3], shortarray[m3+1]+1):
                tframes = tframes + 1
                if m4 < len(x_array):
                    if x_array[m4] > qpt:
                        sframes = sframes + 1
                    if x_array[m4] < cqpt:
                        asframes = asframes + 1
        socialpreference = (sframes - asframes) / tframes

        activitydict[prefix +
                     "socialpreferenceframes"].append(socialpreference)
    for k2, v2 in activitydict.items():
        binnedlist.append(ProcessedData(k2, v2, timebintup))
    return binnedlist

# Calculate the preferences for the center vs edge for bouts and interbouts
def calculate_socialfrac(bout_start, bout_end, interbout_start, interbout_end, xarr, fish_rois):
    # THIS IS HARDCODED BASED ON OUR SOCIAL RIG AND WOULD NEED TO BE MODIFIED IF THE SHAPE OR SETUP CHANGED
    # BASICALLY THE FISH ARE ON THE OUTSIDE, SO WE HAVE TO CONVERT THE XARR TO HAVE LARGE VALUES NEAR THE SOCIAL STIM ON ONE SIDE AND NOT THE OTHER
    # print(xarr)
    # could make this the halfway point of the image or something, but just going with a value I'm sure will work is ok until the code needs to be more flexible.
    if fish_rois[0] < 400:
        # the entire arrangement would need to change anyway if there were a change to the shape
        maxrealx = np.amax(xarr)
        zeroedgearr = xarr - maxrealx
        zeroedgearr = -1 * zeroedgearr
    else:
        minrealx = np.amin(xarr)
        zeroedgearr = xarr + abs(minrealx)

    maxrealx = np.amax(zeroedgearr)
    # deciding that zero is the farthest from stimulus fish
    cqpt = maxrealx * 0.25
    cqpttight = maxrealx * 0.10
    qpt = maxrealx * 0.75
    qpttight = maxrealx * 0.90
    csmoments = 0  # moments near control
    smoments = 0
    csmomentstight = 0
    smomentstight = 0
    intersmoments = 0
    tmoments = 0
    intertmoments = 0
    totalx = 0.0
    intertotalx = 0.0
# print("Qpt: ",qpt)
    for x3 in zeroedgearr[bout_start:(bout_end+1)]:
        totalx = x3 + totalx
        # print(totalx)
        if x3 < cqpt:
            csmoments = csmoments + 1
        if x3 > qpt:
            smoments = smoments + 1
        if x3 < cqpttight:
            csmomentstight = csmomentstight + 1
        if x3 > qpttight:
            smomentstight = smomentstight + 1
        tmoments = tmoments + 1

    csocialfrac = float(csmoments) / float(tmoments)
    socialfrac = float(smoments) / float(tmoments)
    csocialfraccloser = float(csmomentstight) / float(tmoments)
    socialfraccloser = float(smomentstight) / float(tmoments)

    avesocialfrac = (float(totalx) / float(tmoments)) / float(maxrealx)
    for x2 in zeroedgearr[interbout_start:(interbout_end+1)]:
        intertotalx = x2 + intertotalx
        if x2 > qpt:
            intersmoments = intersmoments + 1
        intertmoments = intertmoments + 1
# # Only comes up if there is a bout at end of movie, because then the interbout_end and _start are the same and don't really exist
# # None of this data will end up being analyzed though, because we ignore the last seconds usually anyway
    if intertmoments != 0.0:
        intersocialfrac = float(intersmoments) / float(intertmoments)
        interavesocialfrac = (float(intertotalx) /
                              float(intertmoments)) / float(maxrealx)
    else:
        intersocialfrac = 0.0
        interavesocialfrac = 0.0

    return socialfrac, avesocialfrac, intersocialfrac, interavesocialfrac, socialfraccloser, csocialfrac, csocialfraccloser