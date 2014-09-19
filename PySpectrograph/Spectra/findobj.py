#
# detect objects in a 2-D spectra
#
# This is adopted from apextract
#
import math
import numpy as np
import scipy.ndimage as nd

from PySpectrograph import PySpectrographError


def makeflat(data, method='median', specaxis=1):
    """Comparess a 2-D array along the spectral axis

       data

    """
    if method == 'median':
        return np.median(data, axis=specaxis)
    elif method == 'average':
        return np.average(data, axis=specaxis)
    elif method == 'sum':
        return np.sum(data, axis=specaxis)
    else:
        msg = '%s is not method to flatten array'
        raise PySpectrographError(msg)


def calc_ave_cont(data, method='sigclip', thresh=3.0, niter=5):
    """Calculate the average level of the continuum with sig-clip rejection or via median absolute deviation.
       If method is set to sigclip, it will return the sigma-clipped rejected mean and stardard deviation
       given using the threshold and number of iterations

       If method is mad, it will return the median and median absolute deviation for the data

       Returns mean, std
    """
    if method == 'mad':
        return np.median(data), np.median(abs(data - np.median(data)))

    ave = data.mean()
    std = data.std()
    for i in range(niter):
        mask = abs(data - ave) < std * thresh
        ave = data[mask].mean()
        std = data[mask].std()
    return ave, std


def calc_med_cont(data):
    """Calculate the median level of the continuum.  Std is based on MAD

       Returns median, std
    """
    med = np.median(data)
    std = 1.4826 * np.median(abs(data - med))
    return med, std


def findObjects(data, specaxis=1, method='median', thresh=3.0, niter=5, minsize=3):
    """Detect objects in the 2-D spectra
       data: 2D array of a spectral observations
       specaxis: Axis of the dispersion
       method: method for combining the array.  It can either be median, average, or sum
       thresh:  Threshhold for object detection
       niter:  Number of iterations

       return  list of tuples
    """

    # compress the data
    ldata = makeflat(data, method=method, specaxis=specaxis)

    # median the data
    ldata = nd.filters.median_filter(ldata, size=minsize)

    # determine the continuum values
    cont_mean, cont_std = calc_ave_cont(ldata, method='mad')

    return findLines(ldata, method=method, cont_mean=cont_mean, cont_std=cont_std,
                     thresh=thresh, niter=niter, minsize=minsize)


def findLines(ldata, method='median', cont_mean=None, cont_std=None, thresh=3.0, niter=5, minsize=3):
    """Detect objects in 1-D spectra
       ldatl: 1D array of a spectral observations
       specaxis: Axis of the dispersion
       method: method for combining the array.  It can either be median, average, or sum
       thresh:  Threshhold for object detection
       niter:  Number of iterations

       return  list of tuples
    """

    # set the levels of the continuum
    if cont_mean is None or cont_std is None:
        mean, std = calc_ave_cont(ldata)
    if cont_mean is None:
        cont_mean = mean
    if cont_std is None:
        cont_std = std

    # detect the peakc in the distribution
    obj_arr, obj_num = nd.label((abs(ldata - cont_mean) > thresh * cont_std))

    # determine the distributions
    obj_list = []
    mag_list = []

    # determine the boundries for all objects
    for i in range(1, obj_num + 1):
        ind = np.where(obj_arr == i)[0]
        my1 = ind.min()
        my2 = ind.max()
        if my2 - my1 > minsize and my1 < my2:
            objs = deblendObjects(ldata, my1, my2, thresh, niter, minsize)
            for y1, y2 in objs:
                if 0 < y1 < len(ldata) and 0 < y2 < len(ldata):
                    if y2 < y1:
                        y1, y2 = y2, y1
                    if y2 == y1:
                        y2 = y1 + 1
                    obj_list.append((y1, y2))
                    mag_list.append(ldata[y1:y2].max())

    # sort the objects in magnitude order
    ord_obj_list = []
    mag_arr = np.array(mag_list)
    mag_id = mag_arr.argsort()
    for i in mag_id[::-1]:
        ord_obj_list.append(obj_list[i])

    return ord_obj_list


def deblendObjects(ldata, y1, y2, thresh=3.0, niter=5, minsize=3):
    """Deblend a set of objects.  Deblend produces a list of y1,y2 for
           a an array created by scip.ndimages.label based on a set of data

    """

    # take the gradient of the data
    gdata = np.gradient(ldata[y1:y2])

    # determine if there is more than one object
    try:
        pos_ind = np.where(gdata >= 0)[0].max()
        neg_ind = np.where(gdata <= 0)[0].min()
    except:
        return [(y1, y2)]

    # If this is true, then there is only a single object to extract
    if abs(pos_ind - neg_ind) < minsize:
        return [(y1, y2)]

    # manually go through the points and determine where it starts and stops
    obj_list = []
    dy1 = y1
    neg = False
    for i in range(len(gdata)):
        if gdata[i] <= 0 and neg is False:
            neg = True
        if gdata[i] > 0 and neg is True:
            dy2 = dy1 + i
            obj_list.append((dy1, dy2))
            dy1 = dy2 + 1
            neg = False
    obj_list.append((dy1, y2))

    return obj_list


def plotdata(ldata, obj_arr):
    """Just for debuggin purposes"""
    from PySpectrograph.Utilities import makeplots

    nlen = len(ldata)
    x = np.arange(nlen)
    makeplots.figure(figsize=(6, 6), dpi=72)
    ay = makeplots.axes([0.15, 0.10, 0.8, 0.8])
    #ay.imshow(med_data, cmap=makeplots.cm.gray, aspect='equal', vmin=-5, vmax=50  )
    makeplots.plotline(ay, x, np.gradient(ldata))
    #makeplots.plotline(ay, x, gdata)
    makeplots.plotline(ay, x, obj_arr * 100)
    makeplots.show()
