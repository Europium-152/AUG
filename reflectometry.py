import ipfnpytools.aug_read as aug_read
from ipfnpytools import rps_dump
import matplotlib.pyplot as plt
import numpy as np

def fetch(shot, dump=None, verbose=False):
    """Fetch density data from the RPS shotfile of the O-mode reflectometer
    
    Parameters
    ----------
    shot: int
        Shot number. `shot=0` corresponds to the latest shot available.
    dump: str, optional
        Optionaly you can provide a relative or absolute path to an RPS dump file. In this case, `shot`  is ignored
    verbose: bool, optional
        Show debugging messages. Defaults to `False`
        
    Returns
    -------
    time: ndarray
        Time instants. Shape (T, ).
    lfs_density: ndarray
        2D array with probbing density values for the LFS data. Shape (T, N).
    hfs_density: ndarray
        2D array with probbing density values for the HFS data. Shape (T, N).
    lfs_radius: ndarray
        2D ndarray with radius of plasma surface as a function of time and density. Shape (T, N).
    hfs_radius: ndarray
        2D ndarray with radius of plasma surface as a function of time and density. Shape (T, N).
    """
    
    # Fetching data ---------------------------------------------------

    if dump is None:

        data = aug_read.many_signals(
            diagnostics=(["RPS"] * 2),
            names=['neb_LFS', 'neb_HFS'],
            shots=shot)

        # Time array is equal for LFS and HFS
        time = data.times[0]

        # Low-field side
        lfs_signal = data.signals[0]
        lfs_area = data.areas[0]

        # High-field side
        hfs_signal = data.signals[1]
        hfs_area = data.areas[1]

    else:

        # Instantiate the ShotFile class with the path to the dump file
        shotfile = rps_dump.ShotFile(path)

        ne_lfs = shotfile("neb_LFS")
        ne_hfs = shotfile("neb_HFS")

        # Time array is equal for LFS and HFS
        time = ne_lfs.time

        # Low-field side
        lfs_signal = ne_lfs.data
        lfs_area = ne_lfs.area.data

        # High-field side
        hfs_signal = ne_hfs.data
        hfs_area = ne_hfs.area.data


    if verbose:
        print("time:", time.shape)
        print("lfs_signal:", lfs_signal.shape)
        print("lfs_area:", lfs_area.shape)
        print("hfs_signal:", hfs_signal.shape)
        print("hfs_area:", hfs_area.shape)
        
    return time, lfs_signal, hfs_signal, lfs_area, hfs_area



def update_errorbar(errobj, x, y, xerr=None, yerr=None, lowxerr=None, upxerr=None, lowyerr=None, upyerr=None):
    """Update an errorbar object. Adapted from https://github.com/matplotlib/matplotlib/issues/4556
    
    Parameters
    ----------
    x: ndarray
        x data. Used in Line2D.set_xdata().
    y: ndarray
        y data. Used in Line2D.set_ydata().
    x(y)err: float or ndarray, optional
        Error in x(y). Scalar or array with the same dimension as `x(y)`.
    lowx(y)err: float or ndarray, optional
        Lower bound error in x(y). Scalar or array with the same shape as `x(y)`
    upx(y)err: float or ndarray, optional
        Upper bound error in x(y). Scalar or array with the same shape as `x(y)`  
    """
    ln, caps, bars = errobj

    if len(bars) == 2:
        assert (xerr is not None and yerr is not None) or \
               (lowxerr is not None and upxerr is not None and lowyerr is not None and upyerr is not None),\
               "Your errorbar object has 2 dimension of error bars defined. You must provide xerr and yerr."
        barsx, barsy = bars  # bars always exist (?)
        try:  # caps are optional
            errx_top, errx_bot, erry_top, erry_bot = caps
        except ValueError:  # in case there is no caps
            pass

    elif len(bars) == 1:
        assert (xerr is     None and yerr is not None) or\
               (xerr is not None and yerr is     None) or\
               (lowxerr is None and upxerr is None and lowyerr is not None and upyerr is not None) or\
               (lowxerr is not None and upxerr is not None and lowyerr is None and upyerr is None),  \
               "Your errorbar object has 1 dimension of error bars defined. You must provide xerr or yerr."

        if (xerr is not None) or (lowxerr is not None and upxerr is not None):
            barsx, = bars  # bars always exist (?)
            try:
                errx_top, errx_bot = caps
            except ValueError:  # in case there is no caps
                pass
        else:
            barsy, = bars  # bars always exist (?)
            try:
                erry_top, erry_bot = caps
            except ValueError:  # in case there is no caps
                pass

    ln.set_data(x,y)
    
    if lowxerr is None: lowxerr = xerr
    if upxerr is None: upxerr = xerr
    if lowyerr is None: lowyerr = yerr
    if upyerr is None: upyerr = yerr

    try:
        errx_top.set_xdata(x + upxerr)
        errx_bot.set_xdata(x - lowxerr)
        errx_top.set_ydata(y)
        errx_bot.set_ydata(y)
    except NameError:
        pass
    try:
        barsx.set_segments([np.array([[xt, y], [xb, y]]) for xt, xb, y in zip(x + upxerr, x - lowxerr, y)])
    except NameError:
        pass

    try:
        erry_top.set_xdata(x)
        erry_bot.set_xdata(x)
        erry_top.set_ydata(y + upyerr)
        erry_bot.set_ydata(y - lowyerr)
    except NameError:
        pass
    try:
        barsy.set_segments([np.array([[x, yt], [x, yb]]) for x, yt, yb in zip(x, y + upyerr, y - lowyerr)])
    except NameError:
        pass