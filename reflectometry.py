import types
import ipfnpytools.aug_read as aug_read
from ipfnpytools import rps_dump
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import ipywidgets as widgets
from ipfnpytools.closest import closest

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




def fill_between_x(x1, x2, y, color, alpha):
    
    def _set_data(self, x1, x2, y):
        self.set_xy([[_x, _y] for _x, _y in zip(x1, y)] + [[_x, _y] for _x, _y in reversed(zip(x2, y))])
        
    p = Polygon(xy=[[0, 0]]*3, color=color, alpha=alpha, linewidth=0)
    p.set_data = types.MethodType(_set_data, p)
    p.set_data(x1, x2, y)  
    
    return p

def two_plus_one(x1, y1, x23, y2, y3, z2, z3, ey2=None, ey3=None, number_of_profiles=1, lx1='', ly1='', ly2='', ly3='', lz2='', lz3='', m2=None, m3=None, strip_pts=1, show_same=False, paint_rho=False, paint_wall=False, sharex=False, sharey=False):
    
    
    def lb(x, i):
        return x[max(0, i)]
    def ub(x, i):
        return x[min(len(x), i)]
    
    show_wall = (m2 is not None) and (m3 is not None)
    show_std = (ey2 is not None) and (ey3 is not None)
    
    # Determine which time indexes draw profiles
    draw = np.linspace(0, len(x23), number_of_profiles + 2, dtype=np.int)
    draw = draw[1:-1]
    
    N = int(strip_pts/2)
    
    fig = plt.figure()
    plt.subplots_adjust(hspace=0.6)

    # Supporting plot with colored lines
    ax1 = plt.subplot(212)
    stamp = []
    ax1.plot(x1, y1, color='k')
    for i in draw:
        color = next(ax1._get_lines.prop_cycler)['color']
        if N == 0:
            l = ax1.axvline(x23[i], color=color)
        else:
            l = ax1.axvspan(lb(x23, i-N), ub(x23, i+N), color=color, alpha=0.6)
            
        stamp.append(l)

    if show_same:
        ax2 = plt.subplot(211)
        ax3 = ax2
    else:
        ax2 = plt.subplot(221)
        ax3 = plt.subplot(222, sharey=ax2 if sharey else None, sharex=ax2 if sharex else None)
        
    if paint_rho: 
        ax2.axvspan(0.0, 1.0, color='#FFC0CB')
        if not show_same:
            ax3.axvspan(0.0, 1.0, color='#FFC0CB')
                
    if paint_wall:
        ax2.axvspan(0.0, 1.045, color='#D3D3D3')
        if not show_same:
            ax3.axvspan(2.22, 5, color='#D3D3D3')

        
    l2 = []  # Lines
    p2 = []  # Polygons
    w2 = []  # Walls
    
    l3 = []
    p3 = []
    w3 = []

    legend_lines = []
    
    for i in draw:
        color = next(ax2._get_lines.prop_cycler)['color']

        legend_lines.append(Line2D([0], [0], color=color, linewidth=1.5, linestyle='-'))

        l, = ax2.plot(y2[i], z2[i], color=color)
        l2.append(l)

        l, = ax3.plot(y3[i], z3[i], color=color, linestyle='--' if show_same else '-')
        l3.append(l)

        if show_std:
            p = fill_between_x(y2[i] - ey2[i], y2[i] + ey2[i], z2[i], color=color, alpha=0.5)
            ax2.add_patch(p)
            p2.append(p)    

            p = fill_between_x(y3[i] - ey3[i], y3[i] + ey3[i], z3[i], color=color, alpha=0.5)
            ax3.add_patch(p)
            p3.append(p) 

        if show_wall:
            l = ax2.axvline(m2[i], color=color)
            w2.append(l)
            l = ax3.axvline(m3[i], color=color, linestyle='--' if show_same else '-')
            w3.append(l)

    if show_same:
        solid_line = Line2D([0], [0], color='k', linewidth=1.5, linestyle='-')
        dashed_line = Line2D([0], [0], color='k', linewidth=1.5, linestyle='--')
        legend = ax3.legend(handles=legend_lines + [solid_line, dashed_line],
                            labels = ["%.3lf s" % x23[i] for i in draw] + ['HFS', 'LFS'],
                            loc='center left', bbox_to_anchor=(1, 0, 1, 1))
    else:
        legend = ax3.legend(handles=legend_lines,
                            labels = ["%.3lf s" % x23[i] for i in draw],
                            loc='center left', bbox_to_anchor=(1, 0, 1, 1))

    ax1.set_xlabel(lx1)
    ax2.set_xlabel(ly2)
    if not show_same: ax3.set_xlabel(ly3)

    ax1.set_ylabel(ly1)
    ax2.set_ylabel(lz2)
              

    def update(**kwargs):
        slider_x = kwargs.values()
        slider_x.sort()
        for k, ts in enumerate(slider_x):
            i = closest(x23, ts)
            l2[k].set_xdata(y2[i])
            l3[k].set_xdata(y3[i])
            if show_std:
                p2[k].set_data(y2[i] - ey2[i], y2[i] + ey2[i], z2[i])
                p3[k].set_data(y3[i] - ey3[i], y3[i] + ey3[i], z3[i])
               
            if N==0:
                stamp[k].set_xdata(2*[x23[i]])
            else:
                stamp[k].set_xy(
                [[lb(x23, i-N), 0.        ],
                 [lb(x23, i-N), 1.        ],
                 [ub(x23, i+N), 1.        ],
                 [ub(x23, i+N), 0.        ],
                 [lb(x23, i-N), 0.        ]]
                )

            if show_wall:
                w2[k].set_xdata(2*[m2[i]])
                w3[k].set_xdata(2*[m3[i]])
            legend.get_texts()[k].set_text('%.3lf s' % x23[i])


    sliders = []
    for i in range(number_of_profiles):
        sliders.append(widgets.FloatSlider(
            value=x23[draw[i]],
            min=x1[0],
            max=x1[-1],
            step=x1[1]-x1[0],
            description=lx1,
            disabled=False,
            continuous_update=True,
            orientation='horizontal',
            readout=True,
            readout_format='.6f',
        ))

    kwargs = {'p{0}'.format(i):slider for i, slider in enumerate(sliders)}

    widgets.interact(update, **kwargs)
    
    return fig, ax1, ax2, ax3
    

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