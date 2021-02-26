import re
import numpy as np
import matplotlib.pyplot as plt
import sys
import six
from ipfnpytools.getsig import getsig
import logging
from matplotlib.backend_bases import NavigationToolbar2
import pandas as pd
import dd


_probe_jsat = pd.read_excel("Isat_signals.xlsx", header=None)
_probe_jsat = _probe_jsat.to_numpy()
probe_jsat_channel = {}
for row in _probe_jsat:
    if row[1] == 'Isat':
        probe_jsat_channel[row[0][-3:] + "-ch"] = row[2]
        probe_jsat_channel[row[0][-3:] + "-id"] = row[3] - 1

    
_autoscale_axes = []

# Add functionality to the HOME button
# After resetting to the default view, autoscale a certain group of axes
# that might have been drawn after the first render
home = NavigationToolbar2.home
def new_home(self, *args, **kwargs):
    home(self, *args, **kwargs)
    for ax in _autoscale_axes:
        ax.autoscale()
#         ax.set_xlim((1, 69))
NavigationToolbar2.home = new_home


LOG_FILENAME = 'log.txt'
logging.basicConfig(filename=LOG_FILENAME, level=logging.WARNING)


def hex_to_rgb(h):
    h = h.lstrip('#')
    return np.array(tuple(int(h[i:i+2], 16) for i in (0, 2, 4)))/255.0

def index(array, x):
    return np.where(np.array(array).flatten() == x)[0][0]

def _parse_line(line, rx_dict):
    """
    Do a regex search against all defined regexes and
    return the key and match result of the first matching regex
    """
    
    for key, rx in rx_dict.items():
        match = rx.search(line)
        if match:
            return key, match
    # if there are no matches
    return None, None


def _parse_spectroscopy_los(line):
    """
    Parse spectroscopy LOS from '/afs/ipp/u/sprd/loscoord/LOS_COORD_<YEAR>'
    """
    
    rx_dict = {
        'los': re.compile(r"'(?P<name>.*)' *(?P<R0>[^ ]*) *(?P<PHI0>[^ ]*) *(?P<Z0>[^ ]*) *(?P<R1>[^ ]*) *(?P<PHI1>[^ ]*) *(?P<Z1>[^ ]*) *\n"),
    }

    return _parse_line(line, rx_dict)


def _parse_diaggeom_langmuir(line):
    """
    Parse langmuir probe coordinates from diaggeom
    """
    
    rx_dict = {
        'probe': re.compile(r'  (?P<name>.*) \(Position\)\n'),
        'coordinates': re.compile(r'R=(?P<radius>.*)m, z=(?P<height>.*)m\n'),
        'metadata': re.compile(r'LSD \(Stat. Langmuir probes, (?P<divertor>.*) divertor\) #(?P<shot>.*) (?P<time>.*)s')
    }

    return _parse_line(line, rx_dict)


def get_spectroscopy_los_coordinates(filepath):
    """
    Get spectroscopy los coordinates from a file like '/afs/ipp/u/sprd/loscoord/LOS_COORD_<YEAR>'
    
    Parameters
    ----------
    filepath: str
        Path to the diaggeom-produced `.coords` file with the probes' coordinates.
    
    Returns
    -------
    meta_dict: dict
        Dictionary with the meta-data retrieved from the file. Usefull for debugging.
    names: str list
        Names of the probes as provided by `diaggeom`.
    radius: float list
        Radial position of each probe in meters R(m).
    height: float list
        Height of each probe in meters Z(m).
    """

    probes = {}

#     meta_dict = {'shot':0, 'divertor':'', 'time':0.0}

    # open the file and read through it line by line
    with open(filepath, 'r') as file_object:
        for line in file_object:
            key, match = _parse_spectroscopy_los(line)

#             print(key, match)

            if key == 'los':
                name = match.group('name').strip()
                probes[name] = {}
                for coord in ['R', 'PHI', 'Z']:
                    probes[name][coord] = [float(match.group(coord+'0')), float(match.group(coord+'1'))]

    return probes




def get_probe_coordinates(filepath):
    """
    Get divertor Langmuir probes' coordinates from a `.coords` file provided by `diaggeom`.
    
    Parameters
    ----------
    filepath: str
        Path to the diaggeom-produced `.coords` file with the probes' coordinates.
    
    Returns
    -------
    meta_dict: dict
        Dictionary with the meta-data retrieved from the file. Usefull for debugging.
    names: str list
        Names of the probes as provided by `diaggeom`.
    radius: float list
        Radial position of each probe in meters R(m).
    height: float list
        Height of each probe in meters Z(m).
    """

    names = []
    radius = []
    height = []

    meta_dict = {'shot':0, 'divertor':'', 'time':0.0}

    # open the file and read through it line by line
    with open(filepath, 'r') as file_object:
        for line in file_object:
            key, match = _parse_diaggeom_langmuir(line)

#             print(key, match)

            if key == 'probe':
                name = match.group('name')
                names.append(name)
#                 print("'"+name+"'")
            elif key == 'coordinates':
                r = float(match.group('radius'))
                radius.append(r)
                z = float(match.group('height'))
                height.append(z)
#                 print(r, z)
            elif key == 'metadata':
                meta_dict['divertor'] = match.group('divertor')
                meta_dict['shot'] = int(match.group('shot'))
                meta_dict['time'] = float(match.group('time'))
#                 print(meta_dict)

    return meta_dict, names, radius, height


def running_mean(x, pts, mode='same'):
    try:
        return np.convolve(x, np.ones((pts,))/pts, mode=mode)
    except:
        logging.exception("Error calculating the running mean of" + str(x))
        six.reraise(*sys.exc_info())
        

def post_process(signal, mean_pts):
    
    dt = signal.time[1] - signal.time[0]
    
    # Remove NaN's
    print("Removing NaN's")
    signal.time = signal.time[~np.isnan(signal.data)]
    signal.data = signal.data[~np.isnan(signal.data)]
    
    # Running mean over `mean_pts` points
    print("Computing a running mean over %d points <=> %.3lf ms" % (mean_pts, dt*mean_pts))
    signal.data = running_mean(signal.data, mean_pts)
    
    return signal


def get_jsat(shot, probe_name):
    """Get jsat data from Langmuir probes
    
    Parameters
    ----------
    shot: int
        Shot number, dah.
    probe_name: str
        Three-characters probe identifier according to diaggeom.
        
    Returns
    -------
    signal: ddGroupSignal
        signal.time is the time axis, signal.data is the jsat.
    """
        
    gettem = getsig(shot, 'LSF', 'CH' + str(probe_jsat_channel[probe_name + '-ch']))
    gettem.data = gettem.data[probe_jsat_channel[probe_name + '-id']]
    
    return gettem


def plot_probes(x, y, names, shot, axes_vessel=None, axes_jsat=None, axes_ne=None, axes_te=None, mean_pts=1):
    
    probe_traces = {}
    
    axes = {"jsat-": axes_jsat, "ne-": axes_ne, "te-": axes_te}
    
    for name in names:
        probe_traces["an-" + name] = name
    
    if axes_vessel is None:
        ax = plt.gca()
    else:
        ax = axes_vessel
    
    fig = ax.get_figure()
        
    rgba_colors = np.zeros((len(x), 4))
    rgba_colors[:, 3] = 0.5
    sc = ax.scatter(x, y, c=rgba_colors, zorder=500)
    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"),
                        zorder=1000)
    annot.set_visible(False)
    
    
    def toggle_visibility(artist):
        artist.set_visible(not artist.get_visible())
        
        
    def set_probe_color(probe_name):
        i = index(names, probe_name)
        rgba_colors[i, 0:3] = hex_to_rgb(probe_traces["color-" + probe_name])
        rgba_colors[i, 3] = 1
        sc.set_color(rgba_colors)
        
        
    def toggle_probe_color(probe_name):
        i = index(names, probe_name)
        if rgba_colors[i, 3] == 1:
            rgba_colors[i] = np.array((0, 0, 0, 0.5))
        else:
            rgba_colors[i, 0:3] = hex_to_rgb(probe_traces["color-" + probe_name]) 
            rgba_colors[i, 3] = 1
        sc.set_color(rgba_colors)
        
    def get_langmuir_signal(shot, probe_name, signal_type):
        
        if signal_type == 'jsat':
            return get_jsat(shot, probe_name)
        elif signal_type == 'ne' or signal_type == 'te':
            return getsig(shot, 'LSD', signal_type + '-' + probe_name)
        else:
            raise TypeError("signal_type must belong to {'ne', 'te', 'jsat'}")
        
    def show_probe(probe_name):
        
        i = index(names, probe_name)
        for prefix in ["te-", "ne-", "jsat-"]:
            try:
                probe_traces[prefix + probe_name].set_visible(True)
                rgba_colors[i, 0:3] = hex_to_rgb(probe_traces["color-" + probe_name]) 
                rgba_colors[i, 3] = 1
                sc.set_color(rgba_colors)
            except AttributeError, KeyError:
                pass
            
        probe_traces['plotted-' + probe_name] = True
            
    def hide_probe(probe_name):
        
        i = index(names, probe_name)
        for prefix in ["te-", "ne-", "jsat-"]:
            try:
                probe_traces[prefix + probe_name].set_visible(False)
                rgba_colors[i] = np.array((0, 0, 0, 0.5))
                sc.set_color(rgba_colors)
            except AttributeError, KeyError:
                pass
            
        probe_traces['plotted-' + probe_name] = False

            
    def plot_probe_trace(probe_name):
        
#         probe_already_plotted = False
        
#         for prefix in ["te-", "ne-", "jsat-"]:
#             try:
#                 toggle_visibility(probe_traces[prefix + probe_name])
#             except AttributeError:
#                 probe_already_plotted = True
#             except KeyError:
#                 pass
        
#         if probe_already_plotted:
#             toggle_probe_color(probe_name)
        
#         else:
        probe_traces['plotted-' + probe_name] = True

        for prefix, _ax in zip(["te", "ne", "jsat"], [axes_te, axes_ne, axes_jsat]):
            try:
                x = get_langmuir_signal(shot, probe_name, signal_type=prefix)
                x = post_process(x, mean_pts)

                try:
                    color = probe_traces["color-" + probe_name]
                except KeyError:
                    color = next(axes_te._get_lines.prop_cycler)['color']
                    probe_traces["color-" + probe_name] = color

                probe_traces[prefix + "-" + probe_name] = _ax.plot(x.time, x.data, color=color)[0]
                set_probe_color(probe_name)
            except:
                logging.exception("Shotfile not available")
                probe_traces[prefix + "-" + probe_name] = None
                probe_traces['an-' + probe_name] += ", " + prefix + " NA"


    def update_annot(ind):

        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        
        # This option shows all scatter points (good for overlaps)
#         text = "{}, {}".format(" ".join(list(map(str,ind["ind"]))), 
#                                " ".join([names[n] for n in ind["ind"]]))
        
        # This shows just one of the scatter points
#         text = "%s" % names[ind["ind"][0]] + ("" if success else " Not available")
        text = probe_traces['an-' + names[ind["ind"][0]]]
        
        annot.set_text(text)
    #     annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
    #     annot.get_bbox_patch().set_alpha(0.4)
    
    
    def click(event):
        try:
            vis = annot.get_visible()
            if event.inaxes == ax:
                cont, ind = sc.contains(event)
                if cont:
                    try: 
                        if probe_traces["plotted-" + names[ind["ind"][0]]]:
                            hide_probe(names[ind["ind"][0]])
                        else:
                            show_probe(names[ind["ind"][0]])
                    except KeyError:
                        plot_probe_trace(names[ind["ind"][0]])
                    except ValueError:
                        logging.exception("ValueError raise during plot_probe_trace")
                    update_annot(ind)
#                     plot_probe_trace(ind)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                    if vis:
                        annot.set_visible(False)
        except:
            logging.exception("Uncaught exception at the end of click handler")

#     fig.canvas.mpl_connect("motion_notify_event", hover)
    fig.canvas.mpl_connect("button_press_event", click)
    
    # Add axes for autoscalling once HOME is called
    _autoscale_axes.extend([axes_ne, axes_te])

    plt.show()
    
    