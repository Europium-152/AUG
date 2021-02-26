from ipfnpytools.spec_channels import spec_channels
from ipfnpytools.getsig import getsig
import dd
import re
import datetime as dt


_cached_ne = {'EVL': {}, 'GVL': {}, 'FVL': {}, 'HVL': {}}
_cached_channels = {'EVL': {}, 'GVL': {}, 'FVL': {}, 'HVL': {}}


def get_sbd(shot, channel):
    """
    Get the time dependent density measurements from  spectroscopy.
    
    Parameters
    ----------
    shot: int
        Shot number.
    channel: str
        Spectroscopy channel identifier, e.g., 'ROV-07' 
        
    Returns
    -------
    time: ndarray
        Array of time instants.
    ne: ndarray
        Array of densities.
    """
    spectroscopy_diags = ['EVL', 'GVL', 'FVL', 'HVL']
    
    for d in spectroscopy_diags:
        try:
            try: 
                i = _cached_channels[d][str(shot)].index(channel)
            except KeyError:
                _cached_channels[d][str(shot)] = spec_channels(shot, shotfile=d)
                i = _cached_channels[d][str(shot)].index(channel)
            index = i
            diag = d
            all_good = True
            break
        except (ValueError, dd.PyddError) as e:
            all_good = False
            pass
    if not all_good:
        raise ValueError("%s not available for shot %d" % (channel, shot))
           
    # Try chached ne for requested shot
    try:  
        ne_data = _cached_ne[diag][str(shot)]
    except KeyError:
        try:
            ne_data = getsig(shot, diag, 'Ne')
            _cached_ne[diag][str(shot)] = ne_data
        except dd.PyddError:
            print("NOT FOUND: shotfile(%d, '%s', 'Ne')" % (shot, diag))

    return (ne_data.time, ne_data.data[:,index])


def get_spec_channels_and_los(shot, shotfiles):
    
    
    spectroscopy_shotfiles = []
    for sf in shotfiles:
        spectroscopy_shotfiles.append(sf[:2] + 'S')
        

    for sf in spectroscopy_shotfiles:
        try:
            shot_date = dd.shotfile(sf, shot).date
            print("%s shotfile OK" % sf)
            shot_year = dt.datetime.strptime(
                shot_date if len(shot_date) == 18 else '0' + shot_date, 
                "%d%b%Y;%H:%M:%S"
            ).year
        except dd.PyddError:
            print("%s shotfile NOT FOUND" % sf)

    while True:
        try:
            filepath = "/afs/ipp/u/sprd/loscoord/LOS_COORD_%d" % shot_year
            los = get_spectroscopy_los_coordinates(filepath)
            break
        except IOError:
            shot_year -= 1
        if shot_year < 2005:
            raise TypeError("Did not find LOS coordinates for shot %d. Shot preceeds spectrocopy" % shot)

    print("Shot %d took place on %s" % (shot, shot_date))
    print("Using the LOS stored in file %s" % filepath)

    channels = {}

    print("Shot: %d" % shot)
    for sf in shotfiles:
        try:
            channels[sf] = spec_channels(shot, shotfile=sf)
            print(sf + " correctly oppened with %d channels" % len(channels[sf]))
        except dd.PyddError: 
            print(sf + " shotfile does not exist")

    # Cross channels with coordinates
    for sf in channels:
        for ch in channels[sf]:
            try:
                los[ch]
            except KeyError:
                channels[sf].remove(ch)
                print("%s-%s NOT FOUND inside coordinates file" % (sf, ch))
                
    return channels, los


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


def factors(n): 
    facts = list(set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0))))
    facts.sort()
    pairs = []
    for f1, f2 in zip(facts, facts[::-1]):
        pairs.append([f1, f2])
        
    return pairs[:len(pairs)/2 + (0!=len(pairs)%2)]


def best_column_number(n, form_factor=1.4):
    for i in range(n, 2*n):
        facts = factors(i)[-1]
        if (float(facts[1]) / float(facts[0])) > form_factor:
#             print('foda-se')
            pass
        else:
#             print('caralho')
            return facts[0], facts[1]