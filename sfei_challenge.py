import csv
import pandas as pd
import time, datetime, dateutil
from math import pi, sin, cos, asin, sqrt
import operator
import matplotlib.pyplot as plt
import matplotlib.patches as mp
import numpy as np

LATS = [37.419724, 38.154395]
LONS = [-121.921501, -122.546033]
BAY_RECT = [(37.419724, -121.921501), (38.154395, -122.546033)]
SOUTH_BAY_RECT = [(37, 441954, -121.92515), (37.489369, -122.13698)]

R_EARTH = 3963.0  # in miles
DEG_PER_MILE_LAT = 1.0 / 69.0
DEG_PER_MILE_LON = 1.0 / 54.6

SQUARE_SIZE = 2.5
DLAT = DEG_PER_MILE_LAT * SQUARE_SIZE
DLON = DEG_PER_MILE_LON * SQUARE_SIZE
N_LATS = int((LATS[1] - LATS[0]) / DLAT + 2)
N_LONS = int((LONS[0] - LONS[1]) / DLON + 2)

GRID = list()
for i in range(N_LATS):
    for j in range(N_LONS):
        lat = LATS[0] + i * DLAT
        lon = LONS[0] - j * DLON
        GRID.append((lat, lon))

OLDEST_DATE = dateutil.parser.parse('01/01/1950')
DATE0 = dateutil.parser.parse('01/01/2000')
DATE1 = dateutil.parser.parse('01/01/1985')

UNIT_DICT = {'ng/g': 1000.0, 'ug/g': 1.0, 'mg/kg': 1.0, 'ug/kg': 1000.0}
WET_DICT = {'ww': 1.0, 'dw': 5.0}
# UNIT_DICT = {'ng/g dw':1000.0, 'ug/g dw':5.0, 'ug/g ww':1.0, 'mg/kg dw':1.0}
PIER_DICT = {'WarmCove':(37.754612, -122.382454), 'Oakland':(37.784290, -122.246733), 'SanFran':(37.810345, -122.424451),
             'Berkley':(37.862886, -122.318493), 'Valejo':(38.094715, -122.259788), 'PointPinole':(38.013905, -122.368568),
             'Pacifica':(37.600728, -122.504744), 'CoyoteCreek':(37.444018, -121.958745),
             'CandleStick':(37.708926, -122.373535), 'CoyotePoint':(37.592322, -122.313600),
             'Benicia':(38.048339, -122.165400), 'SturgeonHole':(38.035922, -122.432689)}

RADIUS = 2.5  # miles
THRESH_TIS = 15   # Need at least so many tissue samples per station
THRESH_SED = 3    # Need at least so many sediment samples per station
THRESH = 10
DANGER = 1.0

def translate_unit(unit):
    try:
        u = unit.split()

    except AttributeError:
        return None
    if u[0] not in UNIT_DICT or u[1] not in WET_DICT:
        return None
    else:
        div0 = UNIT_DICT[u[0]]
        div1 = WET_DICT[u[1]]
    return tuple((div0, div1))

def add_station_code(fname):
    tmp = fname.replace('.csv', '')
    out_fname = tmp.replace("_master", "")+'_staID.csv'
    with open(fname, 'r') as f, open(out_fname, 'w') as of:
        reader = csv.DictReader(f)
        l1 = list()
        for row in reader:
            l1.append(row['stationname'])
        l1 = list(set(l1))
        d1 = {l1[i] : i for i in range(len(l1))}

        fieldnames = ['StationID']+reader.fieldnames
        writer = csv.DictWriter(of, fieldnames)
        writer.writeheader()
        f.seek(0)
        reader = csv.DictReader(f)
        for row in reader:
            row['StationID'] = d1[row['stationname']]
            writer.writerow(dict(row))
    
def box_coords(fname, rect = BAY_RECT):
    tmp = fname.replace('.csv', '')
    in_fname = tmp.replace("_master", "")+'_staID.csv'
    out_fname = "bay_"+in_fname
    latlo = rect[0][0]
    lathi = rect[1][0]
    lones = rect[0][1]
    lonws = rect[1][1]
    first = True
    df = pd.read_csv(in_fname)
    for idx in range(len(df)):
        lat = df.iloc[idx]['latitude']
        if lat < latlo or lat > lathi:
            continue
        lon = df.iloc[idx]['longitude']
        if lon > lones or lon < lonws:
            continue
        mode = 'w' if first else 'a'
        header = True if first else False
        df.iloc[idx:idx+1].to_csv(out_fname, index=False, mode=mode, header=header)
        first = False
    return out_fname

''' Create a dictionary with key=stationID, value = (lat, lon) tuple'''
def sta_dict_coord(fname):
    df = pd.read_csv(fname)
    df = df.sort_values('StationID')
    d1 = dict()
    for idx, row in df.iterrows():
        key = row['StationID']
        if key not in d1:
            d1[key] = (row['latitude'], row['longitude'])
    return d1

def get_date(sample_date):
    date = dateutil.parser.parse(sample_date)
    return date

''' Create a dictionary, key=StationID, value = list of results'''
def sta_dict_res(fname, start_date=DATE0, mussel=None):
    df = pd.read_csv(fname)
    df = df.sort_values('StationID')
    d1 = dict()
    for idx, row in df.iterrows():
        date = get_date(row['sampledate'])
        if not date:
            print("Error on idx = %d" % idx)
            print(date)
            return None
        if date < start_date:
            continue
        if np.isnan(row['result']) or row['result'] == 0:
            continue
        if mussel and 'mussel' not in row['commonname'].lower():
            continue
        # Exclude unknown units (they show up very rarely or belong to freshwater etc..)
        div = translate_unit(row['unitname'])
        if not div:
            continue
        # Translate units if needed
        row['result'] /= div[0]
        if 'tissue' in fname:
            row['result'] /= div[1]
        key = row['StationID']
        if key in d1:
            d1[key].append(row['result'])
        else:
            d1[key] = [row['result']]
    return d1

''' Simply calculate averages, discarding stas with too few data
    Inputs are dictionary with stationID, value_list and discard threshold'''
def sta_avg(sdv, thresh):
    d1 = dict()
    for key in sdv:
        if len(sdv[key]) >= thresh:
            d1[key] = sum(sdv[key])/len(sdv[key])
    return d1

''' Create a dictionary key=StationID, value = list of results for given year'''
def sta_year_dict(fname, year):
    df = pd.read_csv(fname)
    df = df.sort_values('StationID')
    d1 = dict()
    for idx, row in df.iterrows():
        date = get_date(row['sampledate'])
        if not date:
            return None
        if date.year != year:
            continue
        # Translate units if needed
        div = translate_unit(row['unitname'])
        if not div:
            continue
        row['result'] /= div[0]
        if 'tissue' in fname:
            row['result'] /= div[1]
        key = row['StationID']
        if key in d1:
            d1[key].append(row['result'])
        else:
            d1[key] = [row['result']]
    return d1

# Calculates distance in miles between two (lat, long) coordinates
def distance(coord1, coord2):
    def toRadians(deg):
        return deg * pi / 180

    lat1 = toRadians(coord1[0])
    lon1 = toRadians(coord1[1])
    lat2 = toRadians(coord2[0])
    lon2 = toRadians(coord2[1])
    dlat = lat1 - lat2
    dlon = lon1 - lon2
    tmp = sin(dlat/2.0)**2 + cos(lat1)*cos(lat2)*sin(dlon/2.0)**2
    return 2.0*asin(sqrt(tmp)) * R_EARTH

# center coords, other coords
def group_coords(ccord, ocord, radius):
    d1 = dict()
    for ckey, center in ccord.iteritems():
        for okey, cord in ocord.iteritems():
            if distance(center, cord) < radius:
                if ckey in d1:
                    d1[ckey].append(okey)
                else:
                    d1[ckey] = [okey]
    return d1

'''
Out: {year: [list of values]} for pier and fish
'''
def by_pier_and_fish(pier, inf, fish=None, start_date=OLDEST_DATE, radius = RADIUS):
    d_coord = sta_dict_coord(inf)
    stas = find_circle(PIER_DICT[pier], d_coord, radius=radius)
    out = dict()
    df = pd.read_csv(inf)
    for idx, row in df.iterrows():
        # Filter by fish
        if fish and fish not in row['commonname']:
            continue
        # Some entries are empty..
        if np.isnan(row['result']) or row['result'] == 0:
            print(row['result'])
            continue
        # Filter by date
        date = get_date(row['sampledate'])
        if date < start_date:
            continue
        if row['StationID'] not in stas:
            continue
        # Exclude unknown units (they show up very rarely or belong to freshwater etc..)
        div = translate_unit(row['unitname'])
        if not div:
            continue
        # Translate units if needed
        row['result'] /= div[0]
        if 'tissue' in inf:
            row['result'] /= div[1]
        if date.year not in out:
            out[date.year] = [row['result']]
        else:
            out[date.year].append(row['result'])
    return out

def process_row(row, use_div1):
    if np.isnan(row['result']):
        return None
    date = get_date(row['sampledate'])
    if date < DATE0:
        return None
    div = translate_unit(row['unitname'])
    if not div:
        return None
    # Translate units if needed
    row['result'] /= div[0]
    if use_div1:
        row['result'] /= div[1]
    return row['result']

def grid_corr(inf1, inf2):
    #{staID: (lat, lon)}
    coord1 = sta_dict_coord(inf1)
    coord2 = sta_dict_coord(inf2)
    # {staID: [list of results]}
    d1 = sta_dict_res(inf1, start_date=DATE0)
    d2 = sta_dict_res(inf2, start_date=DATE0)
    g = list()
    for grid in GRID:
        stas = find_circle(grid, coord1, SQUARE_SIZE/sqrt(2.0))
        val1 = avg_within_circle(stas, d1, 5)
        stas = find_circle(grid, coord2, SQUARE_SIZE/sqrt(2.0))
        val2 = avg_within_circle(stas, d2, 5)
        if val1 and val2:
            g.append((val1, val2))
    return g

''' dictionary {pier: [average_tissue, average_sediment]} '''
def answer_pier_contamination(infa, infb, fish=None, start_date=OLDEST_DATE, radius=RADIUS):
    '''
    add_station_code(infa)
    add_station_code(infb)
    infa = box_coords(infa)
    infb = box_coords(infb)
    '''
    d1 = dict()
    for pier in PIER_DICT:
        print('doing infa, pier=%s' % pier)
        d2 = by_pier_and_fish(pier, infa, fish, start_date, radius)
        avg_a = 0.0
        num = 0
        for year in d2:
            avg_a += sum(d2[year])
            num += len(d2[year])
        if num == 0:
            avg_a = 0.0
        else:
            avg_a = avg_a / (1.0 * num)
        print('doing infa, pier=%s' % pier)
        d2 = by_pier_and_fish(pier, infb, fish=None, start_date=start_date, radius=radius)
        avg_b = 0.0
        num = 0
        for year in d2:
            avg_b += sum(d2[year])
            num += len(d2[year])
        if num == 0:
            avg_b = 0.0
        else:
            avg_b = avg_b / (1.0 * num)
        d1[pier] = [avg_a, avg_b]
    return d1

def plot_barh(d1):
    l = len(d1)
    w = 1.0
    ind1 = np.array([i+i/2 for i in range(0, 2*int(w)*l, 2*int(w))])
    ind2 = np.array([i+i/2 for i in range(1*int(w), 2*int(w)*l, 2*int(w))])
    fig, ax = plt.subplots(figsize=(10.0, 10.0))
    ax.patch.set_facecolor([0.8, 0.5, 0.5])
    y = list()
    x = list()
    for key in d1:
        x.append(key)
        if d1[key][0] > 1.0:
            y.append(1.0)
        else:
            y.append(d1[key][0])
    ax.barh(ind1, y, w, color='b', label='Tissue')
    shim = max(y) / 100.0
    for i, val in enumerate(y):
        text = str(float("{:.4f}".format(val)))
        if val > 0.99:
            text = "> 1.0 (3.4)"
        ax.text(val+shim, ind1[i]+0.25, text, color='k', fontsize='small')
    y = list()
    x = list()
    for key in d1:
        x.append(key)
        if d1[key][1] > 1.0:
            y.append(1.0)
        else:
            y.append(d1[key][1])
    ax.barh(ind2, y, w, color='m', label='Sediment')
    shim = max(y) / 100.0
    for i, val in enumerate(y):
        text = str(float("{:.4f}".format(val)))
        if val > 0.99:
            text = "> 1.0 (1.6)"
        ax.text(val+shim, ind2[i]+0.25, text, color='k', fontsize='small')
    ytics = [l+w for i, l in enumerate(ind1)]
    ax.set_yticks(ytics)
    ax.set_yticklabels(x, minor=False)
    plt.xlabel('PCB, parts per million (ppm)')
    plt.legend(ncol=2, bbox_to_anchor=(.6,.1),
              loc='lower left', fontsize='small')
    plt.title('PCB Concentration in Sediment and Tissue at Popular Piers')
    plt.show()

# Input: center coord, dict of points {stationID : coord}, radius
# Output: list of stationIDs within radius from center
def find_circle(center, points, radius):
    l1 = list()
    for pkey, cord in points.iteritems():
        if distance(center, cord) < radius:
            l1.append(pkey)
    return l1


# Inputs: stationIDs within circle, results dictionary
def avg_within_circle(sta_ids, results, thresh=1):
    res = 0.0
    N = 0
    for id in sta_ids:
        if id in results:
            res += sum(results[id])
            N += len(results[id])
    if N > thresh:
        return res / (1.0 * N)
    else:
        return None

# Inputs: mean values dictionary, other values dictionary and coordinate dictionaries.
def calc_cross(mean_vals, other_vals, coord1, coord2):
    l1 = list()
    for id, mean in mean_vals.iteritems():
        center = coord1[id]
        sta_ids = find_circle(center, coord2, RADIUS)
        other_mean = avg_within_circle(sta_ids, other_vals)
        if other_mean:
            l1.append((mean, other_mean))
    return l1

''' rect = [(lat_lo, lon_e), (lat_hi, lon_w)]
    coord = (lat, lon)
'''
def in_rectangle(rect, coord, radius=None):
    if radius:
        if distance(rect, coord) < radius:
            return True
        return False
    else:
        if coord[0] < rect[0][0] or coord[0] > rect[1][0]:
            return False
        if coord[1] > rect[0][1] or coord[1] < rect[1][1]:
            return False
        return True


''' Within a rectangular area and for a given fish (None for all fish)
    collect values for each year
    {year: [list of values]}
    d1, d2 = fish_by_year_dict(inf, rect=BAY_RECT, fish='Bass')
'''
def fish_by_year_dict(fname, rect=SOUTH_BAY_RECT, fish=None, radius=None):
    df = pd.read_csv(fname)
    d1 = dict()
    for idx, row in df.iterrows():
        # Filter by fish
        if fish:
            for one_fish in fish:
                if one_fish in row['commonname'].lower():
                    break
            else:
                continue
        # Filter by area
        if not in_rectangle(rect, (row['latitude'], row['longitude']), radius=radius):
            continue
        if np.isnan(row['result']) or row['result'] == 0:
            continue
        # Filter and translate units
        div = translate_unit(row['unitname'])
        if not div:
            continue
        val = row['result'] / div[0]
        if 'tissue' in fname:
            val /= div[1]
        if val > 3.0:
            continue
        date = get_date(row['sampledate'])
        if not date:
            print("ERROR: No valid date found")
            return None
        year = date.year
        if year not in d1:
            d1[year] = [val]
        else:
            d1[year].append(val)
    d2 = dict()
    d3 = dict()
    for year in d1:
        vals = np.array(d1[year])
        mean = np.mean(vals)
        median = np.median(vals)
        d2[year] = mean
        d3[year] = median
    return d1, d2, d3

"{pier: {year: [list of vals]}}"
def answer_longitudinal(fname, piers, fish=None, radius=RADIUS):
    pier_d = dict()
    for pier in piers:
        _, d1, _ = fish_by_year_dict(fname, PIER_DICT[pier], fish=fish, radius=radius)
        pier_d[pier] = d1
    for pier in pier_d:
        l = [(year, pier_d[pier][year]) for year in sorted(pier_d[pier])]
        write_pair_vals(l, pier+'.txt')
    return pier_d


# Helper function to write pairs of values to a file
# pairs is either a list of tuples or a dictionary
def write_pair_vals(pairs, fout):
    if isinstance(pairs, dict):
        pairs = sorted(pairs.items())
    with open(fout, 'w') as wf:
        writer = csv.writer(wf, delimiter=',')
        for pair in pairs:
            row = [pair[0], pair[1]]
            writer.writerow(row)

'''Helper function to write dict values in two columns,
   first column is dict key and second is value
'''
def write_dict_vals(di, fout):
    with open(fout, 'w') as wf:
        writer = csv.writer(wf, delimiter=',')
        for key in di:
            for el in di[key]:
                row = [key, el]
                writer.writerow(row)

def square_means(res_tis, res_sed, coord_tis, coord_sed):
    d1 = dict()
    for grid in GRID:
        lat0 = grid[0] - DLAT / 2.0
        lat1 = grid[0] + DLAT / 2.0
        lon0 = grid[1] + DLON / 2.0
        lon1 = grid[1] - DLON / 2.0
        l1 = list()
        for key, coord in coord_tis.iteritems():
            if coord[0] < lat0 or coord[0] > lat1 \
                or coord[1] > lon0 or coord[1] < lon1:
                continue
            l1.append(key)
        l2 = list()
        for key, coord in coord_sed.iteritems():
            if coord[0] < lat0 or coord[0] > lat1 \
                or coord[1] > lon0 or coord[1] < lon1:
                continue
            l2.append(key)
        d1[grid] = (l1, l2)
    return d1


# Extract top perc % from a dictionary of {stationID : value}
def get_top_percent(di, perc):
    # Sort by value
    sorted_l = sorted(di.items(), key = operator.itemgetter(1))
    idx = len(sorted_l) * (100-perc) / 100
    top_di = dict(sorted_l[idx:])
    return top_di


def filter_top_by_area(di, coord, radius):
    do = dict()
    for key in di:
        for key_o in do:
            if distance(coord[key], coord[key_o]) < radius:
                break
        else:
            do[key] = di[key]
    return do

def place_top_percent(di, coord_d, outf):
    di = filter_top_by_area(di, coord_d, 1.2)
    with open(outf, 'w') as wf:
        writer = csv.writer(wf, delimiter = ',')
        for key in di:
            row = [coord_d[key][0], coord_d[key][1], di[key]]
            writer.writerow(row)


# answer_high_contamination('merc_sed_orig.csv', 'merc_sed_top.csv', 25, ('area', 1.0))
# answer_high_contamination('merc_tis_orig.csv', 'merc_tis_top.csv', 25, ('thresh', 10))
def answer_high_contamination(inf, outf, perc, avg_method):
    add_station_code(inf)
    bay_file = box_coords(inf)
    coord = sta_dict_coord(bay_file)
    res_date = sta_dict_res(bay_file, start_date=DATE0)
    avg_d = dict()
    if avg_method[0] == 'area':
        radius = avg_method[1]
        for key in res_date:
            points = find_circle(coord[key], coord, radius)
            assert key in points
            l1 = list()
            for point in points:
                l1 += res_date[point]
            avg_d[key] = sum(l1) / (1.0*len(l1))
    else:
        thresh = avg_method[1]
        for key in res_date:
            if len(res_date[key]) >= thresh:
                avg_d[key] = sum(res_date[key]) / (1.0 * len(res_date[key]))
    top_di = get_top_percent(avg_d, perc)
    place_top_percent(top_di, coord, outf)


# answer_correlation(infa, infb, 'corr_merc_meth.csv', radius=2.0, thresh=10)
def answer_correlation(infa, infb, outf, radius=RADIUS, thresh=THRESH):
    '''
    :param infa: tissue master file
    :param infb: sediment master file
    :param outf: output file
    :param radius:
    :param thresh:
    :return: write file [stationID(tissue), sediment, tissue]
    '''
    '''
    add_station_code(infa)
    add_station_code(infb)
    infa = box_coords(infa)
    infb = box_coords(infb)
    '''
    coord_a = sta_dict_coord(infa)
    coord_b = sta_dict_coord(infb)
    # {stationID: [list of values]}
    res_date_a = sta_dict_res(infa, start_date=DATE1)
    res_date_b = sta_dict_res(infb, start_date=DATE1)
    # {stationID(tissue): [list of stationIDs(sediment) within radius
    grcod = group_coords(coord_a, coord_b, radius)
    mean_vals_a = sta_avg(res_date_a, thresh=thresh)
    l1 = list()
    # For each stationID(tissue) find mean of sediment within radius
    for key in grcod:
        mean_val_b = avg_within_circle(grcod[key], res_date_b, thresh=thresh)
        if mean_val_b and key in mean_vals_a:
            l1.append([key, mean_val_b, mean_vals_a[key]])
    with open(outf, 'w') as wf:
        writer = csv.writer(wf, delimiter=',')
        for row in l1:
            writer.writerow(row)

'''
f_tis = 'merc_tis_orig.csv'
f_sed = 'merc_sed_orig.csv'
add_station_code(f_tis)
add_station_code(f_sed)
f_bay_tis = box_coords(f_tis)
f_bay_sed = box_coords(f_sed)

# Dictionaries: stationID : coordinate
coord_tis = sta_dict_coord(f_bay_tis)
coord_sed = sta_dict_coord(f_bay_sed)
circles = group_coords(coord_sed, coord_tis, RADIUS)

res_tis = sta_dict_res(f_bay_tis, start_date=DATE0)
res_sed = sta_dict_res(f_bay_sed, start_date=DATE0)


# Dictionaries: stationID : mean value
mean_tis = sta_avg(res_tis, thresh=THRESH_TIS)
mean_sed = sta_avg(res_sed, thresh=THRESH_SED)

l1 = calc_cross(mean_tis, res_sed, coord_tis, coord_sed)
with open('cross_circles.csv', 'w') as wf:
    writer = csv.writer(wf, delimiter=',')
    for pair in l1:
        row = [pair[0], pair[1]]
        writer.writerow(row)

sorted_tis = sorted(mean_tis.items(), key = operator.itemgetter(1))
idx = len(sorted_tis) * 2 / 3
top3d_tis = dict(sorted_tis[idx:])
sorted_sed = sorted(mean_sed.items(), key = operator.itemgetter(1))
idx = len(sorted_sed) * 2 / 3
top3d_sed = dict(sorted_sed[idx:])

danger_tis = dict()
for key in mean_tis:
    if mean_tis[key] >= DANGER:
        danger_tis[key] = mean_tis[key]
with open('danger_tis.csv', 'w') as wf:
    writer = csv.writer(wf, delimiter=',')
    for key in danger_tis:
        row = [coord_tis[key][0], coord_tis[key][1], danger_tis[key]]
        writer.writerow(row)

with open('top3d_tis.csv', 'w') as wf:
    writer = csv.writer(wf, delimiter=',')
    for key in top3d_tis:
        row = [coord_tis[key][0], coord_tis[key][1], top3d_tis[key]]
        writer.writerow(row)

with open('top3d_sed.csv', 'w') as wf:
    writer = csv.writer(wf, delimiter=',')
    for key in top3d_sed:
        row = [coord_sed[key][0], coord_sed[key][1], top3d_sed[key]]
        writer.writerow(row)

with open('cross.csv', 'w') as wf:
    writer = csv.writer(wf, delimiter=',')
    for sed_id in mean_sed:
        if sed_id in circles:
            for tis_id in circles[sed_id]:
                if tis_id in mean_tis:
                    row = [mean_sed[sed_id], mean_tis[tis_id]]
                    writer.writerow(row)
'''
