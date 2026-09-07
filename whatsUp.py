
#!/usr/bin/python3
# Program: whatsUp.py
# Version: 20260830
# Author:  Sifan Kahale
# Desc:    Lists targets according to date/time, equipment, moon interference, etc.

from flask import Flask, render_template, session, request, redirect, url_for, flash, current_app
from flask_wtf import FlaskForm, CSRFProtect
import mysql.connector, os
import ephem, math
from datetime import datetime, timezone, UTC
import math
import sys
import ast
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from skyfield.api import wgs84, load, Star
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
# import support functions:
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)
from  whatsUpDefs import *

#------------------------------------------------------------------------------------------------------------

# DSO type
DSOTYPE = {
    0: "Star",
    1: "Catalog Star",
    2: "Planet, Dwarf Planet, or Solar System body",
    3: "Open Cluster",
    4: "Globular Cluster",
    5: "Gaseous Nebula",
    6: "Planetary Nebula",
    7: "Supernova Remnant",
    8: "Galaxy",
    9: "Comet",
    10: "Asteroid",
    11: "Constellation",
    12: "Moon",
    13: "Asterism",
    14: "Galaxy Cluster",
    15: "Dark Nebula",
    16: "Quasar",
    17: "Multiple Star (Binary or Multi-star system)",
    18: "Radio Source",
    19: "Satellite",
    20: "Supernova ",
    99: "Unknown"
    }

# Common Name catalogs to accept (filter to)
okCat = ['ngc ', 'm ', 'name ', 'sh ', 'ic ', 'arp ']

# Search type for manageDB
SEARCHTYPE = {
    "Name": "Name",
    "Type": "type",
    "Magnitude": "magnitude",
    "Size": "apSize",
    "Image Loc": "imageLoc",
    "Priority": "procFlag"
    }

defSelection = {
    "database": "KBcatalog",
    "current_time_utc": datetime.now().astimezone(UTC),
    "TELESCOPE": "C14",
    "CAMERAS": "ASI294",
    "REDUCER": 1,
    "TYPE": ['255'],
    "MINALT": 20,
    "MOONDIST": 20,
    "FRAMESTAT": "Perfect Fit",
    "LATITUDE": 44.8825,
    "LONGITUDE": -124.0339,
    "ELEVATION": 15,
    "medReqHr": 4,
    "visMag": 15,
    "defPlot": "No"
}
# name: focal length
TELESCOPE = {
    "C14": 3910,
    "C11": 2800,
    "ES152": 988,
    "Dwarf3T": 150,
    "Dwarf3W": 6.7
    }

# name: (pixel_size, widthpx, heightpx)
CAMERAS = {
    "ASI294": (4.6, 4144, 2822),
    "ASI1600": (3.8, 4656, 3520),
    "ASI2600": (3.76, 6248, 4172),
    "Dwarf3T": (2.0, 3840, 2160),
    "Dwarf3W": (2.9, 1920, 1080)
    }

# priority types
PROCFLAG = [
    "Need",
    "Done",
    "Redo",
    "Should Redo"
    ]

# default object [source, recid, objectName, commonName, otype, ra, dec, mag, size, image, priority]
defObj = [
    " ",
    0,
    " ",
    " ",
    99,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    " ",
    "Need"
    ]

#ID = 0

app = Flask(__name__)
app.secret_key = 'KOBS Observatory'
app.config['CRITERIA'] = defSelection

# ===================================================
@app.route('/', methods=['GET', 'POST'])
def index():
    #debug(f"Inside index, telescope: { app.config['CRITERIA']['TELESCOPE']}")
    CT = [app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%Y-%m-%d"),  app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%H"),  app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%M")]

    #debug(f"DF Type (in main): {app.config['CRITERIA']['TYPE']}, Types type: {type(app.config['CRITERIA']['TYPE'][0])}")

    return render_template('whatsUpMain.html', DS = app.config['CRITERIA'], CT = CT, DSO = DSOTYPE, TS = TELESCOPE, CM = CAMERAS)

#===============================================================
@app.route('/showlist', methods=['POST'])
def showlist():
    #debug(f"DefSel (in /showlist): {app.config['CRITERIA']['CAMERAS']}")

    telescope = request.form.get('telescope')
    if telescope is None:
        telescope = app.config['CRITERIA']["TELESCOPE"]
    app.config['CRITERIA']["TELESCOPE"] = telescope

    reducer = request.form.get('reducer')
    if reducer is None:
        reducer = app.config['CRITERIA']["REDUCER"]
    app.config['CRITERIA']["REDUCER"] = float(reducer)

    camera = request.form.get('camera')
    if camera is None:
        camera = app.config['CRITERIA']["CAMERAS"]
    app.config['CRITERIA']["CAMERAS"] = camera

    framing = request.form.get('framing')
    if framing is None:
        framing = app.config['CRITERIA']["FRAMESTAT"]
    app.config['CRITERIA']["FRAMESTAT"] = framing

    medHr = request.form.get('medHr')
    if medHr is None:
        medHr =app.config['CRITERIA']["medReqHr"]
    app.config['CRITERIA']["medReqHr"] = int(medHr)

    moonDist = request.form.get('moonDist')
    if moonDist is None:
        moonDist = app.config['CRITERIA']["MOONDIST"]
    app.config['CRITERIA']["MOONDIST"] = int(moonDist)

    minAlt = request.form.get('minAlt')
    if minAlt is None:
        minAlt = app.config['CRITERIA']["MINALT"]
    app.config['CRITERIA']["MINALT"] = int(minAlt)

    visMag = request.form.get('visMag')
    if visMag is None:
        visMag = app.config['CRITERIA']["visMag"]
    app.config['CRITERIA']["visMag"] = int(visMag)

    app.config['CRITERIA']["database"] = 'KBcatalog'

    Longitude = request.form.get('Longitude')
    if Longitude is None:
        Longitude = app.config['CRITERIA']["LONGITUDE"]
    app.config['CRITERIA']["Longitude"] = Longitude

    Latitude = request.form.get('Latitude')
    if Latitude is None:
        Latitude = app.config['CRITERIA']["LATITUDE"]
    app.config['CRITERIA']["Latitude"] = Latitude

    Elevation = request.form.get('Elevation')
    if Elevation is None:
        Elevation = app.config['CRITERIA']["ELEVATION"]
    app.config['CRITERIA']["Elevation"] = Elevation
    #debug(f"Long: {app.config['CRITERIA']['LONGITUDE']}, Lat: {app.config['CRITERIA']['LATITUDE']}, Elev: {app.config['CRITERIA']['ELEVATION']}")

    # ----- Manage date/time strings -----------------------------------

    date_str = request.form.get('user_date')      # e.g., "2026-07-10"
    hour_str = request.form.get('user_hour')      # e.g., "03"
    minute_str = request.form.get('user_minute')  # e.g., "30"

    # 2. Combine individual time fragments into a single string
    time_str = f"{hour_str}:{minute_str}" # e.g., "03:30 PM"
    # 3. Parse individual components into date and time objects
    parsed_date = datetime.strptime(date_str, '%Y-%m-%d').date()
    parsed_time = datetime.strptime(time_str, '%H:%M').time() # %H is for 12-hour clock
    # 4. Merge into a final datetime object
    local_now= datetime.combine(parsed_date, parsed_time).astimezone()
    current_time_utc = local_now.astimezone(UTC)
    app.config['CRITERIA']["current_time_utc"] = current_time_utc

    # ----- Manage type selection -----------------------------------
    #TYPESEL = app.config['CRITERIA']["TYPE"]
    TYPESEL = []
    if request.form.get('All'): TYPESEL.append(request.form.get('All'))
    else:
        if request.form.get('Star'): TYPESEL.append(request.form.get('Star'))
        if request.form.get('Mstar'): TYPESEL.append(request.form.get('Mstar'))
        if request.form.get('OC'): TYPESEL.append(request.form.get('OC'))
        if request.form.get('GC'): TYPESEL.append(request.form.get('GC'))
        if request.form.get('PlanN'): TYPESEL.append(request.form.get('PlanN'))
        if request.form.get('Neb'): TYPESEL.append(request.form.get('Neb'))
        if request.form.get('DarkN'): TYPESEL.append(request.form.get('DarkN'))
        if request.form.get('SM'): TYPESEL.append(request.form.get('SM'))
        if request.form.get('Gal'): TYPESEL.append(request.form.get('Gal'))
        if request.form.get('GalC'): TYPESEL.append(request.form.get('GalC'))
        if request.form.get('Asterism'): TYPESEL.append(request.form.get('Asterism'))
        if request.form.get('PL'): TYPESEL.append(request.form.get('PL'))
        if request.form.get('Moon'): TYPESEL.append(request.form.get('Moon'))

    app.config['CRITERIA']['TYPE'] = TYPESEL
    #debug(f"Type selected: {TYPESEL}")

    visible_data, hdrInfo = createList(app, CAMERAS, TELESCOPE)

    funcSel = request.form.get('action')
    #debug(f"Selected function: {funcSel}")
    if funcSel is None or funcSel == "update":
        return redirect('/')

    #debug(f"VS = {visible_data}")

    return render_template('whatsUpList.html', DS = app.config['CRITERIA'], VS = visible_data, HI = hdrInfo, DSO = DSOTYPE, TS = TELESCOPE, CM = CAMERAS)

#===============================================================
@app.route('/defaults', methods=['POST'])
def defaults():
    app.config['CRITERIA'] = defSelection
    app.config['CRITERIA']['current_time_utc'] = datetime.now().astimezone(UTC)
    CT = [app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%Y-%m-%d"),  app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%H"),  app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%M")]
    return redirect('/')

#===============================================================
@app.route('/atTenpm', methods=['POST'])
def atTenpm():
    app.config['CRITERIA']['current_time_utc'] = timeAt10()
    CT = [app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%Y-%m-%d"),  app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%H"),  app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%M")]
    #debug(f"DefSel (in /atTenpm): {app.config['CRITERIA']['CAMERAS']}")
    return redirect('/')

#===============================================================
@app.route('/kfobs', methods=['POST'])
def kfobs():
    #debug(f"Camera (in /kfobs before): {app.config['CRITERIA']['CAMERAS']}")
    GPSfix, app.config['CRITERIA']['LATITUDE'], app.config['CRITERIA']['LONGITUDE'], app.config['CRITERIA']['ELEVATION']  = gpsFromAllsky()
    #debug(f"Camera (in /kfobs after): {app.config['CRITERIA']['CAMERAS']}")
    return redirect('/')

#--------------------------------------------------------
@app.route('/doesNotExist', methods=['POST'])
def doesNotExit(message):
    return render_template('displayMsg.html', message=message)

#--------------------------------------------------------
@app.route('/refresh', methods=['POST'])
def refresh():
    return redirect('/')

#--------------------------------------------------------
@app.route('/file', methods=['POST'])
def file():
    checked_ids = request.form.getlist('selected_ids')
    #debug(f"returned list: {checked_ids}")
    # Early exit if the user submitted without checking anything
    if not checked_ids:
        return render_template('displayMsg.html', message = "No items selected")

    funcSel = request.form.get('action')
    if funcSel is None or funcSel == "output":
        return render_template('/whatsUpSelected.html', SEL = checked_ids)

    elif funcSel =="cleanup":
        checked_ids = request.form.getlist('selected_ids')
        #debug("=" * 40)
        #debug(f"checked_ids = {checked_ids}")

        # get previous list
        fullList = request.form.getlist('Selected')
        inner_string = fullList[0]
        newList = ast.literal_eval(inner_string)
        #debug(f"newList = {newList}")

        # Retrieve hardware info
        HI= request.form.getlist('HI')
        debug(f"HI = {HI}")
        inner_string = HI[0]
        HI = ast.literal_eval(inner_string)

        # remove rows not in selected
        checked_set = set(checked_ids)
        filtered_list = [item for item in newList if item['name'] in checked_set]

        return render_template('whatsUpList.html', DS = app.config['CRITERIA'], VS = filtered_list, HI = HI, DSO = DSOTYPE, TS = TELESCOPE, CM = CAMERAS)

    #-------------  else plot ------------------------------------------------------
    Selected = request.form.getlist('Selected')
    #debug(f"targets selected: {Selected}")
    checked_ids = request.form.getlist('selected_ids')
    #debug(f"CheckedIDs: {checked_ids}\n")
    # Early exit if the user submitted without checking anything
    if not checked_ids:
        return render_template('displayMsg.html', messasge = "No items selected")

    inner_string = Selected[0]
    #debug(f"plot - selected targets: {Selected[0]}")
    celestial_list = ast.literal_eval(inner_string)
    visible_az, visible_alt, visible_data = [], [], []
    for item in celestial_list:
        if item['name'] in checked_ids:
            visible_data.append(item['name'])
            visible_az.append(item['az'])
            visible_alt.append(90.0 - item['alt'])

    #debug(f"Selected items: {visible_data} {visible_az} {visible_alt}")
    # ---- create plot
    output_dir = os.path.join(app.root_path, 'static')
    image_path = os.path.join(output_dir, 'whatsUpSkyChart.png')

    fig = plt.figure(figsize=(11, 11), facecolor='#0A0000')
    ax = fig.add_subplot(111, polar=True)
    ax.set_facecolor('#050000')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(1)
    for az_deg, alt_offset, obj in zip(visible_az, visible_alt, visible_data):
        az_rad = np.radians(az_deg)
        # plots text
        ax.text(az_rad, alt_offset+5, obj, color='#FF3333', fontsize=10, ha='center')
        # plots points
        ax.scatter(az_rad, alt_offset, color='#FF0000', s=60, marker="*")

    ax.set_ylim(0, 90)
    ax.set_yticklabels([])
    ax.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'], color='#AA0000', fontsize=11, weight='bold')
    ax.grid(color='#800000', linestyle='-', linewidth=0.9)
    plt.title(f"TACTICAL PHOTOGRAPHY PLANNER CHART\nLincoln Beach, OR | Night-Vision Active", color='#FF0000', fontsize=13, weight='bold', pad=25)
    plt.savefig(image_path, facecolor=fig.get_facecolor(), edgecolor='none', dpi=200, bbox_inches='tight')

    return render_template('skyChart.html')

#================Tonight's Hilites (from Telescopius)=================
@app.route('/hilites', methods=['GET', 'POST'])
def hilites():
    status, Session = getHilites(app)
    if not status:
        return render_template('SQLerror.html', err=hilites)

    debug(f"Hilites: {hilites}")
    return render_template('whatsUpDBmain.html', Session=Session, DSO=DSOTYPE, SEARCHTYPE= SEARCHTYPE, SQL='None')

#----------------Manage Lists ----------------------------------------
@app.route('/manageLists', methods=['POST'])
def manageLists():
    return render_template('displayMsg.html', message="Not Implemented Yet")

#----------------Redirect to Obs Sessions ----------------------------------------
@app.route('/Session', methods=['POST'])
def Session():
    return render_template('displayMsg.html', message="Not Implemented Yet")

#----------------Redirect to Obs Targets----------------------------------------
@app.route('/Target', methods=['POST'])
def Target():
    return render_template('displayMsg.html', message="Not Implemented Yet")

#--------------------------------------------------------
if __name__ == '__main__':
    app.run(host='0.0.0.0',port=5010, debug=True)

