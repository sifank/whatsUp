#!/usr/bin/python3
# Program: whatsUp.py
# Version: 20260822
# Author:  Sifan Kahale
# Desc:    Lists targets according to date/time, equipment, moon interference, etc.

# FIXME image loc web field is too wide in searchDB and DBmain templates
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

# Equipment selector choices
EQUIPFOV = {
    1: ("C14 1 ASI294", 16.8, 20.4),
    2: ("C14 .7 ASI294", 24.1, 16.4),
    3: ('ES152 1 ASI2600', 81.6, 54.0),
    4: ("C14 .2 ASI294", 162, 111),
    5: ('Dwarf3 T', 174.0, 102.0),
    6: ('Dwarf3 W', 2700.0, 1500.0)
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
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)
import whatsUpDefs as defs

# ===================================================
@app.route('/', methods=['GET', 'POST'])
def index():
    #defs.debug(f"Inside index, telescope: { app.config['CRITERIA']['TELESCOPE']}")
    CT = [app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%Y-%m-%d"),  app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%H"),  app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%M")]

    #defs.debug(f"DF Type (in main): {app.config['CRITERIA']['TYPE']}, Types type: {type(app.config['CRITERIA']['TYPE'][0])}")

    return render_template('whatsUpMain.html', DS = app.config['CRITERIA'], EQ = EQUIPFOV, CT = CT, DSO = DSOTYPE, TS = TELESCOPE, CM = CAMERAS)

#===============================================================
@app.route('/showlist', methods=['POST'])
def showlist():
    #defs.debug(f"DefSel (in /showlist): {app.config['CRITERIA']['CAMERAS']}")

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

    database = request.form.get('database')
    if database is None:
        database = app.config['CRITERIA']["database"]
    app.config['CRITERIA']["database"] = database

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
    #defs.debug(f"Long: {app.config['CRITERIA']['LONGITUDE']}, Lat: {app.config['CRITERIA']['LATITUDE']}, Elev: {app.config['CRITERIA']['ELEVATION']}")

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
    #defs.debug(f"Type selected: {TYPESEL}")

    visible_data, hdrInfo = defs.createList(app, CAMERAS, TELESCOPE)

    funcSel = request.form.get('action')
    #defs.debug(f"Selected function: {funcSel}")
    if funcSel is None or funcSel == "update":
        return redirect('/')

    #defs.debug(f"VS = {visible_data}")

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
    app.config['CRITERIA']['current_time_utc'] = defs.timeAt10()
    CT = [app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%Y-%m-%d"),  app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%H"),  app.config['CRITERIA']['current_time_utc'].astimezone().strftime("%M")]
    #defs.debug(f"DefSel (in /atTenpm): {app.config['CRITERIA']['CAMERAS']}")
    return redirect('/')

#===============================================================
@app.route('/kfobs', methods=['POST'])
def kfobs():
    #defs.debug(f"Camera (in /kfobs before): {app.config['CRITERIA']['CAMERAS']}")
    GPSfix, app.config['CRITERIA']['LATITUDE'], app.config['CRITERIA']['LONGITUDE'], app.config['CRITERIA']['ELEVATION']  = defs.gpsFromAllsky()
    #defs.debug(f"Camera (in /kfobs after): {app.config['CRITERIA']['CAMERAS']}")
    return redirect('/')

#--------------------------------------------------------
@app.route('/doesNotExist', methods=['POST'])
def doesNotExit(message):
    return render_template('displayMsg.html', message=message)

#--------------------------------------------------------
@app.route('/refresh', methods=['POST'])
def refresh():
    return redirect('/')

#===============================================================
@app.route('/file', methods=['POST'])
def file():
    checked_ids = request.form.getlist('selected_ids')
    #defs.debug(f"returned list: {checked_ids}")
    # Early exit if the user submitted without checking anything
    if not checked_ids:
        return render_template('displayMsg.html', message = "No items selected")

    funcSel = request.form.get('action')
    if funcSel is None or funcSel == "output":
        return render_template('/whatsUpSelected.html', SEL = checked_ids)

    #-------------  else plot ------------------------------------------------------
    Selected = request.form.getlist('Selected')
    #defs.debug(f"targets selected: {Selected}")
    checked_ids = request.form.getlist('selected_ids')
    #defs.debug(f"CheckedIDs: {checked_ids}\n")
    # Early exit if the user submitted without checking anything
    if not checked_ids:
        return render_template('displayMsg.html', messasge = "No items selected")

    inner_string = Selected[0]
    celestial_list = ast.literal_eval(inner_string)
    visible_az, visible_alt, visible_data = [], [], []
    for item in celestial_list:
        if item['name'] in checked_ids:
            visible_data.append(item['name'])
            visible_az.append(item['az'])
            visible_alt.append(90.0 - item['alt'])

    #defs.debug(f"Selected items: {visible_data} {visible_az} {visible_alt}")
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
    status, Session = defs.getHilites(app)
    if not status:
        return render_template('SQLerror.html', err=hilites)

    defs.debug(f"Hilites: {hilites}")
    return render_template('whatsUpDBmain.html', Session=Session, DSO=DSOTYPE, SEARCHTYPE= SEARCHTYPE)

#================Manage DataBase ==============================
@app.route('/managedb', methods=['GET', 'POST'])
def managedb():
    conn, cursor = defs.openDB()
    if conn == "None":
        return render_template('SQLerror.html', err="Can not open Database")

    sql = f"SELECT * FROM {app.config['CRITERIA']['database']} limit 100;"
    try:
        #defs.debug(f"managedb sql: {sql}")
        cursor.execute(sql)
        Session = cursor.fetchall()
    except mysql.connector.Error as err:
        defs.closeDB(conn, cursor)
        return render_template('SQLerror.html', err=err)
    defs.closeDB(conn, cursor)

    return render_template('whatsUpDBmain.html', Session=Session, DSO=DSOTYPE, SEARCHTYPE= SEARCHTYPE, PF = PROCFLAG, SQL = sql)

#--------------------------------------------------------
@app.route('/updatePriority', methods=['POST'])
def updatePriority():
    Priority = request.form.get('Priority')
    row = request.form.get('Row')
    SQL = request.form.get('SQL')

    #defs.debug(f"Updating Row Index {row} to Priority: {Priority}")
    #defs.debug(f"===========================\n{SQL}")

    conn, cursor = defs.openDB()
    if conn == "None":
        return render_template('SQLerror.html', err="Can not open Database")

    #---- update prority on this target ------
    sql = f"UPDATE KBcatalog SET procFlag = '{Priority}' where ID = {row};"
    try:
        #defs.debug(f"update pri - update sql: {sql}")
        cursor.execute(sql)
        conn.commit()
    except mysql.connector.Error as err:
        defs.closeDB(conn, cursor)
        defs.debug(f"sql error: {err}")
        return render_template('SQLerror.html', err=err)

    # ---- reread results list using preivous search info -----
    try:
        #defs.debug(f"update pri - reread sql: {SQL}")
        cursor.execute(SQL)
        Session = cursor.fetchall()
    except mysql.connector.Error as err:
        defs.closeDB(conn, cursor)
        defs.debug(f"sql error: {err}")
        return render_template('SQLerror.html', err=err)

    defs.closeDB(conn, cursor)

    return render_template('whatsUpDBmain.html', Session=Session, DSO=DSOTYPE, SEARCHTYPE= SEARCHTYPE, SQL = SQL)

#--------------------------------------------------------
@app.route('/search', methods=['POST'])
def search():
    sType = request.form.get('stype')
    sRegex = request.form.get('sregex').replace("'", "\\'")
    #defs.debug(f"You entered a search for {sType} looking for: {sRegex}")

    if sType == "Name":
        sql = f"SELECT * FROM {app.config['CRITERIA']['database']} WHERE objectName RLIKE \'{sRegex}\' or commonName RLIKE \'{sRegex}\';"
    elif sType == "Type":
        sql = f"SELECT * FROM {app.config['CRITERIA']['database']} WHERE type RLIKE \'{sRegex}\';"
    else:
        sql = f"SELECT * FROM {app.config['CRITERIA']['database']} WHERE {sType} RLIKE \'{sRegex}\';"

    conn, cursor = defs.openDB()
    try:
        #defs.debug(f"Search SQL: {sql}")
        cursor.execute(sql)
        Session = cursor.fetchall()
    except mysql.connector.Error as err:
        defs.closeDB(conn, cursor)
        return render_template('SQLerror.html', err=err)

    if cursor.rowcount == 0:
        return render_template('SQLerror.html', err = f"'{sType}' not found in {app.config['CRITERIA']['database']}")
        # get the targets associated with this session

    return render_template('whatsUpDBmain.html', Session=Session, DSO=DSOTYPE, SEARCHTYPE= SEARCHTYPE, PF = PROCFLAG, SQL=sql)

#--------------------------------------------------------
@app.route('/add', methods=['POST'])  
def add():
    defObj = {
    "objectName":  "",
    "commonName":  "",
    "otype":  99,
    "rightascension":  0.0,
    "declination":  0.0,
    "magnitude":  0.0,
    "apSize":  0.0,
    "PA":  0.0,
    "imageLoc":  " ",
    "procFlag":  "Need"
    }
    return render_template('whatsUpAdd.html', DF=defObj, DSO=DSOTYPE, SEARCHTYPE= SEARCHTYPE, PF = PROCFLAG)
    
#--------------------------------------------------------
@app.route('/addok', methods=['POST'])
def addok():
    if request.form.get('action') == "Cancel":
        return redirect('/managedb')

    sql = f"INSERT INTO {app.config['CRITERIA']['database']} SET \
    objectName = \"{request.form.get('objectName')}\", \
    commonName = \"{request.form.get('commonName')}\", \
    type = \"{request.form.get('Type')}\", \
    rightascension = \"{request.form.get('ra')}\", \
    declination = \"{request.form.get('dec')}\", \
    magnitude = \"{request.form.get('mag')}\", \
    apSize = \"{request.form.get('Size')}\", \
    PA = \"{request.form.get('PA')}\", \
    imageLoc = \"{request.form.get('imageLoc')}\", \
    procFlag = \"{request.form.get('Priority')}\" \
    ;"
    conn, cursor = defs.openDB()

    try:
        #defs.debug(f"add sql: {sql}")
        cursor.execute(sql)
        conn.commit()
    except mysql.connector.Error as err:
        defs.closeDB(conn, cursor)
        return render_template('SQLerror.html', err=err)
    defs.closeDB(conn, cursor)

    return redirect('/managedb')

#---------------- fill in add/modify form with selected row ----------------------------------------
@app.route('/fillWithSelected', methods=['POST'])
def fillWithSelected():
    SEL = list(ast.literal_eval(request.form.get('SEL')))
    defs.debug(f"SEL as list: {SEL}")
    Session = session.get('onlineSearch', [])

    return render_template('whatsUpSearchDB.html', Session=Session, SEL = SEL, DSO=DSOTYPE, SEARCHTYPE= SEARCHTYPE, PF = PROCFLAG)

#--------------------------------------------------------
# TODO if imgLoc is url, download, place in spc folder and rewrite imgLoc to point to disk
@app.route('/objSearch', methods=['POST'])
def objSearch():
    newObject = request.form.get('lookupObject').replace("'", "\\'")
    OpType = request.form.get('OpType')
    #defs.debug(f"OpType in objSearch: {OpType}")
    if not newObject or newObject.strip() == '':
        return render_template('SQLerror.html', err=f"**{newObject}** not a valid search term")
    Session = []


    defs.debug(f"objSearch - ##{newObject}##")
    # ---- lookup in databases -------------------------------------------------------------------
    sql = f"SELECT 'KBcatalog' AS source_table, t1.* FROM KBcatalog t1 WHERE t1.objectname LIKE '%{newObject}%'  or t1.commonName = '%{newObject}%' \
        UNION ALL SELECT 'IPcatalog' AS source_table, t2.* FROM IPcatalog t2 WHERE t2.objectname LIKE '%{newObject}%'  or t2.commonName = '%{newObject}%' \
        UNION ALL SELECT 'DScatalog' AS source_table, t3.* FROM DScatalog t3 WHERE t3.objectname LIKE '%{newObject}%'  or t3.commonName = '%{newObject}%';"

    conn, cursor = defs.openDB()
    try:
        #defs.debug(f"search find in db sql: {sql}")
        cursor.execute(sql)
        Session = cursor.fetchall()
    except mysql.connector.Error as err:
        defs.closeDB(conn, cursor)
        defs.debug('SQLerror.html', err=err)
    defs.closeDB(conn, cursor)
    #defs.debug(f"Database results: {inSession}")
    if not Session:
        Session = []

    i = 1 #  set id to numbers starting with 1 (db starts above 300)
    # ---- lookup in Simbad ----------------------------------------------------
    simObj = defs.lookup_objectSimbad(newObject)
    if simObj[0] == "Error":
        defs.debug(f"Simbad lookup error: {simObj[1]}")
    else:
        simObj[1] = i
        #defs.debug(f"(whatsUp) Simbad results: {simObj}")
        i += 1
        Session.append(simObj)    # only add if found

    # ---- lookup in Telescopius ----------------------------------------------------------
    terror, numObjects, foundObj = defs.lookup_objectTelescopius(newObject)
    if terror != "Success":
        defs.debug(f"Telescopius error: {terror}")
    elif numObjects == 0:
        defs.debug(f"Telescopius: none found")
    else:
        for obj in foundObj:
            obj = obj[:1] + (i,) + obj[2:]  #setting id to be unique
            i += 1
            Session.append(obj)

    # ---- test if KBcatalog has this entry and fill edit form with it, if not, set to blank form
    #defs.debug(f"objSearch: Session is: {Session}")
    SEL = next((item for item in Session if item[0] == 'KBcatalog'), None)
    if SEL:
        defs.debug(f"Found entry: {SEL}")
    else:
        SEL = ["", 0, "", "", 99, 0.0, 0.0, 0.0, 0.0, 0.0, "", "Need"]

    session['onlineSearch'] = Session
    defs.debug(f"svdSession: {session.get('onlineSearch', [])}")
    return render_template('whatsUpSearchDB.html', Session=Session, SEL = SEL, DSO=DSOTYPE, SEARCHTYPE= SEARCHTYPE, PF = PROCFLAG, OT = OpType)

#--------------------------------------------------------
@app.route('/modify', methods=['POST'])
def modify():
    ID = request.form.get('ID')
    sql = f"SELECT * FROM {app.config['CRITERIA']['database']} WHERE ID={ID};"
    conn, cursor = defs.openDB()
    try:
        #defs.debug(f"Modify SQL: {sql}")
        cursor.execute(sql)
        Session = cursor.fetchone()
    except mysql.connector.Error as err:
        defs.closeDB(conn, cursor)
        return render_template('SQLerror.html', err=err)
    defs.closeDB(conn, cursor)
        
    if cursor.rowcount == 0:
        return render_template('SQLerror.html', err = f"{ID} in {app.config['CRITERIA']['database']} not found")
    
    #defs.debug(f"\nModify for ID: {ID} OTA: {Session[9]} Notes: {Session[8]} Objective: {Session[7]}\n")
    return render_template('whatsUpModify.html', ID = ID, SEL = Session, DSO=DSOTYPE, PF = PROCFLAG)

#--------------------------------------------------------
@app.route('/modifyok', methods=['POST'])
def modifyok():
    ID = request.form.get('ID')
    objectName = request.form.get('objectName')
    if objectName == None:  objectName = request.form.get('objectName_orig')
    commonName = request.form.get('commonName')
    if commonName is None:  commonName = request.form.get('commonName_orig')
    Type = request.form.get('Type')
    if Type is None:  Type = request.form.get('Type_orig')
    ra = request.form.get('ra')
    if ra == None:  ra = request.form.get('ra_orig')
    dec = request.form.get('dec')
    if dec == None:  dec = request.form.get('dec_orig')
    mag = request.form.get('mag')
    if mag == None:  mag = request.form.get('mag_orig')
    Size = request.form.get('Size')
    if Size == None:  Size =request.form.get('Size_orig')
    PA = request.form.get('PA')
    if PA == None:  PA =request.form.get('PA_orig')
    imageLoc = request.form.get('imageLoc')
    if imageLoc == None:  imageLoc = request.form.get('imageLoc_orig')
    Priority = request.form.get('Priority')
    if Priority == None:  Priority = request.form.get('Priority_orig')

    defs.debug(f"ID: {ID}: {objectName} {commonName} {Type} {ra} {dec} {mag} {Size} {PA} {imageLoc} {Priority}")

    sql = f"INSERT INTO {app.config['CRITERIA']['database']} \
            ( objectName, commonName, type, rightascension, declination, magnitude, apSize, PA, imageLoc, procFlag) \
        VALUES ( \
            \"{objectName}\", \
            \"{commonName}\", \
            \"{Type}\", \
            \"{ra}\", \
            \"{dec}\", \
            \"{mag}\", \
            \"{Size}\", \
            \"{PA}\", \
            \"{imageLoc}\", \
            \"{Priority}\" \
            ) \
        ON DUPLICATE KEY UPDATE \
            objectName = \"{objectName}\", \
            commonName = \"{commonName}\", \
            type = \"{Type}\", \
            rightascension = \"{ra}\", \
            declination = \"{dec}\", \
            magnitude = \"{mag}\", \
            apSize = \"{Size}\", \
            PA = \"{PA}\", \
            imageLoc = \"{imageLoc}\", \
            procFlag = \"{Priority}\";"

    conn, cursor = defs.openDB()
    try:
        #defs.debug(f"ModifyOK SQL {sql}")
        cursor.execute(sql)
        conn.commit()
    except mysql.connector.Error as err:
        defs.closeDB(conn, cursor)
        return render_template('SQLerror.html', err=err)
    defs.closeDB(conn, cursor)

    return redirect('/managedb')

#--------------------------------------------------------
@app.route('/removeok', methods=['POST'])
def removeok():
    ID = request.form.get('ID')
    sql = f"DELETE FROM {app.config['CRITERIA']['database']} WHERE ID={ID};"
    conn, cursor = defs.openDB()
    try:
        defs.debug(f"removeok SQL {sql}")
        cursor.execute(sql)
        conn.commit()
    except mysql.connector.Error as err:
        defs.closeDB(conn, cursor)
        return render_template('SQLerror.html', err=err)
    defs.closeDB(conn, cursor)

    return redirect('/managedb')

#--------------------------------------------------------
if __name__ == '__main__':
    app.run(host='0.0.0.0',port=5010, debug=True)

'''
CREATE TABLE `KBcatalog` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `objectName` varchar(50) DEFAULT NULL,
  `commonName` varchar(700) DEFAULT NULL,
  `type` int(11) DEFAULT NULL,
  `rightascension` float DEFAULT NULL,
  `declination` float DEFAULT NULL,
  `magnitude` float DEFAULT NULL,
  `apSize` float DEFAULT NULL,
  `PA` float DEFAULT NULL,
  `imageLoc` varchar(100) DEFAULT NULL,
  `procFlag` varchar(15) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `idx_objectName` (`objectName`)
) ENGINE=InnoDB
'''
