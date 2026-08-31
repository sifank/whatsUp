#!/usr/bin/python3
# Program: whatsUp
# Version: 20260706
# Author:  Sifan Kahale
# Desc:    Lists targets according to date/time, equipment, moon interference, etc.

from flask import Flask, render_template, request, redirect, url_for, flash
from flask_wtf import FlaskForm, CSRFProtect
import mysql.connector, os
import ephem, math
from datetime import datetime, timezone, UTC
import datetime
import math
import sys
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
from skyfield.api import wgs84, load, Star
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time

DEBUG = "True"

#------------------------------------------------------------------------------------------------------------
# Database Configuration
db_config = {
    'host': 'depoe',
    'user': 'sifan',
    'password': 'all4Sky',
    'database': 'kahaleobs'
}

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
    20: "Supernova "
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

# TODO change type back to 255
defSelection = {
    "database": "KBcatalog",
    "current_time_utc": datetime.datetime.now().astimezone(UTC),
    "TELESCOPE": "C14",
    "CAMERAS": "ASI294",
    "REDUCER": 1,
    "TYPE": 4,
    "MINALT": 20,
    "MOONDIST": 20,
    "FRAMESTAT": 3,
    "LATITUDE": 44.8825,
    "LONGITUDE": -124.0339,
    "ELEVATION": 15,
    "medReqLoc": "All",
    "medReqHr": 13,
    "visMag": 15,
    "priFlag": "All",
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

CRITERIA = defSelection

frameMap = ["All", "Too Large", "Too Small", "Perfect Fit"]

linelength = 145
app = Flask(__name__)
ID = 0

#--------------------------------------------------------
def calculate_sensor_dimensions(pixel_size_um, width_px, height_px):
    # Converts sensor pixel dimensions to physical dimensions in millimeters.
    width_mm = (pixel_size_um * width_px) / 1000.0
    height_mm = (pixel_size_um * height_px) / 1000.0
    return width_mm, height_mm

#--------------------------------------------------------
def calculate_fov(tele_fl, reducer_factor, sensor_width_mm, sensor_height_mm):
    #Calculates the Field of View (FOV) in degrees for a given setup.
    # Calculate effective focal length
    eff_fl = tele_fl * reducer_factor
    # Calculate FOV using the arctangent geometric formula
    fov_width_deg = 120 * math.atan(sensor_width_mm / (2 * eff_fl)) * (180 / math.pi)
    fov_height_deg = 120 * math.atan(sensor_height_mm / (2 * eff_fl)) * (180 / math.pi)
    return fov_width_deg, fov_height_deg

#--------------------------------------------------------
# Diag printout if DEBUG is true
def debug(msg):
    if DEBUG == "True":
        print(f"DIAG {msg}")
    return

#--------------------------------------------------------
# Opens the database
def openDB():
    try:
        conn = mysql.connector.connect(**db_config)
        cursor = conn.cursor()
    except mysql.connector.Error as err:
        print(f"Database error: {err}")
        print("Exiting ...")
        sys.exit(1)

    return conn, cursor

#--------------------------------------------------------
# Closes the databas
def closeDB(conn, cursor):
    conn.close()
    cursor.close()
    return

#------------------------------------------------------------------------------------------------------------
# Calculates the specific calendar date where target is at zenith at midnight.
# Earth circles the sun; midnight face aligns opposite to Sun location.
def calculate_midnight_zenith_date(ra_hours):
    target_sun_ra = (ra_hours + 12.0) % 24.0
    days_since_march_21 = target_sun_ra / (24.0 / 365.2422)
    equinox_base = datetime.date(2026, 3, 21)
    zenith_date = equinox_base + datetime.timedelta(days=int(days_since_march_21))
    return zenith_date.strftime("%b %d")


#------------------------------------------------------------------------------------------------------------
# calculate time from Meridian
def calculate_timeFromMeridian(objra, objdec,dt_utc, longitude):
    objCoord = SkyCoord(ra=objra* u.deg * 15, dec=objdec*u.deg)
    time_utc = Time(dt_utc, scale='utc')
    lst = time_utc.sidereal_time('mean', longitude=(longitude*u.deg).to(u.hourangle))
    hour_angle = (lst - objCoord.ra).wrap_at(12 * u.hour)
    deg_from_meridian = abs(hour_angle.degree)
    ha_hours = hour_angle.hour
    medHr = int(abs(ha_hours))
    medMin = int((abs(ha_hours) - medHr) * 60)
    medLoc = "After" if ha_hours > 0 else "Before"
    return(medLoc, medHr, medMin)

#------------------------------------------------------------------------------------------------------------
#Recommends base sub-exposure frame timings to capture details safely.
# Bright clusters require short exposures to preserve color; faint nebulae require longer exposures
def recommend_exposure(magnitude, object_type):
    if object_type == 3:  # Cluster
        base = 60 if magnitude > 6.0 else 30
    elif object_type ==  4:  # Globular
        base = 90 if magnitude > 7.0 else 60
    elif object_type in [6, 5, 15, 7]: # Planetary and other nebs
        base = 120 if magnitude > 9.0 else 60
    elif object_type in [8, 14]: # Galaxy
        base = 300 if magnitude > 8.5 else 180
    else:  # Emission / Reflection Nebulae
        base = 300 if magnitude > 6.5 else 180
    return f"{base}s"

#------------------------------------------------------------------------------------------------------------
# Determines observation conditions based on Sun altitude.
def get_sky_condition(sun_alt):
    if sun_alt >= 0:
        return "Daylight (Not viewable)"
    elif -6 <= sun_alt < 0:
        return "Civil Twilight (Bright sky)"
    elif -12 <= sun_alt < -6:
        return "Nautical Twilight (Sky darkening)"
    elif -18 <= sun_alt < -12:
        return "Astronomical Twilight (Deep sky targets becoming visible)"
    else:
        return "True Night (Ideal viewing conditions)"

#------------------------------------------------------------------------------------------------------------
# get gps from fieldallsky
def gpsFromAllsky():
    conn, cursor = openDB()
    sql = "SELECT * FROM gpsloc"
    try:
        cursor.execute(sql)
        row = cursor.fetchone()
    except mysql.connector.Error as err:
        print(f"Database error: {err}")
        sys.exit(1)

    GPSfix, CTime, UTime, LATITUDE, LONGITUDE, ELEVATION = row
    if GPSfix == "Invalid" or GPSfix == "NO_FIX":
        print(f"GPS issue: {GPSfix}")
        sys.exit(1)
    return(GPSfix, LATITUDE, LONGITUDE, ELEVATION )

#------------------------------------------------------------------------------------------------------------
# Evaluates how well a deep sky object frames within the sensor limits.
def evaluate_framing(target_size_arcmin, fov_w, fov_h):
    limiting_dimension = min(fov_w, fov_h)
    if target_size_arcmin > limiting_dimension:
        return "Too Large"
    elif target_size_arcmin < (limiting_dimension * 0.10):
        return "Too Small"
    else:
        return "Perfect Fit"

#------------------------------------------------------------------------------------------------------------
# Calc datetime obj for 11pm tonight
def timeAt10():
    local_now = datetime.datetime.now().astimezone()
    local_11pm = datetime.datetime.combine(local_now.date(), datetime.time(22, 0))
    local_11pm_aware = local_11pm.replace(tzinfo=local_now.tzinfo)
    current_time_utc = local_11pm_aware.astimezone(datetime.timezone.utc)
    return(current_time_utc)

#============ Creates visible obj list from selection conf  ==================
def createList(CRITERIA):
    # Initialize Skyfield Ephemeris
    ts = load.timescale()
    eph = load('de421.bsp')
    earth, sun, moon = eph['earth'], eph['sun'], eph['moon']
    #observer = earth + wgs84.latlon(CRITERIA["LONGITUDE"], elevation_m=CRITERIA["ELEVATION"])
    observer = earth + wgs84.latlon(CRITERIA["LATITUDE"], CRITERIA["LONGITUDE"], elevation_m=CRITERIA["ELEVATION"])

    # Time Parsing
    t = ts.from_datetime(CRITERIA["current_time_utc"])

    # Tracker System: Sun
    sun_pos = observer.at(t).observe(sun).apparent()
    sun_alt, _, _ = sun_pos.altaz()
    sky_status = get_sky_condition(sun_alt.degrees)

    # Tracker System: Moon Position & Illumination Phase
    moon_pos = observer.at(t).observe(moon).apparent()
    moon_alt, moon_az, _ = moon_pos.altaz()

    # Calculate phase/illumination using earth-centered sun/moon separation
    e_sun = earth.at(t).observe(sun)
    e_moon = earth.at(t).observe(moon)
    phase_angle = e_moon.phase_angle(sun).degrees
    moon_illumination = (1.0 + math.cos(math.radians(phase_angle))) / 2.0 * 100

    # --- Compute FOV --------------pixel_size, widthpx, heightpx
    camPx = CAMERAS[CRITERIA['CAMERAS']][0]
    camWx = CAMERAS[CRITERIA['CAMERAS']][1]
    camHx = CAMERAS[CRITERIA['CAMERAS']][2]
    teleFL = TELESCOPE[CRITERIA['TELESCOPE']]
    #debug(f"FL: {teleFL}, Px: {camPx}, Wx: {camWx}, Hx: {camHx}")
    sensor_w_mm, sensor_h_mm = calculate_sensor_dimensions(camPx, camWx, camHx)
    fov_width, fov_height = calculate_fov(teleFL, CRITERIA['REDUCER'], sensor_w_mm, sensor_h_mm)

    hdrInfo = {
        "time-utc": CRITERIA["current_time_utc"].strftime('%Y-%m-%d %H:%M:%S'),
        "time-loc": CRITERIA["current_time_utc"].astimezone().strftime('%Y-%m-%d %H:%M:%S'),
        "sun_alt": f"{sun_alt.degrees:.1f}°",
        "sky_status": sky_status,
        "moon_illumination": f"{moon_illumination:.1f}%",
        "moon_alt": f"{moon_alt.degrees:.1f}°",
        "moon_az": f"{moon_az.degrees:.1f}°",
        "fov_width": f"{fov_width:.1f}\"",
        "fov_height": f"{fov_height:.1f}\""
        }

    # Visual Canvas Instantiation
    visible_az, visible_alt, visible_data = [], [], [] # FIXME maybe don't need az/alt?

    #------------------ get objects from db --------------------------------------------------------
    conn, cursor = openDB()
    sql = f"SELECT * FROM {CRITERIA['database']}"  # TODO add selection for vMag, pri, and types to sql statement
    try:
        cursor.execute(sql)
        row = cursor.fetchone()

    except mysql.connector.Error as err:
        print(f"Database error: {err}")
        sys.exit(1)

    # for name, data in db:
    while row is not None:
        index, name, common_name, obj_type, ra, dec, mag, size, pa, image, flag = row
        #debug(f"Row {index} {name}")
        target = Star(ra_hours=ra, dec_degrees=dec)
        alt, az, _ = observer.at(t).observe(target).apparent().altaz()

        # TODO add visual mag selection
        # TODO add type selection
        # TODO what to do about PA?
        # ---- only include if above min altitude
        if alt.degrees.item() > CRITERIA["MINALT"]:
            # --- check framing stat, if size not set skip ------
            #debug(f"{name} meets minalt")
            #if size == "None" or size is None:
                #framing_status = "--"
                #psize = "--"
            #else:
            framing_status = evaluate_framing(float(size),   fov_width, fov_height)
            psize = f"{size:1f}"
            #debug(f"{name} framing: {framing_status}, request: {frameMap[CRITERIA['FRAMESTAT']]}")
            #if framing_status == CRITERIA["FRAMESTAT"] or CRITERIA["FRAMESTAT"] == "All" or framing_status == "--":
            if framing_status == frameMap[CRITERIA["FRAMESTAT"]] or frameMap[CRITERIA["FRAMESTAT"]] == "All" or framing_status == "--":
                debug(f"{name} meets framing")
                # ---- check for before or after median
                medLoc, medHr, medMin = calculate_timeFromMeridian(ra, dec, CRITERIA["current_time_utc"], CRITERIA["LONGITUDE"])
                if CRITERIA["medReqLoc"] == medLoc or CRITERIA["medReqLoc"] == 'All':

                    # ----- check hours before or after median
                    if medHr <= CRITERIA["medReqHr"]:
                        visible_az.append(np.radians(az.degrees)) # FIXME maybe don't need?'
                        visible_alt.append(90.0 - alt.degrees) # FIXME maybe don't need?'

                        # Calculate angular separation between target and moon
                        astrometric = observer.at(t).observe(target)
                        target_apparent = astrometric.apparent()
                        sep_deg = target_apparent.separation_from(moon_pos).degrees

                        # Determine capture rating dynamically
                        # Run Mathematical Engine Functions
                        midnight_zenith = calculate_midnight_zenith_date(ra)
                        medLoc, medHr, medMin = calculate_timeFromMeridian(ra, dec, CRITERIA["current_time_utc"], CRITERIA["LONGITUDE"])
                        timeFromMeridian = f"{medHr:>2d}:{medMin:>02d} hrs {medLoc:<6} "
                        # ---- skip if we don't have magnitude'
                        if not isinstance(mag, (int, float)):
                            recommended_sub = "--"
                            pmag = 0.0
                        else:
                            recommended_sub = recommend_exposure(float(mag), obj_type)
                            pmag = f"{float(mag):.1f}"

                        # --- and save for output to table and/or plot
                        visible_data.append({"name": name,
                                            "common_name": common_name,
                                            "alt": alt.degrees,
                                            "az": az.degrees,
                                            "size": psize,
                                            "mag": pmag,
                                            "type": obj_type,
                                            "exposure": recommended_sub,
                                            "framing": framing_status,
                                            "moon_sep": sep_deg,
                                            "meridian": timeFromMeridian,
                                            "zenith_date": midnight_zenith})

        row = cursor.fetchone()

    closeDB(conn, cursor)
    return visible_data, hdrInfo

# ===================================================
@app.route('/', methods=['GET', 'POST'])
def index():
    CT = [CRITERIA['current_time_utc'].astimezone().strftime("%Y-%m-%d"),  CRITERIA['current_time_utc'].astimezone().strftime("%H"),  CRITERIA['current_time_utc'].astimezone().strftime("%M")]

    return render_template('whatsUpMain.html', DS = CRITERIA, EQ = EQUIPFOV, FM = frameMap, CT = CT, DSO = DSOTYPE, TS = TELESCOPE, CM = CAMERAS)

#===============================================================
@app.route('/list', methods=['POST'])
def list():
        # TODO add visual mag selection
        # TODO add type selection
        # TODO what to do about PA?

    telescope = request.form.get('telescope')
    if telescope is None:
        telescope = CRITERIA["TELESCOPE"]
    CRITERIA["TELESCOPE"] = telescope

    reducer = request.form.get('reducer')
    if reducer is None:
        reducer = CRITERIA["REDUCER"]
    CRITERIA["REDUCER"] = float(reducer)

    camera = request.form.get('camera')
    if camera is None:
        camera = CRITERIA["CAMERAS"]
    CRITERIA["CAMERAS"] = camera

    framing = request.form.get('framing')
    if framing is None:
        framing = CRITERIA["FRAMESTAT"]
    CRITERIA["FRAMESTAT"] = framing

    medHr = request.form.get('medHr')
    if medHr is None:
        medHr =CRITERIA["medReqHr"]
    CRITERIA["medReqHr"] = int(medHr)

    medLoc = request.form.get('medLoc')
    if medLoc is None:
        medLoc =CRITERIA["medReqLoc"]
    CRITERIA["medReqLoc"] = medLoc

    moonDist = request.form.get('moonDist')
    if moonDist is None:
        moonDist = CRITERIA["MOONDIST"]
    CRITERIA["MOONDIST"] = int(moonDist)

    minAlt = request.form.get('minAlt')
    if minAlt is None:
        minAlt = CRITERIA["MINALT"]
    CRITERIA["MINALT"] = int(minAlt)

    visMag = request.form.get('visMag')
    if visMag is None:
        visMag = CRITERIA["visMag"]
    CRITERIA["visMag"] = int(visMag)

    priFlag = request.form.get('priFlag')
    if priFlag is None:
        priFlag = CRITERIA["priFlag"]
    CRITERIA["priFlag"] = priFlag

    defPlot = request.form.get('defPlot')
    if defPlot is None:
        defPlot = CRITERIA["defPlot"]
    CRITERIA["defPlot"] = defPlot

    database = request.form.get('database')
    if database is None:
        database = CRITERIA["database"]
    CRITERIA["database"] = database

    # TODO read in long, lat, elev from template
    CRITERIA["LATITUDE"] = 44.8825
    CRITERIA["LONGITUDE"] = -124.0339
    CRITERIA["ELEVATION"] = 15

    # ----- Manage date/time strings -----------------------------------
    date_str = request.form.get('user_date')      # e.g., "2026-07-10"
    hour_str = request.form.get('user_hour')      # e.g., "03"
    minute_str = request.form.get('user_minute')  # e.g., "30"
    # 2. Combine individual time fragments into a single string
    time_str = f"{hour_str}:{minute_str}" # e.g., "03:30 PM"
    # 3. Parse individual components into date and time objects
    parsed_date = datetime.datetime.strptime(date_str, '%Y-%m-%d').date()
    parsed_time = datetime.datetime.strptime(time_str, '%H:%M').time() # %H is for 12-hour clock
    # 4. Merge into a final datetime object
    local_now= datetime.datetime.combine(parsed_date, parsed_time).astimezone()
    #local_now = final_datetime.astimezone()
    current_time_utc = local_now.astimezone(UTC)
    CRITERIA["current_time_utc"] = current_time_utc

    visible_data, hdrInfo = createList(CRITERIA)
    return render_template('whatsUpList.html', DS = CRITERIA, VS = visible_data, HI = hdrInfo, DSO = DSOTYPE, FM = frameMap, TS = TELESCOPE, CM = CAMERAS)

#===============================================================
@app.route('/defaults', methods=['POST'])
def defaults():
    CRITERIA = defSelection
    CRITERIA['current_time_utc'] = datetime.datetime.now().astimezone(UTC)
    CT = [CRITERIA['current_time_utc'].astimezone().strftime("%Y-%m-%d"),  CRITERIA['current_time_utc'].astimezone().strftime("%H"),  CRITERIA['current_time_utc'].astimezone().strftime("%M")]
    return redirect(url_for('index'))

#===============================================================
@app.route('/atTenpm', methods=['POST'])
def atTenpm():
    CRITERIA['current_time_utc'] = timeAt10()
    CT = [CRITERIA['current_time_utc'].astimezone().strftime("%Y-%m-%d"),  CRITERIA['current_time_utc'].astimezone().strftime("%H"),  CRITERIA['current_time_utc'].astimezone().strftime("%M")]
    return redirect(url_for('index'))

#===============================================================
@app.route('/kfobs', methods=['POST'])
def kfobs():  # FIXME this is resetting things like vMag!
    GPSfix, CRITERIA['LATITUDE'], CRITERIA['LONGITUDE'], CRITERIA['ELEVATION']  = gpsFromAllsky()
    return redirect(url_for('index'))











#===============================================================
#===============================================================
#--------------------------------------------------------
@app.route('/details', methods=['POST'])
def details():
    ID = request.form.get('ID')
    sql = f"SELECT * FROM obsSession WHERE ID={ID};"

    conn, cursor = openDB()
    try:
        cursor.execute(sql)
        Session = cursor.fetchone()
        #return render_template('SQLerror.html', err=sql)
    except mysql.connector.Error as err:
        return render_template('SQLerror.html', err=err)
        
    if cursor.rowcount == 0:
        return render_template('DoesNotExist.html', ID = ID, Type = "Session")

    sql = f"SELECT * FROM obsTarget WHERE SessionID = {Session[1]};"
    try:
        cursor.execute(sql)
        print(f"Targets sql: {sql}")
        Targets = cursor.fetchall()
    except mysql.connector.Error as err:
        return render_template('SQLerror.html', err=err)
    closeDB(conn, cursor)

    dtDate = Session[2]
    #debug(f"dtDate before test: {dtDate}, Type: {type(dtDate)}")
    if dtDate is None:
        dtDate = datetime.strptime(Session[1],"%Y%m%d").astimezone(timezone.utc)

    #debug(f"dtDate @ render: {dtDate}, Type: {type(dtDate)}")
    envInfo = refreshData(dtDate)
    #debug(f"envinfo date: {dtDate}, {envInfo}")
    return render_template('/sessionDetails.html', ID = ID, Session = Session, Targets = Targets, envInfo = envInfo)

#--------------------------------------------------------
@app.route('/search', methods=['POST'])
def search():
    sType = request.form.get('stype')
    if sType == "Location": sType = "locName"
    sRegex = request.form.get('sregex')

    debug(f"You entered a search for {sType} looking for: {sRegex}")

    sql = f"SELECT * FROM obsSession WHERE {sType} RLIKE \'{sRegex}\';"
    debug(f"Search SQL: {sql}")

    conn, cursor = openDB()
    try:
        cursor.execute(sql)
        Session = cursor.fetchall()
    except mysql.connector.Error as err:
        return render_template('SQLerror.html', err=err)

    if cursor.rowcount == 0:
        return render_template('DoesNotExist.html', ID = ID, Type = "Session")
        # get the targets associated with this session

    return render_template('/sessionSearch.html', ID = ID, Session = Session)

#--------------------------------------------------------
@app.route('/add', methods=['POST'])  
def add():
    ctime = datetime.now().strftime("%Y%m%d")
    reqDate = datetime.now().astimezone(timezone.utc)
    envInfo = refreshData(reqDate)
    return render_template('whatsUpAdd.html', envInfo=envInfo, ctime=ctime, OTAItems = OTAItems, ReducerItems = ReducerItems, MountItems = MountItems, ControllerItems = ControllerItems, FocuserItems = FocuserItems, RotatorItems = RotatorItems, GuiderItems = GuiderItems,  FilterItems = FilterItems, CameraItems = CameraItems)
    
#--------------------------------------------------------
@app.route('/addok', methods=['POST'])
def addok():
    if request.form['action'] == 'Add Session':
        sql = f"INSERT INTO obsSession set \
        Date = \"{request.form.get('Date')}\", \
        dtDate = \"{request.form.get('dtDate')}\", \
        locName = \"{request.form.get('Name')}\", \
        locLong = \"{request.form.get('Longitude')}\", \
        locLat = \"{request.form.get('Latitude')}\", \
        locHeight = \"{request.form.get('Altitude')}\", \
        OTA = \"{request.form.get('OTA')}\", \
        Reducer = \"{request.form.get('Reducer')}\", \
        Mount = \"{request.form.get('Mount')}\", \
        Controller = \"{request.form.get('Controller')}\", \
        Focuser = \"{request.form.get('Focuser')}\", \
        Rotator = \"{request.form.get('Rotator')}\", \
        Guider = \"{request.form.get('Guider')}\", \
        Filter = \"{request.form.get('Filter')}\", \
        Camera = \"{request.form.get('Camera')}\", \
        Objective = \"{request.form.get('Objectives')}\", \
        Notes = \"{request.form.get('Notes')}\" \
        ;"
        conn, cursor = openDB()
        try:
            #debug(f"add sql stmt: {sql}")
            cursor.execute(sql)
            conn.commit()
        except mysql.connector.Error as err:
            return render_template('SQLerror.html', err=err)
        closeDB(conn, cursor)

    else:
        pass

    return redirect('/')

#--------------------------------------------------------
@app.route('/modify', methods=['POST'])
def modify():
    ID = request.form.get('ID')
    sql = f"SELECT * FROM obsSession WHERE ID={ID};"
    conn, cursor = openDB()
    try:
        cursor.execute(sql)
        Session = cursor.fetchone()
    except mysql.connector.Error as err:
        return render_template('SQLerror.html', err=err)
    closeDB(conn, cursor)
        
    if cursor.rowcount == 0:
        return render_template('DoesNotExist.html', ID = ID, Type = "Session")
    
    debug(f"\nModify for ID: {ID} OTA: {Session[9]} Notes: {Session[8]} Objective: {Session[7]}\n")
    return render_template('whatsUpModify.html', ID = ID, Session = Session, OTAItems = OTAItems, ReducerItems = ReducerItems, MountItems = MountItems, ControllerItems = ControllerItems, FocuserItems = FocuserItems, RotatorItems = RotatorItems, GuiderItems = GuiderItems,  FilterItems = FilterItems, CameraItems = CameraItems)

#--------------------------------------------------------
@app.route('/modifyok', methods=['POST'])
def modifyok():
    if request.form['action'] == 'Submit':
        ID = request.form.get('ID')
        sql = f"SELECT * FROM obsSession WHERE ID={ID};"
        conn, cursor = openDB()
        try:
            cursor.execute(sql)
            Session = cursor.fetchone()
        except mysql.connector.Error as err:
            return render_template('SQLerror.html', err=err)
        closeDB(conn, cursor)

        OTA = request.form.get('OTA')
        if OTA == None:
            OTA = Session[9]

        Reducer = request.form.get('Reducer')
        if Reducer == None:
            Reducer = Session[10]

        Mount = request.form.get('Mount')
        if Mount == None:
            Mount = Session[11]

        Focuser = request.form.get('Focuser')
        if Focuser == None:
            Focuser = Session[13]

        Controller = request.form.get('Controller')
        if Controller == None:
            Controller = Session[14]

        Rotator = request.form.get('Rotator')
        if Rotator == None:
            Rotator = Session[14]

        Guider = request.form.get('Guider')
        if Guider == None:
            Guider =Session[15]

        Filter = request.form.get('Filter')
        if Filter == None:
            Filter =Session[16]

        Camera = request.form.get('Camera')
        if Camera == None:
            Camera = Session[17]

        ID = request.form.get('ID')
        sql = f"UPDATE obsSession set \
        Date = \"{request.form.get('Date')}\", \
        dtDate = \"{request.form.get('dtDate')}\", \
        locName = \"{request.form.get('Name')}\", \
        locLong = \"{request.form.get('Longitude')}\", \
        locLat = \"{request.form.get('Latitude')}\", \
        locHeight = \"{request.form.get('Altitude')}\", \
        Objective = \"{request.form.get('Objectives')}\", \
        Notes = \"{request.form.get('Notes')}\", \
        OTA = \"{OTA}\", \
        Reducer = \"{Reducer}\", \
        Mount = \"{Mount}\", \
        Controller = \"{Controller}\", \
        Focuser = \"{Focuser}\", \
        Rotator = \"{Rotator}\", \
        Guider = \"{Guider}\", \
        Filter = \"{Filter}\", \
        Camera = \"{Camera}\" \
        WHERE ID=\"{ID}\";"
        debug(f"\nModifyOk form for ID: {ID} OTA: {request.form.get('OTA')} Notes: {request.form.get('Notes')} Objective: {request.form.get('Objectives')}\n\n")
        conn, cursor = openDB()
        try:
            cursor.execute(sql)
            conn.commit()
            #debug(f"\nsql statement: {sql}\n")
        except mysql.connector.Error as err:
            return render_template('SQLerror.html', err=err)
        closeDB(conn, cursor)
        
    else:
        pass
    
    return redirect('/')

#--------------------------------------------------------
@app.route('/remove', methods=['POST'])
def remove():
    ID = request.form.get('ID')
    sql = f"SELECT * FROM obsSession WHERE ID={ID};"
    conn, cursor = openDB()
    #debug(f"Delete sql: {sql}")
    try:
        cursor.execute(sql)
        Session = cursor.fetchone()
    except mysql.connector.Error as err:
        return render_template('SQLerror.html', err=err)
    closeDB(conn, cursor)
        
    if cursor.rowcount == 0:
        return render_template('DoesNotExist.html', ID = ID, Type = "Session")
    
    return render_template('whatsUpRemove.html', ID = ID, Session = Session)

#--------------------------------------------------------
@app.route('/removeok', methods=['POST'])
def removeok():
    ID = request.form.get('ID')
    sql = f"DELETE FROM obsSession WHERE ID={ID};"
    conn, cursor = openDB()
    try:
        cursor.execute(sql)
        conn.commit()
    except mysql.connector.Error as err:
        return render_template('SQLerror.html', err=err)
    closeDB(conn, cursor)

    return redirect('/')

#--------------------------------------------------------
@app.route('/doesNotExist', methods=['POST'])
def doesNotExit(message):
    return redirect('/displayMsg', message=message)

#--------------------------------------------------------
@app.route('/refresh', methods=['POST'])
def refresh():
    return redirect('/')

#--------------------------------------------------------
@app.route('/rtnmain', methods=['POST'])
def rtnmain():
    return redirect('<a href="http://depoe:5000"> </a>')
    #return redirect('/')

#--------------------------------------------------------
if __name__ == '__main__':
    app.run(host='0.0.0.0',port=5010, debug=True)

'''
CREATE TABLE `obsSession` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `Date` varchar(20) DEFAULT NULL,
  `locName` varchar(40) DEFAULT NULL,
  `locLong` double DEFAULT NULL,
  `locLat` double DEFAULT NULL,
  `locHeight` double DEFAULT NULL,
  `Objective` varchar(150) DEFAULT NULL,
  `Notes` varchar(100) DEFAULT NULL,
  `OTA` varchar(20) DEFAULT NULL,
  `Reducer` varchar(20) DEFAULT NULL,
  `Mount` varchar(20) DEFAULT NULL,
  `Controller` varchar(20) DEFAULT NULL,
  `Focuser` varchar(20) DEFAULT NULL,
  `Rotator` varchar(20) DEFAULT NULL,
  `Guider` varchar(25) DEFAULT NULL,
  `Filter` varchar(20) DEFAULT NULL,
  `Camera` varchar(25) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `Date` (`Date`)
) ENGINE=InnoDB
'''
