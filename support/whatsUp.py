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

DEBUG = "False"

#------------------------------------------------------------------------------------------------------------
# Database Configuration
db_config = {
    'host': 'depoe',
    'user': 'sifan',
    'password': 'all4Sky',
    'database': 'kahaleobs'
}

# Equipment selector choices
EQUIPFOV = {
    1: ("C14 .2 ASI294", 162, 111),
    2: ("C14 .7 ASI294", 24.1, 16.4),
    3: ("C14 1 ASI294", 16.8, 20.4),
    4: ('Dwarf3 W', 35.0, 35.0),
    5: ('Dwarf3 T', 2.03, 1.47),
    6: ('ES152 1 ASI2600', 81.6, 54.0)
}

frameMap = ["All", "Too Large", "Too Small", "Perfect Fit"]

linelength = 145
app = Flask(__name__)
ID = 0

#--------------------------------------------------------
# Diag printout if DEBUG is true
def debug(msg):
    if DEBUG == "True":
        print(f"DIAG {msg}")

    return

#--------------------------------------------------------
# Opens the whatsUp database
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
    if "Cluster" in object_type:
        base = 60 if magnitude > 6.0 else 30
    elif "Globular" in object_type:
        base = 90 if magnitude > 7.0 else 60
    elif "Planetary" in object_type:
        base = 120 if magnitude > 9.0 else 60
    elif "Galaxy" in object_type:
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
def timeAt11():
    local_now = datetime.datetime.now().astimezone()
    local_11pm = datetime.datetime.combine(local_now.date(), datetime.time(23, 0))
    local_11pm_aware = local_11pm.replace(tzinfo=local_now.tzinfo)
    current_time_utc = local_11pm_aware.astimezone(datetime.timezone.utc)
    return(current_time_utc)


#============ Creates visible obj list from selection conf  ==================
def createList(DEFSEL):
    # Initialize Skyfield Ephemeris
    ts = load.timescale()
    eph = load('de421.bsp')
    earth, sun, moon = eph['earth'], eph['sun'], eph['moon']
    #observer = earth + wgs84.latlon(DEFSEL["LONGITUDE"], elevation_m=DEFSEL["ELEVATION"])
    observer = earth + wgs84.latlon(DEFSEL["LATITUDE"], DEFSEL["LONGITUDE"], elevation_m=DEFSEL["ELEVATION"])

    # Time Parsing
    t = ts.from_datetime(DEFSEL["current_time_utc"])

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

    hdrInfo = {
        "time-utc": DEFSEL["current_time_utc"].strftime('%Y-%m-%d %H:%M:%S'),
        "time-loc": DEFSEL["current_time_utc"].astimezone().strftime('%Y-%m-%d %H:%M:%S'),
        "sun_alt": f"{sun_alt.degrees:.1f}°",
        "sky_status": sky_status,
        "moon_illumination": f"{moon_illumination:.1f}%",
        "moon_alt": f"{moon_alt.degrees:.1f}°",
        "moon_az": f" {moon_az.degrees:.1f}°"
        }

    # Visual Canvas Instantiation
    visible_az, visible_alt, visible_data = [], [], []

    #------------------ get objects from db --------------------------------------------------------
    conn, cursor = openDB()
    sql = "SELECT * FROM whatsUp"
    try:
        cursor.execute(sql)
        row = cursor.fetchone()

    except mysql.connector.Error as err:
        print(f"Database error: {err}")
        sys.exit(1)

    # for name, data in db:
    while row is not None:
        index, name, common_name, obj_type, ra, dec, mag, size, image, flag = row
        debug(f"Row {index} {name}")
        target = Star(ra_hours=ra, dec_degrees=dec)
        alt, az, _ = observer.at(t).observe(target).apparent().altaz()

        # TODO add visual mag selection
        # ---- only include if above min altitude
        #print(f"DIAG type of alt.degrees: {type(alt.degrees)}, minalt: {type(DEFSEL['MINALT'])}")
        if alt.degrees.item() > DEFSEL["MINALT"]:

            # --- check framing stat, if size not set skip ------
            if size == "--":
                framing_status = "--"
                psize = "--"
            else:
                framing_status = evaluate_framing(float(size),  EQUIPFOV[(DEFSEL['confNumber'])][1], EQUIPFOV[(DEFSEL['confNumber'])][2])
                #framing_status = evaluate_framing(float(size),  DEFSEL["fov_w_arcmin"], DEFSEL["fov_h_arcmin"])
                psize = f"{float(size):.1f}"
            if framing_status == DEFSEL["FRAMESTAT"] or DEFSEL["FRAMESTAT"] == "All" or framing_status == "--":

                # ---- check for before or after median
                medLoc, medHr, medMin = calculate_timeFromMeridian(ra, dec, DEFSEL["current_time_utc"], DEFSEL["LONGITUDE"])
                if DEFSEL["medReqLoc"] == medLoc or DEFSEL["medReqLoc"] == 'All':

                    # ----- check hours before or after median
                    if medHr <= DEFSEL["medReqHr"]:
                        visible_az.append(np.radians(az.degrees))
                        visible_alt.append(90.0 - alt.degrees)

                        # Calculate angular separation between target and moon
                        astrometric = observer.at(t).observe(target)
                        target_apparent = astrometric.apparent()
                        sep_deg = target_apparent.separation_from(moon_pos).degrees

                        # Determine capture rating dynamically
                        # Run Mathematical Engine Functions
                        midnight_zenith = calculate_midnight_zenith_date(ra)
                        medLoc, medHr, medMin = calculate_timeFromMeridian(ra, dec, DEFSEL["current_time_utc"], DEFSEL["LONGITUDE"])
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
    return visible_data, sun_alt, sky_status, moon_illumination, moon_alt, moon_az, hdrInfo

# ===================================================
@app.route('/', methods=['GET', 'POST'])
def index():
    #confName, fov_w_arcmin, fov_h_arcmin = EQUIPFOV["2"]
    #current_time_utc = timeAt11()

    # ----- extract current date/time for selection boxes
    now = datetime.datetime.now()
    current_date = now.strftime('%Y-%m-%d')          # e.g., "2026-07-10"
    current_hour = now.strftime('%H')                # 24-hour format with leading zero (e.g., "03")
    current_minute = now.minute
    currDateTime = [now.strftime('%Y-%m-%d'), now.strftime('%H'), str(now.minute)]
    current_time_utc = now.astimezone(UTC)
    #print(f"DIAG Data/time now: {now}")
    #print(f"DIAG Data/time utc now: {current_time_utc}")

    DEFSEL = {
        "current_time_utc": current_time_utc,
        "confNumber": 2,
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
        "defPlot": False
    }
    #print(f"DIAG:  medReqLoc: {DEFSEL['medReqLoc']}")

    visible_data, sun_alt, sky_status, moon_illumination, moon_alt, moon_az, hdrInfo = createList(DEFSEL)
    return render_template('whatsUpMain.html', DS = DEFSEL, HI = hdrInfo, EQ = EQUIPFOV, FM = frameMap, CT = currDateTime)

#===============================================================
@app.route('/list', methods=['POST'])
def list():
    confNumber = int(request.form.get('confNumber'))
    print(f"DIAG: raw confNumber returned: {confNumber}")
    if confNumber is None:
        confNumber = 2  # FIXME get this from defaults (need to set global)
    print(f"DIAG: equipment conf # returned: {EQUIPFOV[confNumber]}")

    framing = request.form.get('framing')
    if framing is None:
        framing = frameMap[3]  # FIXME get this from defaults (need to set global)
    print(f"DIAG: framing after: {framing}")

    medHr = request.form.get('medHr')
    if medHr is None:
        medHr =12  # FIXME get this from defaults (need to set global)
    print(f"DIAG: Median Hr after: {medHr}")

    medLoc = request.form.get('medLoc')
    if medLoc is None:
        medLoc ="All"  # FIXME get this from defaults (need to set global)
    print(f"DIAG: Median Loc after: {medLoc}")

    moonDist = request.form.get('moonDist')
    if moonDist is None:
        moonDist =20  # FIXME get this from defaults (need to set global)
    print(f"DIAG: Moon Distance after: {moonDist}")

    minAlt = request.form.get('minAlt')
    if minAlt is None:
        minAlt =20  # FIXME get this from defaults (need to set global)
    print(f"DIAG: Minimum Alt after: {minAlt}")

    # Manage date/time strings
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
    print(f"DIAG: date/time after: {local_now.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"DIAG: date/time utc: {current_time_utc.strftime('%Y-%m-%d %H:%M:%S')}")

    DEFSEL = {
        "current_time_utc": current_time_utc,
        "confNumber": confNumber,
        "MINALT": minAlt,
        "MOONDIST": moonDist,
        "FRAMESTAT": framing,
        "LATITUDE": 44.8825,
        "LONGITUDE": -124.0339,
        "ELEVATION": 15,
        "medReqLoc": medLoc,
        "medReqHr": medHr,
        "visMag": 15,
        "priFlag": "All",
        "defPlot": False
    }

    visible_data, sun_alt, sky_status, moon_illumination, moon_alt, moon_az, hdrInfo = createList(DEFSEL)
    return render_template('whatsUpList.html', DS = DEFSEL, VS = visible_data, HI = hdrInfo)

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
            #debug(f"DIAG add sql stmt: {sql}")
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
    #print(f"DIAG: Delete sql: {sql}")
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
