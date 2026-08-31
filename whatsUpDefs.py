#!/usr/bin/python3
# Program: whatsUpDefs.py
# Version: 20260830
# Author:  Sifan Kahale
# Desc:    contains most of the def functions for whatsUp.py

DEBUG = True   # TODO Turn off when all done

import mysql.connector, os
from datetime import datetime, timezone, UTC, timedelta, time
import ephem, math, sys
import matplotlib.pyplot as plt
import numpy as np
from skyfield.api import wgs84, load, Star, Topos
from skyfield import almanac
from skyfield.magnitudelib import planetary_magnitude
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
import urllib.parse, socket
import requests
from telescopius import (
    TelescopiusClient,
    TelescopiusError,
    TelescopiusAuthError,
    TelescopiusBadRequestError,
    TelescopiusRateLimitError,
    TelescopiusNotFoundError,
    TelescopiusNetworkError,
)

# Common Name catalogs to accept (filter to)
okCat = ['ngc ', 'm ', 'name ', 'sh ', 'ic ', 'arp ']

# Client key for Telescopius
Telescopius_api_key="4c4a21528fcfc7141bca5b3d6e8b60a4"

# Map Simbad otype to KSTARS type
Telescopius_TO_KSTARS = {
    'star': 0,
    'dstar': 0,
    'mstar': 17,
    'gxy': 8,
    'eneb': 5,
    'rneb': 5,
    'dineb': 5,
    'pneb': 6,
    'snr': 7,
    'gcl': 4,
    'glcl': 4,
    'ocl': 3,
    'opcl': 3,
    'stcl': 3,
    'planet': 2,
    'asteroid': 10,
    'comet': 9
    }

ttypes = ['star', 'dstar', 'mstar', 'gxy', 'eneb', 'rneb', 'pneb', 'snr', 'gcl', 'glcl', 'ocl', 'opcl', 'stcl', 'planet', 'asteroid', 'comet']

# Map Simbad otype to KSTARS type
SIMBAD_TO_KSTARS = {
    # --- 0: Stars ---
    "Star": 0, "*": 0, "V*": 0, "Variable*": 0, "Double*": 0, "PM*": 0,
    "YSO": 0, "Be*": 0, "Carbon*": 0, "WhiteDwarf": 0, "Neutron*": 0,
    "SB*": 0,  "BY*": 0, "Em*": 0, "WR*": 0, "Be*": 0, "pA*": 0, "reg": 1,
    "RS*": 0, "s*b": 0, "Y*O": 0, "TT*": 0, "Or*": 0, "RR*": 0,

    # --- 3: Open Clusters ---
    "OpenCluster": 3, "Cl*": 3, "Assoc*": 3, "PartofCl": 3,  "As*": 3,
    "OpC": 3, "MGr": 3, "*cl": 3,  "*Cl": 3,

     # --- 4: Globular Clusters ---
    "GlobCluster": 4, "GlbCl": 4, "GlC": 4,

    # --- 5: Gaseous Nebulae ---
    "Nebula": 5, "EmNebula": 5, "ReflNebula": 5, "DarkNebula": 5,
    "HII": 5, "Cloud": 5, "RfNeb": 5, "AsyNeb": 5,  "Neb": 5,  "PlNeb": 6,
    "DrkNeb": 15, "MolCld": 5, "RNe": 5, "DNe": 5, "ISM": 5,

    # --- 6: Planetary Nebulae ---
    "PlanetaryNeb": 6, "PN": 6, "PND": 6, "PN?": 6,

    # --- 7: Supernova Remnants ---
    "SNRemnant": 7, "SNR": 7,

    # --- 8: Galaxies ---
    "Galaxy": 8, "G": 8, "BClG": 8, "Seyfert": 8, "Seyfert_1": 8,
    "Seyfert_2": 8, "HII_G": 8, "InteractG": 8, "PairG": 8, "GroupG": 8,
    "GTrpl": 8, "GlClG": 14, "GGroup": 14, "GClstr": 14, "AGN": 8,
    "ClG": 8, "EmG": 8, "GrG": 8, "IG": 8, "CGC": 8, "Sy1": 8, "GiG": 8,
    "CGG": 8,

    # --- 18: Radio Sources ---
    "Radio": 18, "RadioS": 18, "Maser": 18,

    # Other / Fallback
    "QSO": 16, "PN": 6, "Unknown": 99
}

#--------------------------------------------
# Database Configuration
db_config = {
    'host': 'depoe',
    'user': 'sifan',
    'password': 'all4Sky',
    'database': 'kahaleobs'
}

#--------------------------------------------
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
    if DEBUG:
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
    #equinox_base = datetime.date(2026, 3, 21).date()
    equinox_base = datetime(2026, 3, 21, tzinfo=UTC).date()
    zenith_date = equinox_base + timedelta(days=int(days_since_march_21))
    return zenith_date.strftime("%b %d")

#------------------------------------------------------------------------------------------------------------
# calculate time from Meridian
def calculate_timeFromMeridian(ra_deg, longitude_deg, dt_utc):
    """Calculates Hour Angle (HA) and returns it formatted as hours and minutes."""
    # Extract UTC date and time components
    year, month, day = dt_utc.year, dt_utc.month, dt_utc.day
    hour = dt_utc.hour + dt_utc.minute / 60.0 + dt_utc.second / 3600.0
    if month <= 2:
        year -= 1
        month += 12
    # Calculate Julian Date
    A = int(year / 100)
    B = 2 - A + int(A / 4)
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5 + hour / 24.0
    # Centuries from J2000.0
    d = jd - 2451545.0
    t = d / 36525.0
    # Greenwich Sidereal Time (GMST) in degrees
    gmst = (280.46061837 + 360.98564736629 * d + 0.000387933 * t**2 - (t**3) / 38710000.0) % 360.0
    # Local Sidereal Time (LST) in degrees
    lst_deg = (gmst + longitude_deg) % 360.0
    # Hour Angle (HA = LST - RA) in degrees
    ha_deg = (lst_deg - ra_deg) % 360.0
    # Normalize HA to range from -180 to +180 degrees
    if ha_deg > 180.0:
        ha_deg -= 360.0
    # Convert HA from degrees to total hours (15 degrees = 1 hour)
    ha_hours_total = ha_deg / 15.0
    # Extract hours and minutes
    is_negative = ha_hours_total < 0
    abs_hours = abs(ha_hours_total)
    hours = int(abs_hours)
    minutes = int(round((abs_hours - hours) * 60))
    # Handle minute overflow during rounding (e.g., 60m becomes +1h)
    if minutes == 60:
        hours += 1
        minutes = 0
    # Apply negative sign back if needed
    if is_negative:
        hours = -hours
        if hours == 0:
            return (-0, minutes)

    return(hours, minutes)

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
        return render_template('SQLerror.html', err=err)
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
    local_now = datetime.now().astimezone()
    local_11pm = datetime.combine(local_now.date(), time(22, 0))
    local_11pm_aware = local_11pm.replace(tzinfo=local_now.tzinfo)
    current_time_utc = local_11pm_aware.astimezone(timezone.utc)
    return(current_time_utc)

#------------------------------------------------------------------------------------------------------------
# find new target using Telescopius
def lookup_objectTelescopius(newObject):
    client = TelescopiusClient(api_key=Telescopius_api_key)
    terror = "Success"
    numObjects = 0

    try:
        results = client.search_targets(
            name=newObject,
            name_exact=True,
            lat=44.8825,
            lon=-124.0339,
            timezone="America/Los_Angeles",
            #types="GXY,ENEB",
            #min_alt=30,
            #mag_max=18,
            #datetime="2026-08-09 21:00:00",
            #results_per_page=120
        )
    except TelescopiusAuthError:
        terror = "Invalid API key"
        return(terror, 0, [])
    except TelescopiusRateLimitError:
        terror = "Rate limit exceeded - please wait before retrying"
        return(terror, 0, [])
    except TelescopiusBadRequestError as e:
        terror = f"Invalid parameters: {e.message}"
        return(terror, 0, [])
    except TelescopiusNetworkError:
        terror = "Network error - check your connection"
        return(terror, 0, [])
    except TelescopiusError as e:
        terror = f"API error: {e.message}"
        return(terror, 0, [])

    numObjects = results['matched']
    Session = []

    for item in results['page_results']:
        obj = item['object']
        # ['main_id', 'main_name', 'ids', 'names', 'family', 'types', 'url', 'main_image_url', 'thumbnail_url', 'alt_ids', 'ra', 'dec', 'con', 'con_name', 'visual_mag']
        obj.get('alt_ids').append(obj.get('main_name') or obj('main_id'))
        #debug(f"Telescopius type: {obj.get('types')}")
        coord = SkyCoord(ra=obj.get('ra'), dec=0, unit=(u.hourangle, u.deg))
        try:
            otype = Telescopius_TO_KSTARS[[item for item in obj.get('types') if item in set(ttypes)][0]]
        except:
            otype = f"**{obj.get('types')}**"

        target = (
            'Telescopius',
            0,
            newObject,
            ', '.join(obj.get('alt_ids')),
            otype,
            round(coord.ra.deg, 4),
            round(obj.get('dec'), 4),
            obj.get('visual_mag'),
            round(obj.get('major_axis')/60, 1),
            0.0,
            obj.get('thumbnail_url'),
            "Need"
            )

        Session.append(target)

    return(terror, numObjects, Session)

    # ====================================

#------------------------------------------------------------------------------------------------------------
# find new target using the Sesame online DB
def lookup_objectSimbad(object_name):
    """Looks up an object using the CDS Sesame Resolver with proper headers and mirrors."""
    encoded_name = urllib.parse.quote(object_name)
    # Alternative, high-availability mirror to bypass the connection block
    url = f"https://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-oI/A?{encoded_name}"
    # Headers fake a regular web browser profile so the server doesn't reject you
    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
        )
    }

    try:
        # Pass the safety headers explicitly
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
        lines = response.text.splitlines()

        simObj = []
        ra = None
        dec = None
        otype = "Unknown"
        commonName = []

        decFound = False
        for line in lines:
            #debug(f"Line: {line}")
            if line.startswith("%J "):
                parts = line.split()
                if len(parts) >= 3:
                    ra = parts[1]
                    if not decFound:
                        dec = parts[2]
                        decFound = True

            elif line.startswith("%C"):
                otype = line.replace("%C", "").strip().removeprefix('.0 ')
                try:
                    otype = SIMBAD_TO_KSTARS[otype]
                except:
                    otype = f"##{otype}##"
            elif line.startswith("%I "):
                canidate = (line.removeprefix('%I ').removeprefix('NAME ')).replace("[", "").replace("]", "")
                #if any(sub.lower() in canidate.lower() for sub in okCat): commonName.append(canidate)
                commonName.append(canidate)

        if ra is not None and dec is not None:
            # simObj:  [source, recid, objectName, commonName, otype, ra, dec, mag, size, image, priority]
            simObj = [
                "Simbad",
                0,
                object_name,
                ", ".join(commonName),
                otype,
                ra,
                dec,
                0.0,
                0.0,
                0.0,
                " ",
                "Need"
            ]
            #debug(f"(defs) Found: {simObj}")
            return simObj

        return ["Error", "No RA/Dec"]

    except requests.exceptions.RequestException as e:
        return ["Error", f"HTTP Request failed: {e}"]
    except Exception as e:
        return ["Error", f"Parsing failed: {e}"]

#============ Creates visible obj list from selection conf  ==================
def createList(app, CAMERAS, TELESCOPE):
    # Initialize Skyfield Ephemeris
    ts = load.timescale()
    base_dir = os.path.abspath(os.path.dirname(__file__))
    bsp_path = os.path.join(base_dir, 'de421.bsp')
    #eph = load('de421.bsp')
    eph = load(bsp_path)
    earth, sun, moon = eph['earth'], eph['sun'], eph['moon']
    #observer = earth + wgs84.latlon(app.config['CRITERIA']["LONGITUDE"], elevation_m=app.config['CRITERIA["ELEVATION"])
    observer = earth + wgs84.latlon(app.config['CRITERIA']["LATITUDE"], app.config['CRITERIA']["LONGITUDE"], elevation_m=app.config['CRITERIA']["ELEVATION"])

    # Time Parsing
    t = ts.from_datetime(app.config['CRITERIA']["current_time_utc"])

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
    camPx = CAMERAS[app.config['CRITERIA']['CAMERAS']][0]
    camWx = CAMERAS[app.config['CRITERIA']['CAMERAS']][1]
    camHx = CAMERAS[app.config['CRITERIA']['CAMERAS']][2]
    teleFL = TELESCOPE[app.config['CRITERIA']['TELESCOPE']]
    #debug(f"FL: {teleFL}, Px: {camPx}, Wx: {camWx}, Hx: {camHx}")
    sensor_w_mm, sensor_h_mm = calculate_sensor_dimensions(camPx, camWx, camHx)
    fov_width, fov_height = calculate_fov(teleFL, app.config['CRITERIA']['REDUCER'], sensor_w_mm, sensor_h_mm)

    hdrInfo = {
        "time-utc": app.config['CRITERIA']["current_time_utc"].strftime('%Y-%m-%d %H:%M:%S'),
        "time-loc": app.config['CRITERIA']["current_time_utc"].astimezone().strftime('%Y-%m-%d %H:%M:%S'),
        "sun_alt": f"{sun_alt.degrees:.1f}°",
        "sky_status": sky_status,
        "moon_illumination": f"{moon_illumination:.1f}%",
        "moon_alt": f"{moon_alt.degrees:.1f}°",
        "moon_az": f"{moon_az.degrees:.1f}°",
        "fov_width": f"{fov_width:.1f}\"",
        "fov_height": f"{fov_height:.1f}\""
        }

    # Visual Canvas Instantiation
    visible_data =  []

    #------------------ get objects from db --------------------------------------------------------
    conn, cursor = openDB()
    #debug(f"Type is: { int(app.config['CRITERIA']['TYPE'][0]) }")
    if int(app.config['CRITERIA']['TYPE'][0]) == 255:
        sql = f"SELECT * FROM {app.config['CRITERIA']['database']}"
    else:
        typesel = ",".join(app.config['CRITERIA']['TYPE'])
        #debug(f"sql type selection: {typesel}")
        sql = f"SELECT * FROM {app.config['CRITERIA']['database']} where type in ({typesel})"

    try:
        #debug(f"sql: {sql}")
        cursor.execute(sql)
        row = cursor.fetchone()

    except mysql.connector.Error as err:
        return render_template('SQLerror.html', err=err)

    # for name, data in db:
    while row is not None:
        index, name, common_name, obj_type, ra, dec, mag, size, pa, image, flag = row
        #debug(f"Row {index} {name}")
        target = Star(ra_hours=ra / 15, dec_degrees=dec)
        alt, az, _ = observer.at(t).observe(target).apparent().altaz()

        #  ----------- only select targets brighter than visMag -----------------------------
        if  mag <= app.config['CRITERIA']["visMag"]:

            # ---- only include if above min altitude --------------------------------------------
            if alt.degrees.item() > app.config['CRITERIA']["MINALT"]:
                framing_status = evaluate_framing(float(size),   fov_width, fov_height)
                psize = round(size, 2)

                # ----- only include targets matching framing criteria --------------------
                if framing_status == app.config['CRITERIA']["FRAMESTAT"] or app.config['CRITERIA']["FRAMESTAT"] == "All" or obj_type <= 1 :

                    # ---- only include targets within xhrs of meridian ----------------------
                    medHr, medMin = calculate_timeFromMeridian(ra, app.config['CRITERIA']["LONGITUDE"], app.config['CRITERIA']["current_time_utc"])
                    if abs(medHr + medMin/60) <= app.config['CRITERIA']["medReqHr"]:
                        #debug(f"made it: {name}")

                        # ----Calculate angular separation between target and moon-------
                        astrometric = observer.at(t).observe(target)
                        target_apparent = astrometric.apparent()
                        sep_deg = target_apparent.separation_from(moon_pos).degrees
                        if app.config['CRITERIA']['MOONDIST'] <= sep_deg:

                            # ------- Determine exposure -----------------------------------------------------
                            midnight_zenith = calculate_midnight_zenith_date(ra/15)
                            timeFromMeridian = f"{medHr:>2d}:{medMin:>02d}"
                            # ---- skip exp determination if we don't have magnitude'----------
                            if not isinstance(mag, (int, float)):
                                recommended_sub = "--"
                                pmag = 0.0
                            else:
                                recommended_sub = recommend_exposure(float(mag), obj_type)
                                pmag = f"{float(mag):.1f}"

                            # --- and save for output to table and/or plot ---------------------------
                            visible_data.append({"name": name,
                                                "common_name": common_name,
                                                "alt": alt.degrees,
                                                "az": az.degrees,
                                                "size": psize,
                                                "mag": pmag,
                                                "type": obj_type,
                                                "exposure": recommended_sub,
                                                "framing": framing_status,
                                                "PA": pa,
                                                "moon_sep": sep_deg,
                                                "priority": flag,
                                                "meridian": timeFromMeridian,
                                                "zenith_date": midnight_zenith,
                                                "image": image})

        row = cursor.fetchone()

    closeDB(conn, cursor)

    # ------------ add planets ------------------------------------------------
    if  '2' in app.config['CRITERIA']["TYPE"]:
        for planet_name in ['mercury barycenter', 'venus barycenter', 'mars barycenter', 'jupiter barycenter', 'saturn barycenter', 'uranus barycenter', 'neptune barycenter']:
            planet = eph[planet_name]

            # Compute apparent position and alt/az
            astrometric = observer.at(t).observe(planet)
            apparent = astrometric.apparent()
            alt, az, distance = apparent.altaz()

            # Apparent magnitude
            mag = planetary_magnitude(astrometric)

            # Apparent angular size in arcseconds
            # (approximate angular diameter formula: size_at_1_au / distance_in_au)
            # Standard radii values in km or angular size routines can be applied here
            au_km = 149597870.7
            planet_radius_km = {
                'mercury barycenter': 2439.7, 'venus barycenter': 6051.8, 'mars barycenter': 3389.5, 'jupiter barycenter': 69911.0, 'saturn barycenter': 58232.0, 'uranus barycenter': 25362.0, 'neptune barycenter': 24622.0
            }[planet_name]
            angular_size_arcsec = 2.0 * np.arctan(planet_radius_km / (distance.au * au_km)) * 206264.806

            # Hours before or after meridian (Hour Angle calculation)
            # Local Sidereal Time minus Right Ascension gives Hour Angle (HA)
            _, ra, _ = apparent.radec()
            lst = observer.to_scalars(t).apparent().radec()[0] if hasattr(observer, 'to_scalars') else apparent.hadec()[0] # simplified HA extraction
            ha, _, _ = observer.at(t).observe(planet).apparent().hadec()
            hours_from_meridian = ha.hours
            alt_decimal = alt.degrees
            az_decimal = az.degrees

            # ----- apply selection criteria -------------------------------------------------
            #  ----------- only select planets brighter than visMag -----------------------------
            if  mag <= app.config['CRITERIA']["visMag"]:

                # ---- only include if above min altitude --------------------------------------------
                if alt.degrees.item() > app.config['CRITERIA']["MINALT"]:
                    framing_status = evaluate_framing(float(angular_size_arcsec),   fov_width, fov_height)
                    psize = round(angular_size_arcsec, 2)

                    # ---- only include planets within xhrs of meridian ----------------------
                    if (ha.hours) <= app.config['CRITERIA']["medReqHr"]:

                        # ----Calculate angular separation between planet and moon-------
                        astrometric = observer.at(t).observe(planet)
                        planet_apparent = astrometric.apparent()
                        sep_deg = planet_apparent.separation_from(moon_pos).degrees

                        # ------- Determine exposure -----------------------------------------------------
                        recommended_sub = recommend_exposure(float(mag), 2)

                        #debug(f"Planet: {planet_name.split()[0].capitalize()} Alt: {alt_decimal:.4f}° Az: {az_decimal:.4f}° appMag: {mag:.2f} appSize: {angular_size_arcsec:.2f} arcseconds Meridian: {ha.hours:.2f}hrs")


                        # --- and save for output to table and/or plot ---------------------------
                        visible_data.append({"name": planet_name.split()[0].capitalize(),
                                            "common_name": planet_name.split()[0].capitalize(),
                                            "alt": alt_decimal,
                                            "az": az_decimal,
                                            "size": psize,
                                            "mag": mag.round(2),
                                            "type": 2,
                                            "exposure": recommended_sub,
                                            "framing": framing_status,
                                            "PA": 0,
                                            "moon_sep": round(sep_deg, 2),
                                            "priority": "Need",
                                            "meridian": round(ha.hours, 2),
                                            "zenith_date": 0})

    # ------------   moon ------------------------------------------------
    #debug(f"Image type: {app.config['CRITERIA']['TYPE']}")
    if  '12' in app.config['CRITERIA']["TYPE"]:
        # Define observer location (Lincoln Beach, OR)
        observer_loc = Topos(latitude_degrees=44.8105, longitude_degrees=-124.0418)
        # Dynamically calculate UTC time for 10 PM local tonight
        dt_utc = timeAt10()
        t = ts.from_datetime(dt_utc)
        # Setup bodies
        earth, moon, sun = eph['earth'], eph['moon'], eph['sun']
        observer = earth + observer_loc
        apparent_pos = observer.at(t).observe(moon).apparent()
        # 1. Altitude and Azimuth (Decimal degrees)
        alt, az, distance = apparent_pos.altaz()
        alt_deg = alt.degrees
        az_deg = az.degrees
        # 2. Angular Size (Apparent diameter in arcminutes)
        moon_radius_km = 1737.4
        angular_diameter_rad = 2 * np.arctan(moon_radius_km / distance.km)
        angular_diameter_arcmin = np.degrees(angular_diameter_rad) * 60
        # 3. Magnitude and Illuminated Fraction
        fraction = almanac.fraction_illuminated(eph, 'moon', t)
        phase_angle = moon.at(t).phase_angle(sun).degrees
        approx_mag = -12.74 + 0.026 * abs(phase_angle) + 4e-4 * (phase_angle ** 3)
        # 4. Hours from Meridian (Hour Angle)
        ha, dec, _ = apparent_pos.hadec()
        ha_hours = ha.hours
        # 5. Find next Meridian Transit — FIXED FOR LOCAL TIME WITHOUT SECONDS
        t0 = ts.utc(dt_utc.year, dt_utc.month, dt_utc.day)
        t1 = ts.utc(dt_utc.year, dt_utc.month, dt_utc.day + 2)
        f = almanac.meridian_transits(eph, moon, observer_loc)
        times, events = almanac.find_discrete(t0, t1, f)
        transit_times = times[events == 1]
        if len(transit_times) > 0:
            # 1. Extract raw naive UTC datetime object from Skyfield
            raw_utc_dt = transit_times[0].utc_datetime().replace(tzinfo=timezone.utc)
            # 2. Convert to system local timezone target
            local_transit_dt = raw_utc_dt.astimezone()
            # 3. Format string to hide seconds (YYYY-MM-DD HH:MM)
            formatted_transit = local_transit_dt.strftime('%Y-%m-%d %H:%M')
        else:
            formatted_transit = "None found in window"

        visible_data.append({"name": "Moon",
                    "common_name": "Moon",
                    "alt": alt_deg.round(4),
                    "az": az_deg.round(4),
                    "size": angular_diameter_arcmin.round(1),
                    "mag": approx_mag.round(1),
                    "type": 3,
                    "exposure": 0.0,
                    "framing": " ",
                    "PA": phase_angle.round(1) ,
                    "moon_sep": 0.0,
                    "priority": "Need",
                    "meridian": ha_hours.round(2),
                    "zenith_date": formatted_transit})

    return visible_data, hdrInfo

# -------------- get tonight's hilites from Telescopius ----------------------------------------'
def getHilites(app):
    client = TelescopiusClient(api_key=Telescopius_api_key)
    dt_this = app.config['CRITERIA']['current_time_utc']
    debug(f"timezone: {str(dt_this.tzinfo)}")

    highlights = []
    status = False
    try:
        highlights = client.get_target_highlights(
            lat=app.config['CRITERIA']['LATITUDE'],
            lon=app.config['CRITERIA']['LONGITUDE'],
            timezone="America/Los_Angeles", # FIXME (?)
            types="DEEP_SKY_OBJECT", # FIXME set to selections (?)
            min_alt=app.config['CRITERIA']['MINALT'],
            #mag_max=app.config['CRITERIA']['visMag'],
            datetime=app.config['CRITERIA']['current_time_utc'].strftime("%Y-%m-%d %H:%M:%S")
            #results_per_page=120
        )
    except TelescopiusAuthError:
        return(status, "Telescopius: Invalid API key")
    except TelescopiusRateLimitError:
        return(status, "Telescopius: Rate limit exceeded - please wait before retrying")
    except TelescopiusBadRequestError as e:
        return(status, f"Telescopius: Invalid parameters: {e.message}")
    except TelescopiusNetworkError:
        return(status, "Telescopius: Network error - check your connection")
    except TelescopiusError as e:
        return(status, f"Telescopius: API error: {e.message}")
    else:
        status = True

    Session = []
    for item in highlights['page_results']:
        obj = item['object']
        # ['main_id', 'main_name', 'ids', 'names', 'family', 'types', 'url', 'main_image_url', 'thumbnail_url', 'alt_ids', 'ra', 'dec', 'con', 'con_name', 'visual_mag']
        obj.get('alt_ids').append(obj.get('main_name') or obj('main_id'))
        #debug(f"Telescopius type: {obj.get('types')}")
        coord = SkyCoord(ra=obj.get('ra'), dec=0, unit=(u.hourangle, u.deg))
        try:
            otype = Telescopius_TO_KSTARS[[item for item in obj.get('types') if item in set(ttypes)][0]]
        except:
            otype = f"**{obj.get('types')}**"

        target = (
            0,
            f"{obj.get('main_name') or obj['main_id']}",
            ', '.join(obj.get('alt_ids')),
            otype,
            round(coord.ra.deg, 4),
            round(obj.get('dec'), 4),
            obj.get('visual_mag'),
            round(obj.get('major_axis')/60, 1),
            0.0,
            obj.get('thumbnail_url'),
            "Need"
            )

        Session.append(target)

    client.close()
    return(status, Session)
