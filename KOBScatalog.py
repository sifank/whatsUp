#!/usr/bin/python3
# Program: KOBScatalog.py
# Version: 20260906
# Author:  Sifan Kahale
# Desc:    manages the KOBS object database

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

# Search type for manageDB
SEARCHTYPE = {
    "Name": "Name",
    "Type": "type",
    "Magnitude": "magnitude",
    "Size": "apSize",
    "Image Loc": "imageLoc",
    "Priority": "procFlag"
    }

# priority types
PROCFLAG = [
    "Need",
    "Done",
    "Redo",
    "Should Redo"
    ]

#ID = 0

app = Flask(__name__)
app.secret_key = 'KOBS Catalog'

#================Manage DataBase ==============================
@app.route('/', methods=['GET', 'POST'])
def index():
    conn, cursor = openDB()
    if conn == "None":
        return render_template('SQLerror.html', err="Can not open Database")

    sql = f"SELECT * FROM KBcatalog limit 50;"
    try:
        #debug(f"/ sql: {sql}")
        cursor.execute(sql)
        Session = cursor.fetchall()
    except mysql.connector.Error as err:
        closeDB(conn, cursor)
        return render_template('SQLerror.html', err=err)
    closeDB(conn, cursor)

    '''
    # add best month, dt, alt to end of each row of session
    for row in Session:
        debug(f"type of row: {type(row)}, row: {row}")
        month, medianTime, targetAlt = medianXfer(row[1])
        debug(f"Target: {row[1]}, {month}, {medianTime}, {targetAlt}")
    '''

    return render_template('whatsUpDBmain.html', Session=Session, DSO=DSOTYPE, SEARCHTYPE= SEARCHTYPE, PF = PROCFLAG, SQL = sql)

#--------------------------------------------------------
@app.route('/updatePriority', methods=['POST'])
def updatePriority():
    Priority = request.form.get('Priority')
    row = request.form.get('Row')
    SQL = request.form.get('SQL')
    conn, cursor = openDB()
    if conn == "None":
        return render_template('SQLerror.html', err="Can not open Database")

    #---- update prority on this target ------
    sql = f"UPDATE KBcatalog SET procFlag = '{Priority}' where ID = {row};"
    try:
        #debug(f"update pri - update sql: {sql}")
        cursor.execute(sql)
        conn.commit()
    except mysql.connector.Error as err:
        closeDB(conn, cursor)
        debug(f"sql error: {err}")
        return render_template('SQLerror.html', err=err)

    # ---- reread results list using preivous search info -----
    try:
        #debug(f"update pri - reread sql: {SQL}")
        cursor.execute(SQL)
        Session = cursor.fetchall()
    except mysql.connector.Error as err:
        closeDB(conn, cursor)
        debug(f"sql error: {err}")
        return render_template('SQLerror.html', err=err)

    closeDB(conn, cursor)

    return render_template('whatsUpDBmain.html', Session=Session, DSO=DSOTYPE, SEARCHTYPE= SEARCHTYPE, SQL = SQL)

#--------------------------------------------------------
@app.route('/search', methods=['POST'])
def search():
    sType = request.form.get('stype')
    sRegex = request.form.get('sregex').replace("'", "\\'")
    #debug(f"You entered a search for {sType} looking for: {sRegex}")

    if sType == "Name":
        sql = f"SELECT * FROM KBcatalog WHERE objectName RLIKE \'{sRegex}\' or commonName RLIKE \'{sRegex}\';"
    elif sType == "Type":
        sql = f"SELECT * FROM KBcatalog WHERE type RLIKE \'{sRegex}\';"
    else:
        sql = f"SELECT * FROM KBcatalog WHERE {sType} RLIKE \'{sRegex}\';"

    conn, cursor = openDB()
    try:
        #debug(f"Search SQL: {sql}")
        cursor.execute(sql)
        Session = cursor.fetchall()
    except mysql.connector.Error as err:
        closeDB(conn, cursor)
        return render_template('SQLerror.html', err=err)

    if cursor.rowcount == 0:
        return render_template('SQLerror.html', err = f"'{sType}' not found in KBcatalog")
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
        return redirect('//')

    sql = f"INSERT INTO KBcatalog SET \
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
    conn, cursor = openDB()

    try:
        #debug(f"add sql: {sql}")
        cursor.execute(sql)
        conn.commit()
    except mysql.connector.Error as err:
        closeDB(conn, cursor)
        return render_template('SQLerror.html', err=err)
    closeDB(conn, cursor)

    return redirect('//')

#---------------- fill in add/modify form with selected row ----------------------------------------
@app.route('/fillWithSelected', methods=['POST'])
def fillWithSelected():
    debug(f"fillWithSelected returned Sel: {request.form.get('SEL')}")
    SEL = list(ast.literal_eval(request.form.get('SEL')))
    debug(f"SEL as list: {SEL}")
    Session = session.get('onlineSearch', [])

    return render_template('whatsUpSearchDB.html', Session=Session, SEL = SEL, DSO=DSOTYPE, SEARCHTYPE= SEARCHTYPE, PF = PROCFLAG)

#--------------------------------------------------------
# TODO if imgLoc is url, download, place in spc folder and rewrite imgLoc to point to disk
@app.route('/objSearch', methods=['POST'])
def objSearch():
    newObject = request.form.get('lookupObject').replace("'", "\\'")
    OpType = request.form.get('OpType')
    #debug(f"OpType in objSearch: {OpType}")
    if not newObject or newObject.strip() == '':
        return render_template('SQLerror.html', err=f"**{newObject}** not a valid search term")
    Session = []
    #debug(f"objSearch - ##{newObject}##")
    # ---- lookup in databases -------------------------------------------------------------------
    sql = f"SELECT 'KBcatalog' AS source_table, t1.* FROM KBcatalog t1 WHERE t1.objectname LIKE '%{newObject}%'  or t1.commonName = '%{newObject}%' \
        UNION ALL SELECT 'IPcatalog' AS source_table, t2.* FROM IPcatalog t2 WHERE t2.objectname LIKE '%{newObject}%'  or t2.commonName = '%{newObject}%' \
        UNION ALL SELECT 'DScatalog' AS source_table, t3.* FROM DScatalog t3 WHERE t3.objectname LIKE '%{newObject}%'  or t3.commonName = '%{newObject}%';"
    conn, cursor = openDB()
    try:
        #debug(f"search find in db sql: {sql}")
        cursor.execute(sql)
        Session = cursor.fetchall()
    except mysql.connector.Error as err:
        closeDB(conn, cursor)
        #debug(f"objSearch sql error: {err}")
        return render_template('SQLerror.html', err=err)

    closeDB(conn, cursor)
    #debug(f"Database results: {inSession}")
    if not Session:
        Session = []

    i = 1 #  set id to numbers starting with 1 (db starts above 300)
    # ---- lookup in Simbad ----------------------------------------------------
    simObj = lookup_objectSimbad(newObject)
    if simObj[0] == "Error":
        debug(f"Simbad error: {simObj[1]}")
        #return render_template('displayMsg.html', err=f"Simbad can not find: {simObj[1]}")
    else:
        simObj[1] = i
        #debug(f"(whatsUp) Simbad results: {simObj}")
        i += 1
        Session.append(simObj)    # only add if found

    # ---- lookup in Telescopius ----------------------------------------------------------
    terror, numObjects, foundObj = lookup_objectTelescopius(newObject)
    if terror != "Success":
        debug(f"Telescopius error: {terror}")
        #return render_template('displayMsg.html', err=f"Telescopius error: {terror}")
    elif numObjects == 0:
        debug(f"Telescopius error: no objects found")
        #return render_template('displayMsg.html', err=f"Telescopius: none found")
    else:
        for obj in foundObj:
            obj = obj[:1] + (i,) + obj[2:]  #setting id to be unique
            i += 1
            Session.append(obj)

    # ---- test if KBcatalog has this entry and fill edit form with it, if not, set to blank form
    #debug(f"objSearch: Session is: {Session}")
    SEL = next((item for item in Session if item[0] == 'KBcatalog'), None)
    if not SEL:
        SEL = ["", 0, "", "", 99, 0.0, 0.0, 0.0, 0.0, 0.0, "", "Need"]

    session['onlineSearch'] = Session
    #debug(f"svdSession: {session.get('onlineSearch', [])}")
    return render_template('whatsUpSearchDB.html', Session=Session, SEL = SEL, DSO=DSOTYPE, SEARCHTYPE= SEARCHTYPE, PF = PROCFLAG, OT = OpType)

#--------------------------------------------------------
@app.route('/modify', methods=['POST'])
def modify():
    ID = request.form.get('ID')
    sql = f"SELECT * FROM KBcatalog WHERE ID={ID};"
    conn, cursor = openDB()
    try:
        #debug(f"Modify SQL: {sql}")
        cursor.execute(sql)
        Session = cursor.fetchone()
    except mysql.connector.Error as err:
        closeDB(conn, cursor)
        return render_template('SQLerror.html', err=err)
    closeDB(conn, cursor)

    if cursor.rowcount == 0:
        return render_template('SQLerror.html', err = f"{ID} in KBcatalog not found")

    #debug(f"\nModify for ID: {ID} OTA: {Session[9]} Notes: {Session[8]} Objective: {Session[7]}\n")
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

    #debug(f"ID: {ID}: {objectName} {commonName} {Type} {ra} {dec} {mag} {Size} {PA} {imageLoc} {Priority}")

    sql = f"INSERT INTO KBcatalog \
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

    conn, cursor = openDB()
    try:
        #debug(f"ModifyOK SQL {sql}")
        cursor.execute(sql)
        conn.commit()
    except mysql.connector.Error as err:
        closeDB(conn, cursor)
        return render_template('SQLerror.html', err=err)
    closeDB(conn, cursor)

    return redirect('/')

#--------------------------------------------------------
@app.route('/removeok', methods=['POST'])
def removeok():
    ID = request.form.get('ID')
    sql = f"DELETE FROM KBcatalog WHERE ID={ID};"
    conn, cursor = openDB()
    try:
        #debug(f"removeok SQL {sql}")
        cursor.execute(sql)
        conn.commit()
    except mysql.connector.Error as err:
        closeDB(conn, cursor)
        return render_template('SQLerror.html', err=err)
    closeDB(conn, cursor)

    return redirect('/')

#============================================================
# ---- fav lists db manage -----------------------------------

#--------------------------------------------------------
if __name__ == '__main__':
    app.run(host='0.0.0.0',port=5011, debug=True)

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
  `bestMonth` varchar(10) DEFAULT NULL,
  `bestDT` varchar(30) DEFAULT NULL,
  `DTalt` float DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `idx_objectName` (`objectName`)
) ENGINE=InnoDB
'''
