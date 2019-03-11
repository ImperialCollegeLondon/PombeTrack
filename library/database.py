#!/usr/bin/env python3

import os
import sqlite3


def executeQuery(query, args=None, commit=False, fetchone=False, fetchmany=False):
    db_path = os.path.join(".analysis_store", "pombetrack.db")
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        if args:
            cursor.execute(query, args)
        else:
            cursor.execute(query)

        if commit:
            conn.commit()
       
        if fetchone:
            result = cursor.fetchone()
        elif fetchmany:
            result = cursor.fetchall()
        else:
            result = cursor.lastrowid

    return result

def checkTable(table_name):
    query = """
    SELECT name
    FROM sqlite_master
    WHERE type='table'
    AND name=?;
    """
    args = (table_name,)
    return executeQuery(query, args, fetchone=True)
    
def createExperimentsTable():
    query = """
    CREATE TABLE experiments
    (experiment_id   INTEGER PRIMARY KEY,
     date_year       INTEGER,
     date_month      INTEGER,
     date_day        INTEGER,
     medium          TEXT,
     strain          TEXT,
     image_path      TEXT,
     channel_green   INTEGER,
     channel_red     INTEGER);
    """
    executeQuery(query, commit=True)

def checkExperimentDuplicate(date_year, date_month, date_day, medium, strain):
    query = """
    SELECT experiment_id
    FROM experiments
    WHERE date_year = ?
    AND date_month = ?
    AND date_day = ?
    AND medium = ?
    AND strain = ?;
    """
    args = (date_year, date_month, date_day, medium, strain)
    experiments = executeQuery(query, args, fetchone=True)
    if experiments:
        return experiments[0]
    else:
        return False

def insertExperiment(date_year, date_month, date_day, medium, strain, image_path, channel_green, channel_red):
    query = """
    INSERT INTO experiments
    (date_year, date_month, date_day, medium, strain, image_path, channel_green, channel_red)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?);
    """
    args = (
        date_year, date_month, date_day,
        medium, strain, image_path,
        channel_green, channel_red,
    )
    new_id = executeQuery(query, args, commit=True)
    return new_id

def updateExperimentById(experiment_id, image_path, channel_green, channel_red):
    query = """
    UPDATE experiments
    SET image_path = ?,
        channel_green = ?,
        channel_red = ?
    WHERE experiment_id = ?;
    """
    args = (image_path, channel_green, channel_red, experiment_id)
    executeQuery(query, args, commit=True)

def getExperiments():
    query = """
    SELECT *
    FROM experiments;
    """
    results = executeQuery(query, fetchmany=True)
    experiments = [
        dict(zip(
            ["experiment_id", "date", "medium", "strain", "image_path", "channels"],
            [x[0], "{0}-{1:02d}-{2:02d}".format(x[1], x[2], x[3]),
             x[4], x[5], x[6], {"green": bool(x[7]), "red": bool(x[8])}]
        ))
        for x in results 
    ]
    return experiments


if __name__ == "__main__":
    print("This should be imported")
    print("In the future, I might add some export functionality here or something")
