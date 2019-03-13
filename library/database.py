#!/usr/bin/env python3

import os
import sqlite3


def executeQuery(query, args=None, commit=False, fetchone=False, fetchmany=False):
    db_path = os.path.join("data", "pombetrack.db")
    if not os.path.exists(os.path.dirname(db_path)):
        os.makedirs(os.path.dirname(db_path))

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
     channel_red     INTEGER,
     outlined        INTEGER DEFAULT 0,
     verified        INTEGER DEFAULT 0,
     analysed        INTEGER DEFAULT 0);
    """
    executeQuery(query, commit=True)

def checkExperimentDuplicate(date_year, date_month, date_day, medium, strain, image_path):
    query = """
    SELECT experiment_id
    FROM experiments
    WHERE date_year = ?
    AND date_month = ?
    AND date_day = ?
    AND medium = ?
    AND strain = ?
    AND image_path = ?;
    """
    args = (date_year, date_month, date_day, medium, strain, image_path)
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

def updateExperimentById(experiment_id, **kwargs):
    args = []
    set_statement = []

    permitted_columns = [
        "medium", "strain", "image_path", "channel_green", "channel_red",
        "outlined", "verified", "analysed",
    ]
    for kw, val in kwargs.items():
        if kw not in permitted_columns:
            raise ValueError("Column name {0} is illegal".format(kw))
        set_statement.append("{0} = ?".format(kw))
        args.append(val)

    set_string = ", ".join(set_statement)
    args.append(experiment_id)

    query = """
    UPDATE experiments
    SET {0}
    WHERE experiment_id = ?;
    """.format(set_string)
    executeQuery(query, args, commit=True)

def getExperiments():
    query = """
    SELECT *
    FROM experiments
    ORDER BY experiment_id;
    """
    try:
        results = executeQuery(query, fetchmany=True)
    except sqlite3.OperationalError:
        return []
    else:
        experiments = [_parseExperimentRow(x)
                       for x in results]
        return experiments

def _parseExperimentRow(row):
    return {
        "experiment_id": row[0],
        "date": "{0}-{1:02d}-{2:02d}".format(row[1], row[2], row[3]),
        "medium": row[4],
        "strain": row[5],
        "image_path": row[6],
        "channels": {"green": bool(row[7]),
                     "red": bool(row[8])},
        "outlined": bool(row[9]),
        "verified": bool(row[10]),
        "analysed": bool(row[11]),
    }


if __name__ == "__main__":
    print("This should be imported")
    print("In the future, I might add some export functionality here or something")
