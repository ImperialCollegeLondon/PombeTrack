#!/usr/bin/env python3

import os
import sqlite3
import uuid


class OutlineRow(dict):
    COLS = [
        ("outline_num", "INTEGER PRIMARY KEY", int),
        ("outline_id", "TEXT", str),
        ("cell_id", "TEXT", str),
        ("experiment_id", "INTEGER", int),
        ("experiment_hash", "TEXT", str),
        ("image_path", "TEXT", str),
        ("frame_idx", "INTEGER", int),
        ("coords_path", "TEXT", str),
        ("offset_left", "INTEGER", int),
        ("offset_top", "INTEGER", int),
        ("parent_id", "TEXT", str),
        ("child_id1", "TEXT DEFAULT ''", str),
        ("child_id2", "TEXT DEFAULT ''", str),
    ]
    def __init__(self, table_row=None):
        if table_row is not None:
            parsed_row = self.parseOutlineRow(table_row)
            self.update(parsed_row)
            self.__dict__.update(parsed_row)

    def parseOutlineRow(self, r):
        row = {}
        for (col_name, _, caster), r in zip(self.COLS, r):
            casted = caster(r)
            if not casted:
                casted = None
            row[col_name] = casted

        return row

class ExperimentRow(dict):
    COLS = [
        ("experiment_id", "INTEGER PRIMARY KEY", int),
        ("experiment_hash", "TEXT", str),
        ("date_year", "INTEGER", int),
        ("date_month", "INTEGER", int),
        ("date_day", "INTEGER", int),
        ("medium", "TEXT", str),
        ("strain", "TEXT", str),
        ("image_path", "TEXT", str),
        ("channel_green", "INTEGER", bool),
        ("channel_red", "INTEGER", bool),
        ("outlined", "INTEGER", bool),
        ("verified", "INTEGER", bool),
        ("analysed", "INTEGER", bool),
    ]
    def __init__(self, table_row=None):
        if table_row is not None:
            parsed_row = self.parseExperimentRow(table_row)
            self.update(parsed_row)
            self.__dict__.update(parsed_row)

    def parseExperimentRow(self, r):
        row = {}
        for (col_name, _, caster), r in zip(self.COLS, r):
            row[col_name] = caster(r)

        return {
            "experiment_id": row["experiment_id"],
            "experiment_hash": row["experiment_hash"],
            "date": "{date_year}-{date_month:02d}-{date_day:02d}".format(**row),
            "medium": row["medium"],
            "strain": row["strain"],
            "image_path": row["image_path"],
            "channels": {"green": bool(row["channel_green"]),
                         "red": bool(row["channel_red"])},
            "outlined": bool(row["outlined"]),
            "verified": bool(row["verified"]),
            "analysed": bool(row["analysed"]),
        }


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
    
def createOutlinesTable():
    query = "CREATE TABLE outlines ({0});".format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in OutlineRow.COLS
    ]))
    executeQuery(query, commit=True)

def getOutlinesByFrameIdx(frame_idx, experiment_hash):
    query = """
    SELECT * FROM outlines
    WHERE experiment_hash = ?
      AND frame_idx = ?
    """
    args = (experiment_hash, frame_idx)
    try:
        results = executeQuery(query, args, fetchmany=True)
    except sqlite3.OperationalError:
        return []
    else:
        return [OutlineRow(x) for x in results]

def getOutlinesByExperimentId(experiment_id):
    query = """
    SELECT * FROM outlines
    WHERE experiment_id = ?
    """
    args = (experiment_id,)
    results = executeQuery(query, args, fetchmany=True)
    return [OutlineRow(x) for x in results]

def insertOutline(outline_id, cell_id, experiment_id, experiment_hash, image_path,
               frame_idx, coords_path, offset_left, offset_top, parent_id):
    query = """
    INSERT INTO outlines
    (outline_id, cell_id, experiment_id, experiment_hash, image_path,
     frame_idx, coords_path, offset_left, offset_top, parent_id)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    """
    args = (
        outline_id, cell_id, experiment_id, experiment_hash, image_path,
        frame_idx, coords_path, offset_left, offset_top, parent_id,
    )
    new_id = executeQuery(query, args, commit=True)
    return new_id

def addOutlineChild(outline_id, child1, child2=None):
    if child2:
        query = """
        UPDATE outlines
        SET child_id1 = ?, child_id2 = ?
        WHERE outline_id = ?;
        """
        args = (child1, child2, outline_id)
    else:
        query = """
        UPDATE outlines
        SET child_id1 = ?
        WHERE outline_id = ?;
        """
        args = (child1, outline_id)
    
    executeQuery(query, args, commit=True)

def deleteOutlineById(outline_id):
    query = """
    DELETE FROM outlines
    WHERE outline_id = ?;
    """
    args = (outline_id,)
    executeQuery(query, args, commit=True)

    query = """
    UPDATE outlines
    SET parent_id = ''
    WHERE parent_id = ?;
    """
    executeQuery(query, args, commit=True)

    query = """
    UPDATE outlines
    SET child_id1 = ''
    WHERE child_id1 = ?;
    """
    executeQuery(query, args, commit=True)

    query = """
    UPDATE outlines
    SET child_id2 = ''
    WHERE child_id2 = ?;
    """
    executeQuery(query, args, commit=True)

def getOutlineById(outline_id):
    query = """
    SELECT * FROM outlines
    WHERE outline_id = ?;
    """
    args = (outline_id,)
    result = executeQuery(query, args, fetchone=True)
    return OutlineRow(result)

def createExperimentsTable():
    query = "CREATE TABLE experiments ({0});".format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in ExperimentRow.COLS
    ]))
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

def insertExperiment(date_year, date_month, date_day, medium, strain,
                     image_path, channel_green, channel_red):
    query = """
    INSERT INTO experiments
    (experiment_hash, date_year, date_month, date_day, medium, strain, image_path, channel_green, channel_red)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);
    """
    experiment_hash = str(uuid.uuid4())
    args = (
        experiment_hash,
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
        experiments = [ExperimentRow(x)
                       for x in results]
        return experiments

def getExperimentById(experiment_id):
    query = """
    SELECT *
    FROM experiments
    WHERE experiment_id = ?
    """
    args = (experiment_id,)
    result = executeQuery(query, args, fetchone=True)
    return ExperimentRow(result)

def deleteExperimentById(experiment_id):
    query = """
    DELETE FROM experiments
    WHERE experiment_id = ?
    """
    args = (experiment_id,)
    executeQuery(query, args, commit=True)

if __name__ == "__main__":
    print("This should be imported")
    print("In the future, I might add some export functionality here or something")
