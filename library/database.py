#!/usr/bin/env python3

import os
import sqlite3
import uuid


class Row(dict):
    COLS = []
    def __init__(self, table_row=None):
        if table_row is not None:
            parsed_row = self.parseRow(table_row)
            self.update(parsed_row)
            self.__dict__.update(parsed_row)

    def parseRow(self, r):
        row = {}
        for (col_name, _, caster), r in zip(self.COLS, r):
            casted = caster(r)
            if caster is str and (not casted or casted == "None"):
                casted = None
            row[col_name] = casted

        return row


class CellRow(Row):
    COLS = [
        ("cell_num", "INTEGER PRIMARY KEY", int),
        ("cell_id", "TEXT", str),
        ("experiment_id", "TEXT", str),
        ("start_frame_idx", "INTEGER", int),
        ("end_frame_idx", "INTEGER", int),
        ("birth_observed", "INTEGER DEFAULT 0", bool),
        ("division_observed", "INTEGER DEFAULT 0", bool),
        ("is_wildtype", "INTEGER", bool),
        ("first_outline_id", "TEXT", str),
        ("last_outline_id", "TEXT", str),
        ("parent_cell_id", "TEXT", str),
        ("child_cell_id1", "TEXT", str),
        ("child_cell_id2", "TEXT", str),
    ]


class OutlineRow(Row):
    COLS = [
        ("outline_num", "INTEGER PRIMARY KEY", int),
        ("outline_id", "TEXT", str),
        ("cell_id", "TEXT", str),
        ("experiment_num", "INTEGER", int),
        ("experiment_id", "TEXT", str),
        ("image_path", "TEXT", str),
        ("frame_idx", "INTEGER", int),
        ("coords_path", "TEXT", str),
        ("offset_left", "INTEGER", int),
        ("offset_top", "INTEGER", int),
        ("parent_id", "TEXT", str),
        ("child_id1", "TEXT DEFAULT ''", str),
        ("child_id2", "TEXT DEFAULT ''", str),
    ]


class ExperimentRow(Row):
    COLS = [
        ("experiment_num", "INTEGER PRIMARY KEY", int),
        ("experiment_id", "TEXT", str),
        ("date_year", "INTEGER", int),
        ("date_month", "INTEGER", int),
        ("date_day", "INTEGER", int),
        ("medium", "TEXT", str),
        ("strain", "TEXT", str),
        ("image_path", "TEXT", str),
        ("channel_green", "INTEGER", bool),
        ("channel_red", "INTEGER", bool),
        ("outlined", "INTEGER DEFAULT 0", bool),
        ("verified", "INTEGER DEFAULT 0", bool),
        ("analysed", "INTEGER DEFAULT 0", bool),
    ]
    def parseRow(self, r):
        row = {}
        for (col_name, _, caster), r in zip(self.COLS, r):
            row[col_name] = caster(r)

        return {
            "experiment_num": row["experiment_num"],
            "experiment_id": row["experiment_id"],
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

def createCellsTable():
    query = "CREATE TABLE cells ({0});".format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in CellRow.COLS
    ]))
    executeQuery(query, commit=True)

def insertCell(cell_id, experiment_id, start_frame_idx, end_frame_idx,
               birth_observed, division_observed, is_wildtype,
               first_outline_id, last_outline_id, parent_cell_id,
               child_cell_id1, child_cell_id2):
    query = """
    INSERT INTO cells
    (cell_id, experiment_id, start_frame_idx, end_frame_idx, birth_observed,
    division_observed, is_wildtype, first_outline_id, last_outline_id,
    parent_cell_id, child_cell_id1, child_cell_id2)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    """
    args = (
        cell_id, experiment_id, start_frame_idx, end_frame_idx, birth_observed,
        division_observed, is_wildtype, first_outline_id, last_outline_id,
        parent_cell_id, child_cell_id1, child_cell_id2
    )
    executeQuery(query, args, commit=True)

def getCellById(cell_id):
    query = """
    SELECT *
    FROM cells
    WHERE cell_id = ?;
    """
    args = (cell_id,)
    r = executeQuery(query, args, fetchone=True)
    return CellRow(r)

def deleteCellById(cell_id):
    query = """
    DELETE
    FROM cells
    WHERE cell_id = ?;
    """
    args = (cell_id,)
    executeQuery(query, args, commit=True)

def updateCellById(cell_id, **kwargs):
    args = []
    set_statement = []
    for kw, val in kwargs.items():
        if kw not in [x[0] for x in CellRow.COLS]:
            raise ValueError("Column name {0} is illegal".format(kw))
        set_statement.append("{0} = ?".format(kw))
        args.append(val)

    set_string = ", ".join(set_statement)
    args.append(cell_id)

    query = """
    UPDATE cells
    SET {0}
    WHERE cell_id = ?;
    """.format(set_string)
    executeQuery(query, args, commit=True)

def createOutlinesTable():
    query = "CREATE TABLE outlines ({0});".format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in OutlineRow.COLS
    ]))
    executeQuery(query, commit=True)

def getOutlinesByFrameIdx(frame_idx, experiment_id):
    query = """
    SELECT * FROM outlines
    WHERE experiment_id = ?
      AND frame_idx = ?;
    """
    args = (experiment_id, str(frame_idx))
    results = executeQuery(query, args, fetchmany=True)
    return [OutlineRow(x) for x in results]

def getOutlinesByExperimentId(experiment_id):
    query = """
    SELECT * FROM outlines
    WHERE experiment_id = ?
    """
    args = (experiment_id,)
    results = executeQuery(query, args, fetchmany=True)
    return [OutlineRow(x) for x in results]

def insertOutline(outline_id, cell_id, experiment_num, experiment_id, image_path,
               frame_idx, coords_path, offset_left, offset_top, parent_id):
    query = """
    INSERT INTO outlines
    (outline_id, cell_id, experiment_num, experiment_id, image_path,
     frame_idx, coords_path, offset_left, offset_top, parent_id)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    """
    args = (
        outline_id, cell_id, experiment_num, experiment_id, image_path,
        frame_idx, coords_path, offset_left, offset_top, parent_id,
    )
    new_id = executeQuery(query, args, commit=True)
    return new_id

def updateOutlineById(outline_id, **kwargs):
    args = []
    set_statement = []

    permitted_columns = [
        "cell_id", "parent_id", "child_id1", "child_id2",
    ]
    for kw, val in kwargs.items():
        if kw not in permitted_columns:
            raise ValueError("Column name {0} is illegal".format(kw))
        set_statement.append("{0} = ?".format(kw))
        args.append(val)

    set_string = ", ".join(set_statement)
    args.append(outline_id)

    query = """
    UPDATE outlines
    SET {0}
    WHERE outline_id = ?;
    """.format(set_string)
    executeQuery(query, args, commit=True)

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

def getOutlinesByCellId(cell_id):
    query = """
    SELECT * FROM outlines
    WHERE cell_id = ?;
    """
    args = (cell_id,)
    results = executeQuery(query, args, fetchmany=True)
    return [OutlineRow(x) for x in results]

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
    (experiment_id, date_year, date_month, date_day, medium, strain, image_path, channel_green, channel_red)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);
    """
    experiment_id = str(uuid.uuid4())
    args = (
        experiment_id,
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
    ORDER BY experiment_num;
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
