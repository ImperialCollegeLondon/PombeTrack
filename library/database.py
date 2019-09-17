#!/usr/bin/env python3

import os
import pathlib
import shutil
import sqlite3
import sys
import time
import uuid

from . import loader


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


class VersionRow(Row):
    COLS = [
        ("major_version", "INTEGER", int),
        ("minor_version", "INTEGER", int),
    ]


class AssociationRow(Row):
    COLS = [
        ("association_num", "INTEGER PRIMARY KEY", int),
        ("association_id", "TEXT", str),
        ("experiment_id", "TEXT", str),
        ("associated_experiment_id", "TEXT", str),
        ("association_type", "TEXT", str),
    ]


class NucleusRow(Row):
    COLS = [
        ("nucleus_num", "INTEGER PRIMARY KEY", int),
        ("nucleus_id", "TEXT", str),
        ("outline_id", "TEXT", str),
        ("cell_id", "TEXT", str),
        ("experiment_id", "TEXT", str),
        ("coords_path", "TEXT", str),
    ]
    def parseRow(self, r):
        row = super().parseRow(r)
        if "\\" in row["coords_path"]:
            row["coords_path"] = pathlib.PureWindowsPath(row["coords_path"]).as_posix()

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
    def parseRow(self, r):
        row = super().parseRow(r)
        if "\\" in row["coords_path"]:
            row["coords_path"] = pathlib.PureWindowsPath(row["coords_path"]).as_posix()

        return row


class SignalRow(Row):
    BACKGROUND = 1
    COLS = [
        ("signal_num", "INTEGER PRIMARY KEY", int),
        ("signal_id", "TEXT", str),
        ("outline_id", "TEXT", str),
        ("cell_id", "TEXT", str),
        ("experiment_id", "TEXT", str),
        ("channel", "INTEGER", int),
        ("signal_type", "INTEGER DEFAULT 0", int),
        ("signal_value", "REAL", float),
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
        ("image_mode", "TEXT", str),
        ("num_channels", "INTEGER", int),
        ("num_slices", "INTEGER", int),
        ("num_frames", "INTEGER", int),
        ("file_mode", "TEXT", str),
        ("outlined", "INTEGER DEFAULT 0", bool),
        ("verified", "INTEGER DEFAULT 0", bool),
        ("analysed", "INTEGER DEFAULT 0", bool),
    ]
    def parseRow(self, r):
        row = {}
        for (col_name, _, caster), r in zip(self.COLS, r):
            if r is None:
                row[col_name] = None
            else:
                row[col_name] = caster(r)

        if "\\" in row["image_path"]:
            row["image_path"] = pathlib.PureWindowsPath(row["image_path"]).as_posix()

        return {
            "experiment_num": row["experiment_num"],
            "experiment_id": row["experiment_id"],
            "date": "{date_year}-{date_month:02d}-{date_day:02d}".format(**row),
            "medium": row["medium"],
            "strain": row["strain"],
            "image_path": row["image_path"],
            "image_mode": row["image_mode"] == 1 and "movie" or "static",
            "num_channels": row["num_channels"],
            "num_slices": row["num_slices"],
            "num_frames": row["num_frames"],
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

def createVersionTable():
    query = """
    CREATE TABLE version
    (major_version INTEGER, minor_version INTEGER);
    """
    query = """
    CREATE TABLE version
    ({0});
    """.format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in VersionRow.COLS
    ]))
    executeQuery(query, commit=True)

def insertVersion(major, minor):
    query = """
    INSERT INTO version
    (major_version, minor_version)
    VALUES
    (?, ?);
    """
    args = (major, minor)
    executeQuery(query, args, commit=True)

def updateVersion(major, minor):
    query = """
    UPDATE version
    SET major_version = ?,
        minor_version = ?;
    """
    args = (major, minor)
    executeQuery(query, args, commit=True)

def getVersion():
    query = """
    SELECT *
    From version
    """
    result = executeQuery(query, fetchone=True)
    version = VersionRow(result)
    return (version["major_version"], version["minor_version"])

def createSignalsTable():
    query = """
    CREATE TABLE signals
    ({0});
    """.format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in SignalRow.COLS
    ]))
    executeQuery(query, commit=True)

def insertSignal(signal_id, outline_id, cell_id, experiment_id, channel, signal_type, signal_value):
    query = """
    INSERT INTO signals
    (signal_id, outline_id, cell_id, experiment_id, channel, signal_type, signal_value)
    VALUES (?, ?, ?, ?, ?, ?, ?);
    """
    args = (signal_id, outline_id, cell_id, experiment_id, channel, signal_type, signal_value)
    executeQuery(query, args, commit=True)

def getSignalById(signal_id):
    query = """
    SELECT *
    FROM signals
    WHERE signal_id = ?;
    """
    args = (signal_id,)
    r = executeQuery(query, args, fetchone=True)
    return SignalRow(r)

def getSignalsByCriteria(**kwargs):
    selectors = []
    args = []
    for k, v in kwargs.items():
        selectors.append("{0} = ?".format(k))
        args.append(v)

    query = """
    SELECT *
    FROM signals
    WHERE {0};
    """.format("\nAND ".join(selectors))
    r = executeQuery(query, args, fetchmany=True)
    return [SignalRow(x) for x in r]

def deleteSignalsByCriteria(**kwargs):
    selectors = []
    args = []
    for k, v in kwargs.items():
        selectors.append("{0} = ?".format(k))
        args.append(v)

    query = """
    DELETE
    FROM signals
    WHERE {0};
    """.format("\nAND ".join(selectors))
    executeQuery(query, args, commit=True)

def createAssociationsTable():
    query = "CREATE TABLE associations ({0});".format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in AssociationRow.COLS
    ]))
    executeQuery(query, commit=True)

def insertAssociation(association_id, experiment_id, associated_id, association_type):
    query = """
    INSERT INTO associations
    (association_id, experiment_id, associated_experiment_id, association_type)
    VALUES (?, ?, ?, ?);
    """
    args = (association_id, experiment_id, associated_id, association_type)
    executeQuery(query, args, commit=True)

def deleteAssociationById(association_id):
    query = """
    DELETE FROM associations
    WHERE association_id = ?;
    """
    args = (association_id,)
    executeQuery(query, args, commit=True)

def getAssociationsByExperimentId(experiment_id, association_type):
    query = """
    SELECT *
    FROM associations
    WHERE experiment_id = ?
      AND association_type = ?;
    """
    args = (experiment_id, association_type)
    r = executeQuery(query, args, fetchmany=True)
    return [AssociationRow(x) for x in r]

def updateAssociationById(association_id, **kwargs):
    args = []
    set_statement = []
    for kw, val in kwargs.items():
        if kw not in [x[0] for x in AssociationRow.COLS]:
            raise ValueError("Column name {0} is illegal".format(kw))
        set_statement.append("{0} = ?".format(kw))
        args.append(val)

    set_string = ", ".join(set_statement)
    args.append(association_id)

    query = """
    UPDATE cells
    SET {0}
    WHERE association_id = ?;
    """.format(set_string)
    executeQuery(query, args, commit=True)

def createNucleiTable():
    query = "CREATE TABLE nuclei ({0});".format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in NucleusRow.COLS
    ]))
    executeQuery(query, commit=True)

def insertNucleus(nucleus_id, outline_id, cell_id, experiment_id, coords_path):
    query = """
    INSERT INTO nuclei
    (nucleus_id, outline_id, cell_id, experiment_id, coords_path)
    VALUES (?, ?, ?, ?, ?);
    """
    args = (nucleus_id, outline_id, cell_id, experiment_id, coords_path)
    executeQuery(query, args, commit=True)

def deleteNucleusById(nucleus_id):
    query = """
    DELETE
    FROM nuclei
    WHERE nucleus_id = ?;
    """
    args = (nucleus_id,)
    executeQuery(query, args, commit=True)

def getNucleusById(nucleus_id):
    query = """
    SELECT *
    FROM nuclei
    WHERE nucleus_id = ?;
    """
    args = (nucleus_id,)
    r = executeQuery(query, args, fetchone=True)
    return NucleusRow(r)

def getNucleiByExperimentId(experiment_id):
    query = """
    SELECT *
    FROM nuclei
    WHERE experiment_id = ?;
    """
    args = (experiment_id,)
    results = executeQuery(query, args, fetchmany=True)
    return [NucleusRow(x) for x in results]

def getNucleiByCellId(cell_id):
    query = """
    SELECT *
    FROM nuclei
    WHERE cell_id = ?;
    """
    args = (cell_id,)
    results = executeQuery(query, args, fetchmany=True)
    return [NucleusRow(x) for x in results]

def getNucleiByOutlineId(outline_id):
    query = """
    SELECT *
    FROM nuclei
    WHERE outline_id = ?;
    """
    args = (outline_id,)
    results = executeQuery(query, args, fetchmany=True)
    return [NucleusRow(x) for x in results]

def deleteNucleiByExperimentId(experiment_id):
    query = """
    DELETE
    FROM nuclei
    WHERE experiment_id = ?;
    """
    args = (experiment_id,)
    executeQuery(query, args, commit=True)

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

def getCellsByExperimentId(experiment_id, birth_observed=None, division_observed=None, is_wildtype=None):
    added_str = ""
    args = [experiment_id,]
    if birth_observed is not None:
        added_str += " AND birth_observed = ?"
        args.append(birth_observed)

    if division_observed is not None:
        added_str += " AND division_observed = ?"
        args.append(division_observed)

    if is_wildtype is not None:
        added_str += " AND is_wildtype = ?"
        args.append(is_wildtype)

    query = """
    SELECT *
    FROM cells
    WHERE experiment_id = ?{0};
    """.format(added_str)
    r = executeQuery(query, args, fetchmany=True)
    return [CellRow(x) for x in r]

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
                     image_path, image_mode, num_channels, num_slices,
                     num_frames, file_mode):
    query = """
    INSERT INTO experiments
    (experiment_id, date_year, date_month, date_day, medium, strain,
     image_path, image_mode, num_channels, num_slices, num_frames, file_mode)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    """
    experiment_id = str(uuid.uuid4())
    args = (
        experiment_id,
        date_year, date_month, date_day, medium, strain, image_path,
        image_mode, num_channels, num_slices, num_frames, file_mode,
    )
    new_id = executeQuery(query, args, commit=True)
    return new_id, experiment_id

def updateExperimentById(experiment_id, **kwargs):
    args = []
    set_statement = []

    permitted_columns = [
        "medium", "strain", "image_path", "image_mode", "num_channels",
        "num_slices", "num_frames", "file_mode", "outlined", "verified", "analysed",
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


def backup_tables():
    backup_str = "{0}-pombetrack.db.backup".format(time.strftime("%Y-%m-%d"))
    backup_path = os.path.join("data", "backups", backup_str)
    backup_num = 1
    while os.path.exists(backup_path):
        backup_path = os.path.join("data", "backups", "{0}-{1}".format(
            backup_str, backup_num
        ))
        backup_num += 1

    db_path = os.path.join("data", "pombetrack.db")
    if not os.path.exists(os.path.dirname(backup_path)):
        os.makedirs(os.path.dirname(backup_path))

    shutil.copyfile(db_path, backup_path)
    print("Backed up database to {0}".format(backup_path))


def run_database_updates(from_version, to_version):
    update_sequence = [
        ((0, 0), (0, 1), _update1),
        ((0, 1), (0, 2), _update2),
        ((0, 2), (0, 3), _update3),
    ]
    for seq_prev, seq_next, update_func in update_sequence:
        if seq_prev == from_version and seq_prev != to_version:
            print("Updating from version {0} to {1}".format(seq_prev, seq_next))
            backup_tables()
            update_func()
            from_version = seq_next

def _update1():
    print("Adding image_mode, num_frames, num_channels, num_slices columns to "
          "experiments table and removing channels column")
    old_columns = ",".join([
        "experiment_num", "experiment_id", "date_year", "date_month",
        "date_day", "medium", "strain", "image_path", "outlined", "verified",
        "analysed",
    ])
    new_cols = [
        ("experiment_num", "INTEGER PRIMARY KEY", int),
        ("experiment_id", "TEXT", str),
        ("date_year", "INTEGER", int),
        ("date_month", "INTEGER", int),
        ("date_day", "INTEGER", int),
        ("medium", "TEXT", str),
        ("strain", "TEXT", str),
        ("image_path", "TEXT", str),
        ("image_mode", "INTEGER DEFAULT 1", int),
        ("num_channels", "INTEGER", int),
        ("num_slices", "INTEGER", int),
        ("num_frames", "INTEGER", int),
        ("outlined", "INTEGER DEFAULT 0", bool),
        ("verified", "INTEGER DEFAULT 0", bool),
        ("analysed", "INTEGER DEFAULT 0", bool),
    ]

    create_query = "CREATE TABLE experiments ({0});".format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in new_cols
    ]))
    queries = [
        "CREATE TABLE _backup({0});".format(old_columns),
        """
        INSERT INTO _backup
        SELECT {0} FROM experiments;
        """.format(old_columns),
        "DROP TABLE experiments;",
        create_query,
        "INSERT INTO experiments ({0}) SELECT {0} FROM _backup;".format(old_columns),
        "DROP TABLE _backup;",
    ]
    db_path = os.path.join("data", "pombetrack.db")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    for query in queries:
        cursor.execute(query)
        conn.commit()

    print("Inspecting image files for channel information")
    cursor.execute("SELECT * FROM experiments;")
    experiments = cursor.fetchall()

    for experiment in experiments:
        # load image for movie/static, num_frames, num_channels, num_slices
        image_path = experiment[7]
        if not os.path.exists(image_path):
            print("! File {0} is not accessible, skipping it".format(
                image_path
            ))
            num_frames, num_channels, num_slices = 1, 1, 1
        else:
            meta = loader.ImageLoader(image_path).im_metadata
            if "frames" in meta:
                num_frames = meta["frames"]
            else:
                num_frames = 1

            if "channels" in meta:
                num_channels = meta["channels"]
            else:
                num_channels = 1

            if "slices" in meta:
                num_slices = meta["slices"]
            else:
                num_slices = 1

            if num_frames == 1:
                image_mode = 2
            else:
                image_mode = 1

        query = """
        UPDATE experiments
        SET image_mode = ?,
            num_frames = ?,
            num_channels = ?,
            num_slices = ?
        WHERE experiment_id = ?;
        """
        args = (image_mode, num_frames, num_channels, num_slices, experiment[1])
        cursor.execute(query, args)
        conn.commit()

    query = "UPDATE version SET major_version = ?, minor_version = ?;"
    args = (0, 1)
    cursor.execute(query, args)
    conn.commit()
    conn.close()

def _update2():
    print("Changing image_mode column to TEXT")
    new_cols = [
        ("experiment_num", "INTEGER PRIMARY KEY", int),
        ("experiment_id", "TEXT", str),
        ("date_year", "INTEGER", int),
        ("date_month", "INTEGER", int),
        ("date_day", "INTEGER", int),
        ("medium", "TEXT", str),
        ("strain", "TEXT", str),
        ("image_path", "TEXT", str),
        ("image_mode", "TEXT", str),
        ("num_channels", "INTEGER", int),
        ("num_slices", "INTEGER", int),
        ("num_frames", "INTEGER", int),
        ("outlined", "INTEGER DEFAULT 0", bool),
        ("verified", "INTEGER DEFAULT 0", bool),
        ("analysed", "INTEGER DEFAULT 0", bool),
    ]
    col_names = [x[0] for x in new_cols]
    col_subset = [x[0] for x in new_cols if x[0] != "image_mode"]
    db_path = os.path.join("data", "pombetrack.db")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    query = "CREATE TABLE _backup({0});".format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in new_cols
    ]))
    cursor.execute(query)
    conn.commit()

    query = "INSERT INTO _backup ({0}) SELECT {0} FROM experiments;".format(
        ",".join(col_subset)
    )
    cursor.execute(query)
    conn.commit()

    query = "SELECT experiment_id, image_mode FROM experiments;"
    cursor.execute(query)
    for experiment in cursor.fetchall():
        query = "UPDATE _backup SET image_mode = ? WHERE experiment_id = ?"
        if experiment[1] == 1:
            args = ("movie", experiment[0])
        if experiment[1] == 2:
            args = ("static", experiment[0])

        cursor.execute(query, args)
    conn.commit()

    query = "DROP TABLE experiments;"
    cursor.execute(query)
    conn.commit()

    create_query = "CREATE TABLE experiments ({0});".format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in new_cols
    ]))
    cursor.execute(create_query)
    conn.commit()

    query = "INSERT INTO experiments ({0}) SELECT {0} from _backup;".format(
        ",".join(col_names),
    )
    cursor.execute(query)

    query = "DROP TABLE _backup;"
    cursor.execute(query)
    conn.commit()

    query = "UPDATE version SET major_version = ?, minor_version = ?;"
    args = (0, 2)
    cursor.execute(query, args)
    conn.commit()
    conn.close()

def _update3():
    print("Adding file_mode column to experiments table")
    new_cols = [
        ("experiment_num", "INTEGER PRIMARY KEY", int),
        ("experiment_id", "TEXT", str),
        ("date_year", "INTEGER", int),
        ("date_month", "INTEGER", int),
        ("date_day", "INTEGER", int),
        ("medium", "TEXT", str),
        ("strain", "TEXT", str),
        ("image_path", "TEXT", str),
        ("image_mode", "TEXT", str),
        ("num_channels", "INTEGER", int),
        ("num_slices", "INTEGER", int),
        ("num_frames", "INTEGER", int),
        ("file_mode", "TEXT", str),
        ("outlined", "INTEGER DEFAULT 0", bool),
        ("verified", "INTEGER DEFAULT 0", bool),
        ("analysed", "INTEGER DEFAULT 0", bool),
    ]
    col_names = [x[0] for x in new_cols]
    col_subset = [x[0] for x in new_cols if x[0] != "file_mode"]
    db_path = os.path.join("data", "pombetrack.db")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    query = "CREATE TABLE _backup({0});".format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in new_cols
    ]))
    cursor.execute(query)
    conn.commit()

    query = "INSERT INTO _backup ({0}) SELECT {0} FROM experiments;".format(
        ",".join(col_subset)
    )
    cursor.execute(query)
    conn.commit()

    query = "UPDATE _backup SET file_mode = 'single';"
    cursor.execute(query)
    conn.commit()

    query = "DROP TABLE experiments;"
    cursor.execute(query)
    conn.commit()

    create_query = "CREATE TABLE experiments ({0});".format(",".join([
        "{0} {1}".format(x[0], x[1])
        for x in new_cols
    ]))
    cursor.execute(create_query)
    conn.commit()

    query = "INSERT INTO experiments ({0}) SELECT {0} from _backup;".format(
        ",".join(col_names),
    )
    cursor.execute(query)

    query = "DROP TABLE _backup;"
    cursor.execute(query)
    conn.commit()

    query = "UPDATE version SET major_version = ?, minor_version = ?;"
    args = (0, 3)
    cursor.execute(query, args)
    conn.commit()
    conn.close()

if __name__ == "__main__":
    print("This should be imported")
    print("In the future, I might add some export functionality here or something")
