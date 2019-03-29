#!/usr/bin/env python3

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.widgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.figure
import numpy as np
import os
import pandas as pd
import PyQt5.QtWidgets as QtWidgets
import PyQt5.QtGui as QtGui
import PyQt5.QtCore as QtCore
import PyQt5.Qt as Qt
import scipy.spatial
import skimage.draw
import skimage.measure
import skimage.morphology
import uuid

from . import database

class Analyser:
    def __init__(self, experiment_view, experiment_data, image_loader):
        self.ASSIGNING = False
        self.experiment_view = experiment_view
        self._data = experiment_data
        self.image_loader = image_loader
        if not database.checkTable("nuclei"):
            database.createNucleiTable()

    def set_screen_res(self, max_width_px, max_height_px, screen_dpi):
        self.max_width_px = max_width_px
        self.max_height_px = max_height_px
        self.screen_dpi = screen_dpi

    def construct_box(self, layout):
        self.get_outline_objects()
        self.main_layout = layout

        self.grid_layout = QtWidgets.QGridLayout()
        self.grid_layout.setHorizontalSpacing(100)

        self.grid_layout.addWidget(QtWidgets.QLabel("Cells observed from birth to division:"), 0, 0)
        self.grid_layout.addWidget(QtWidgets.QLabel("{0}".format(
            len(self.cells),
        )), 0, 1)

        self.grid_layout.addWidget(QtWidgets.QLabel("Outlines:"), 1, 0)
        self.grid_layout.addWidget(QtWidgets.QLabel("{0}".format(
            len(self.outlines),
        )), 1, 1)

        background_signals = self.calculate_signal(self.wildtype_outlines, stat=None, channel=2)
        self.background = np.mean(background_signals)

        self.grid_layout.addWidget(QtWidgets.QLabel("Wildtype outlines:"), 2, 0)
        self.grid_layout.addWidget(QtWidgets.QLabel("{0} ({1})".format(
            len(self.wildtype_outlines),
            len(background_signals),
        )), 2, 1)
        btn_row = QtWidgets.QHBoxLayout()
        if len(background_signals) != len(self.wildtype_outlines):
            recalc = QtWidgets.QPushButton("Re-calculate background")
            recalc.clicked.connect(lambda: self.recalculate_background())
            btn_row.addWidget(recalc)

        wildtype_associate = QtWidgets.QPushButton("Import wildtype outlines")
        wildtype_associate.clicked.connect(lambda: self.associate_wildtype())
        btn_row.addWidget(wildtype_associate)
        self.grid_layout.addLayout(btn_row, 2, 2)

        self.grid_layout.addWidget(QtWidgets.QLabel("Nuclei assigned:"), 3, 0)
        nuclei = pd.DataFrame(database.getNucleiByExperimentId(self._data.experiment_id))
        if len(nuclei) == 0:
            nuclei = pd.DataFrame(columns=[x[0] for x in database.NucleusRow.COLS])

        nuclei = nuclei[nuclei.outline_id.isin(self.outlines.outline_id.unique())]
        num_nuclei = len(nuclei.outline_id.unique())
        num_nuclei_label = QtWidgets.QLabel("{0}".format(
            num_nuclei,
        ))
        self.grid_layout.addWidget(num_nuclei_label, 3, 1)
        if num_nuclei < len(self.outlines):
            num_nuclei_label.setText("<font color='red'>{0}</font>".format(
                num_nuclei,
            ))
            nuclei_btn = QtWidgets.QPushButton("Assign")
            nuclei_btn.clicked.connect(lambda: self.assign_nuclei(nuclei))
            self.grid_layout.addWidget(nuclei_btn, 3, 2)
        else:
            btn_row = QtWidgets.QHBoxLayout()
            reassign = QtWidgets.QPushButton("Re-assign")
            reassign.clicked.connect(lambda: self.assign_nuclei())
            btn_row.addWidget(reassign)
            nuclei_btn = QtWidgets.QPushButton("Verify")
            nuclei_btn.clicked.connect(lambda: self.verify_nuclei(nuclei))
            btn_row.addWidget(nuclei_btn)
            self.grid_layout.addLayout(btn_row, 3, 2)

        export_btn = QtWidgets.QPushButton("Process data")
        export_btn.clicked.connect(lambda: self.export_data())
        self.grid_layout.addWidget(export_btn, 4, 2)

        self.main_layout.addLayout(self.grid_layout)

    def get_outline_objects(self):
        self.cells = pd.DataFrame(database.getCellsByExperimentId(
            self._data.experiment_id,
            birth_observed=True,
            division_observed=True,
            is_wildtype=False,
        ))
        if len(self.cells) == 0:
            self.cells = pd.DataFrame(columns=[x[0] for x in database.CellRow.COLS])

        wildtype = pd.DataFrame(database.getCellsByExperimentId(
            self._data.experiment_id,
            is_wildtype=True,
        ))
        if len(wildtype) == 0:
            wildtype = pd.DataFrame(columns=[x[0] for x in database.CellRow.COLS])

        all_outlines = pd.DataFrame(database.getOutlinesByExperimentId(
            self._data.experiment_id,
        ))
        self.outlines = all_outlines[
            all_outlines.cell_id.isin(self.cells.cell_id.unique())
        ]
        self.wildtype_outlines = all_outlines[
            all_outlines.cell_id.isin(wildtype.cell_id.unique())
        ]

        # get associations
        associations = database.getAssociationsByExperimentId(self._data.experiment_id, "wildtype")
        for association in associations:
            wt = pd.DataFrame(database.getCellsByExperimentId(
                association.associated_experiment_id,
                is_wildtype=True,
            ))
            if len(wt) == 0:
                wt = pd.DataFrame(columns=[x[0] for x in database.CellRow.COLS])

            assoc_outlines = pd.DataFrame(database.getOutlinesByExperimentId(
                association.associated_experiment_id
            ))
            wt_outlines = assoc_outlines[
                assoc_outlines.cell_id.isin(wt.cell_id.unique())
            ]
            self.wildtype_outlines = self.wildtype_outlines.append(wt_outlines, ignore_index=True)

    def export_data(self):
        self.c2_background = self.calculate_signal(self.wildtype_outlines, channel=1)
        self.c3_background = self.background
        self.outlines["num_nuclei"] = 0
        self.outlines["time_h"] = 0
        self.outlines["mitosis_h"] = 0
        self.outlines["birth_h"] = 0
        self.outlines["cell_area"] = 0

        self.outlines["total_cell_c2"] = 0
        self.outlines["total_cell_c3"] = 0
        self.outlines["norm_cell_c2"] = 0
        self.outlines["norm_cell_c3"] = 0
        self.outlines["total_nuclear1_c2"] = 0
        self.outlines["total_nuclear1_c3"] = 0
        self.outlines["total_nuclear2_c2"] = 0
        self.outlines["total_nuclear2_c3"] = 0
        self.outlines["norm_nuclear1_c2"] = 0
        self.outlines["norm_nuclear1_c3"] = 0
        self.outlines["norm_nuclear2_c2"] = 0
        self.outlines["norm_nuclear2_c3"] = 0
        self.outlines["total_nuclear_area"] = 0
        self.outlines["nuclear_area1"] = 0
        self.outlines["nuclear_area2"] = 0

        spacer = QtWidgets.QLabel("")
        self.main_layout.addWidget(spacer)

        progress_bar = QtWidgets.QProgressBar()
        progress_bar.setMaximum(100)
        progress_bar.setMinimum(0)
        progress_bar.setValue(0)
        self.main_layout.addWidget(progress_bar)

        progress_text = QtWidgets.QLabel("Exporting cell data ({0}/{1})".format(
            0, len(self.cells),
        ))
        self.main_layout.addWidget(progress_text)
        progress_text.show()
        progress_bar.show()
        spacer.show()

        ordered_outlines = []
        cell_data = []
        i = 1
        for _, cell in self.cells.iterrows():
            cell_outlines = self.get_cell_outlines(cell.cell_id)
            ordered_outlines.extend(list(cell_outlines.outline_id))

            out_data = cell.to_dict()
            out_data["start_frame"] = cell.start_frame_idx + 1
            out_data["end_frame"] = cell.end_frame_idx + 1

            for var in ["start_frame_idx", "end_frame_idx", "birth_observed",
                        "division_observed", "cell_num", "is_wildtype"]:
                del out_data[var]

            out_data["growth_rate"] = self.get_growth_rate(cell_outlines)
            out_data["interdivision_time"] = self.get_interdivision_time(cell_outlines)
            out_data["birth_area"] = cell_outlines.iloc[0].cell_area
            out_data["division_area"] = cell_outlines.iloc[-1].cell_area
            cell_data.append(out_data)

            progress_bar.setValue(100 * i / len(self.cells))
            progress_text.setText("Exporting cell data ({0}/{1})".format(
                i, len(self.cells),
            ))
            i += 1

        cell_data = pd.DataFrame(
            cell_data,
            columns=[
                "cell_id", "experiment_id", "parent_cell_id", "child_cell_id1",
                "child_cell_id2", "first_outline_id", "last_outline_id",
                "start_frame", "end_frame", "interdivision_time", "birth_area",
                "division_area", "growth_rate",
            ]
        )
        cell_data_path = os.path.join(
            "data", "output", "{0}-{1}-{2}-cell-data-{3}.xlsx".format(
                self._data.date,
                self._data.strain,
                self._data.medium,
                self._data.experiment_id,
            )
        )
        if not os.path.exists(os.path.dirname(cell_data_path)):
            os.makedirs(os.path.dirname(cell_data_path))

        cell_data.to_excel(cell_data_path)

        outline_data = []
        progress_bar.setValue(0)
        progress_text.setText("Exporting outline data ({0}/{1})".format(
            0, len(self.outlines),
        ))
        for outline_id in ordered_outlines:
            outline = self.outlines[self.outlines.outline_id == outline_id].iloc[0]

            out_data = outline.to_dict()
            for unwanted_col in [
                "coords_path", "experiment_num", "frame_idx", "offset_left",
                "offset_top", "outline_num",
            ]:
                del out_data[unwanted_col]
            out_data["frame"] = outline.frame_idx + 1
            outline_data.append(out_data)

        outline_data = pd.DataFrame(
            outline_data,
            columns=[
                "outline_id", "cell_id", "experiment_id", "parent_id",
                "child_id1", "child_id2", "image_path", "coords_path",
                "time_h", "birth_h", "mitosis_h", "frame", "cell_area",
                "total_cell_c2", "total_cell_c3", "norm_cell_c2",
                "norm_cell_c3", "num_nuclei", "total_nuclear_area",
                "total_nuclear1_c2", "total_nuclear1_c3", "norm_nuclear1_c2",
                "norm_nuclear1_c3", "nuclear_area1", "total_nuclear2_c2",
                "total_nuclear2_c3", "norm_nuclear2_c2", "norm_nuclear2_c3",
                "nuclear_area2",
            ]
        )

        outline_data_path = os.path.join(
            "data", "output", "{0}-{1}-{2}-outline-data-{3}.xlsx".format(
                self._data.date,
                self._data.strain,
                self._data.medium,
                self._data.experiment_id,
            )
        )
        if not os.path.exists(os.path.dirname(outline_data_path)):
            os.makedirs(os.path.dirname(outline_data_path))

        outline_data.to_excel(outline_data_path)
        self.experiment_view._refreshLayout()

    def get_cell_outlines(self, cell_id):
        cell_outlines = self.outlines[self.outlines.cell_id == cell_id].sort_values("frame_idx")
        for _, outline in cell_outlines.iterrows():
            nuclei = database.getNucleiByOutlineId(outline.outline_id)
            num_nuclei = len(nuclei)
            time_h = outline.frame_idx / 6
            birth_h = (outline.frame_idx - cell_outlines.iloc[0].frame_idx) / 6
            cell_area = self.get_cell_area(outline)

            c2 = self.image_loader.load_frame(outline.frame_idx, 1) - self.c2_background
            c3 = self.image_loader.load_frame(outline.frame_idx, 2) - self.c3_background
            outline_coords = np.load(outline.coords_path) + np.array([
                outline.offset_left,
                outline.offset_top,
            ])
            c2_signal = self.get_signal_from_coords(outline_coords, c2)[0]
            c3_signal = self.get_signal_from_coords(outline_coords, c3)[0]
            total_cell_c2 = c2_signal
            total_cell_c3 = c3_signal
            norm_cell_c2 = c2_signal / cell_area
            norm_cell_c3 = c3_signal / cell_area

            total_nuclear_area = 0
            for i, nucleus in enumerate(nuclei):
                this_nuclear_coords = np.load(nucleus.coords_path)
                this_nuclear_area = self.get_cell_area(this_nuclear_coords, is_coords=True)
                total_nuclear_area += this_nuclear_area
                c2_signal = self.get_signal_from_coords(this_nuclear_coords, c2)[0]
                c3_signal = self.get_signal_from_coords(this_nuclear_coords, c3)[0]
                if i == 0:
                    total_nuclear1_c2 = c2_signal
                    total_nuclear1_c3 = c3_signal
                    norm_nuclear1_c2 = c2_signal / this_nuclear_area
                    norm_nuclear1_c3 = c3_signal / this_nuclear_area
                    nuclear_area1 = this_nuclear_area
                elif i == 1:
                    total_nuclear2_c2 = c2_signal
                    total_nuclear2_c3 = c3_signal
                    norm_nuclear2_c2 = c2_signal / this_nuclear_area
                    norm_nuclear2_c3 = c3_signal / this_nuclear_area
                    nuclear_area2 = this_nuclear_area

            if len(nuclei) == 1:
                total_nuclear2_c2 = None
                total_nuclear2_c3 = None
                norm_nuclear2_c2 = None
                norm_nuclear2_c3 = None
                nuclear_area2 = None

            for target in [cell_outlines, self.outlines]:
                stuff = [
                    ("num_nuclei", num_nuclei),
                    ("time_h", time_h),
                    ("birth_h", birth_h),
                    ("cell_area", cell_area),
                    ("total_cell_c2", total_cell_c2),
                    ("total_cell_c3", total_cell_c3),
                    ("norm_cell_c2", norm_cell_c2),
                    ("norm_cell_c3", norm_cell_c3),
                    ("total_nuclear1_c2", total_nuclear1_c2),
                    ("total_nuclear1_c3", total_nuclear1_c3),
                    ("norm_nuclear1_c2", norm_nuclear1_c2),
                    ("norm_nuclear1_c3", norm_nuclear1_c3),
                    ("nuclear_area1", nuclear_area1),
                    ("total_nuclear2_c2", total_nuclear2_c2),
                    ("total_nuclear2_c3", total_nuclear2_c3),
                    ("norm_nuclear2_c2", norm_nuclear2_c2),
                    ("norm_nuclear2_c3", norm_nuclear2_c3),
                    ("nuclear_area2", nuclear_area2),
                    ("total_nuclear_area", total_nuclear_area),
                ]
                target.loc[_, [x[0] for x in stuff]] = [x[1] for x in stuff]

        two_nucl = cell_outlines[cell_outlines.num_nuclei >= 2].time_h.min()
        if time_h - two_nucl > 1.5:
            two_nucl = cell_outlines[(
                (cell_outlines.num_nuclei >= 2) &
                (cell_outlines.time_h > cell_outlines.time_h.iloc[-1] - 1.5)
            )].time_h.min()

        mitosis_h = cell_outlines.time_h - two_nucl
        cell_outlines["mitosis_h"] = mitosis_h
        for _, outline in cell_outlines.iterrows():
            self.outlines.loc[_, ["mitosis_h"]] = outline.mitosis_h

        return cell_outlines

    def get_growth_rate(self, outlines):
        pre_mitotic = outlines[outlines.mitosis_h <= 0]
        try:
            gamma, log_A0 = np.polyfit(pre_mitotic.birth_h, np.log(pre_mitotic.cell_area), 1)
        except:
            print("ERROR:")
            print("pre_mitotic:")
            print(pre_mitotic)
            print("birth_h:", pre_mitotic.birth_h)
            print("cell_area:", pre_mitotic.cell_area)
            return None
        else:
            return gamma

    def get_interdivision_time(self, outlines):
        return outlines.birth_h.iloc[-1]

    def get_cell_area(self, outline, is_coords=False):
        if not is_coords:
            coords = np.load(outline.coords_path)
        else:
            coords = outline

        px = self.image_loader.get_pixel_conversion()
        area = 0.5 * np.abs(
            np.dot(coords[:, 0], np.roll(coords[:, 1], 1)) -
            np.dot(coords[:, 1], np.roll(coords[:, 0], 1))
        ) * px * px
        return area

    def recalculate_background(self):
        self.background = self.calculate_signal(self.wildtype_outlines, channel=2, replace=True, show_progress=True)
        self.experiment_view._refreshLayout()

    def get_signal_from_coords(self, coords, img):
        rr, cc = skimage.draw.polygon(
            coords[:, 0], coords[:, 1]
        )
        blank = np.zeros_like(img)
        blank[rr, cc] = img[rr, cc]
        signal = blank.sum()
        num_pixels = len(rr)
        return signal, num_pixels

    def calculate_signal(self, outlines, stat=np.mean, channel=2, replace=False, show_progress=False, save=True):
        bg_path = os.path.join("data", "signals", "C{0}-{1}.npy".format(
            channel,
            self._data.experiment_id
        ))
        if os.path.exists(bg_path) and not replace:
            signals = np.load(bg_path)
        else:
            if show_progress:
                progress_bar = QtWidgets.QProgressBar()
                progress_bar.setMaximum(100)
                progress_bar.setMinimum(0)
                progress_bar.setValue(0)
                self.main_layout.addWidget(QtWidgets.QLabel(""))
                self.main_layout.addWidget(progress_bar)
                progress_bar.show()
                comment = QtWidgets.QLabel("Determining background ({0}/{1})".format(
                    0, len(outlines),
                ))
                self.main_layout.addWidget(comment)
                comment.show()

            if not os.path.exists(os.path.dirname(bg_path)):
                os.makedirs(os.path.dirname(bg_path))

            signals = np.zeros(len(outlines))

            i = 0
            for _, outline in outlines.iterrows():
                img = self.image_loader.load_frame(outline.frame_idx, channel)
                c = np.load(outline.coords_path) + np.array([
                    outline.offset_left,
                    outline.offset_top,
                ])
                signal, num_pixels = self.get_signal_from_coords(c, img)
                signals[i] = signal / num_pixels
                i += 1
                if show_progress:
                    progress_bar.setValue(100 * i / len(outlines))
                    comment.setText("Determining background ({0}/{1})".format(
                        i, len(outlines)
                    ))

            if show_progress:
                progress_bar.hide()
                progress_bar.deleteLater()
                comment.hide()
                comment.deleteLater()

            if len(signals) == 0:
                signals = np.array([0])

            if save:
                np.save(bg_path, signals)

        if stat is not None:
            return stat(signals)
        else:
            return signals

    def assign_nuclei(self, existing=None):
        if self.ASSIGNING:
            return
        else:
            self.ASSIGNING = True

        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setMaximum(100)
        self.progress_bar.setMinimum(0)
        self.progress_bar.setValue(0)
        self.main_layout.addWidget(QtWidgets.QLabel(""))
        self.main_layout.addWidget(self.progress_bar)
        self.progress_bar.show()

        comment = QtWidgets.QLabel("Defining nuclei ({0}/{1})".format(
            0, len(self.outlines),
        ))
        self.main_layout.addWidget(comment)
        comment.show()

        if existing is None:
            database.deleteNucleiByExperimentId(self._data.experiment_id)

        i = 1
        for _, outline in self.outlines.iterrows():
            if existing is None or (existing is not None and len(existing[existing.outline_id == outline.outline_id]) == 0):
                self.define_nucleus(outline)

            comment.setText("Defining nuclei ({0}/{1})".format(
                i, len(self.outlines),
            ))
            self.progress_bar.setValue(100 * i / len(self.outlines))
            i += 1

        self.ASSIGNING = False
        self.experiment_view._refreshLayout()

    def define_nucleus(self, outline):
        c3 = self.image_loader.load_frame(outline.frame_idx, 2) - self.background
        c = np.load(outline.coords_path) + np.array([
            outline.offset_left,
            outline.offset_top,
        ])
        cell_rr, cell_cc = skimage.draw.polygon(
            c[:, 0], c[:, 1]
        )

        blank = np.zeros_like(c3)
        blank[cell_rr, cell_cc] = 1
        threshold = c3[cell_rr, cell_cc].mean() + 1.1 * c3[cell_rr, cell_cc].std()
        blank[c3 < threshold] = 0
        blank = skimage.morphology.dilation(blank)
        labels = skimage.measure.label(blank)
        labels = skimage.morphology.remove_small_objects(labels, 100)
        px = self.image_loader.get_pixel_conversion()

        nuclei = []
        for label in np.unique(labels):
            if label == 0:
                continue
            image = np.zeros_like(labels)
            image[labels == label] = 1
            image_coords = np.transpose(np.nonzero(image))
            hull = scipy.spatial.ConvexHull(image_coords)
            coords = image_coords[hull.vertices]
            nucleus_id = str(uuid.uuid4())
            coords_path = os.path.join(
                "data", "nuclei",
                self._data.experiment_id,
                outline.cell_id,
                outline.outline_id,
                "{0}.npy".format(nucleus_id),
            )
            nuclei.append({
                "nucleus_id": nucleus_id,
                "outline_id": outline.outline_id,
                "cell_id": outline.cell_id,
                "experiment_id": outline.experiment_id,
                "coords": coords,
                "coords_path": coords_path,
                "nuclear_area": 0.5 * np.abs(
                    np.dot(coords[:, 0], np.roll(coords[:, 1], 1)) -
                    np.dot(coords[:, 1], np.roll(coords[:, 0], 1))
                ) * px * px,
            })

        if len(nuclei) > 2:
            # take the largest 2
            nuclei = sorted(
                nuclei,
                key=lambda x: x["nuclear_area"],
                reverse=True,
            )[:2]


        if not os.path.exists(os.path.dirname(coords_path)):
            os.makedirs(os.path.dirname(coords_path))

        for nucleus in nuclei:
            np.save(nucleus["coords_path"], nucleus["coords"])
            database.insertNucleus(
                nucleus["nucleus_id"],
                nucleus["outline_id"],
                nucleus["cell_id"],
                nucleus["experiment_id"],
                nucleus["coords_path"],
            )

    def verify_nuclei(self, nuclei):
        verifier = NuclearVerifier(
            self.experiment_view.window,
            nuclei,
            self.max_width_px,
            self.max_height_px,
            self.screen_dpi,
            self.image_loader,
            self.background
        )
        verifier.create_layout()

    def associate_wildtype(self):
        if not database.checkTable("associations"):
            database.createAssociationsTable()

        experiments = database.getExperiments()
        self.associations = database.getAssociationsByExperimentId(self._data.experiment_id, "wildtype")
        self.associated_experiment_ids = [x.associated_experiment_id for x in self.associations]

        assoc_window = QtWidgets.QDialog(self.experiment_view.window)
        assoc_window.setModal(True)
        assoc_window.setWindowTitle("Associate wildtype outlines from another experiment")
        assoc_window.setGeometry(0, 60, 700, 200)

        table = QtWidgets.QTableWidget()
        table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        table.setRowCount(len(experiments) - 1)
        table_headers = ["Date", "Medium", "Strain", "Wildtype outlines", "Included"]
        table.setColumnCount(len(table_headers))
        table.setHorizontalHeaderLabels(table_headers)

        row_num = 0
        table_items = []
        for exp in experiments:
            row_items = []
            if exp.experiment_id == self._data.experiment_id:
                continue

            for col_num, var in enumerate(["date", "medium", "strain", "wildtype"]):
                item = QtWidgets.QTableWidgetItem()
                item.setFlags(QtCore.Qt.ItemIsEnabled)

                if var == "wildtype":
                    # calculate number of wildtype cells
                    wildtype = set([
                        x.cell_id
                        for x in database.getCellsByExperimentId(
                            exp.experiment_id,
                            is_wildtype=True,
                        )
                    ])
                    outlines = database.getOutlinesByExperimentId(
                        exp.experiment_id,
                    )
                    wt = 0
                    for out in outlines:
                        if out.cell_id in wildtype:
                            wt += 1
                    item.setText(str(wt))
                else:
                    item.setText(str(exp[var]))

                item.setTextAlignment(QtCore.Qt.AlignCenter)
                if exp.experiment_id in self.associated_experiment_ids:
                    item.setData(
                        QtCore.Qt.BackgroundRole,
                        QtGui.QBrush(QtGui.QColor("green"))
                    )
                else:
                    item.setData(
                        QtCore.Qt.BackgroundRole,
                        QtGui.QBrush(QtGui.QColor("white"))
                    )

                table.setItem(row_num, col_num, item)
                row_items.append(item)

            include_btn = QtWidgets.QPushButton(table)
            include_btn._experiment_id = exp.experiment_id
            include_btn._row_num = row_num
            include_btn.setText("Include")

            def _set_association():
                sender = self.experiment_view.window.sender()
                experiment_id = sender._experiment_id
                row_num = sender._row_num
                if experiment_id in self.associated_experiment_ids:
                    association_id = [
                        x.association_id
                        for x in self.associations
                        if x.associated_experiment_id == experiment_id
                    ][0]
                    database.deleteAssociationById(association_id)
                    for col_num in range(len(table_headers) - 1):
                        table_items[row_num][col_num].setData(
                            QtCore.Qt.BackgroundRole,
                            QtGui.QBrush(QtGui.QColor("white"))
                        )
                else:
                    association_id = str(uuid.uuid4())
                    database.insertAssociation(
                        association_id,
                        self._data.experiment_id,
                        experiment_id,
                        "wildtype",
                    )
                    for col_num in range(len(table_headers) - 1):
                        table_items[row_num][col_num].setData(
                            QtCore.Qt.BackgroundRole,
                            QtGui.QBrush(QtGui.QColor("green"))
                        )

                self.associations = database.getAssociationsByExperimentId(self._data.experiment_id, "wildtype")
                self.associated_experiment_ids = [x.associated_experiment_id for x in self.associations]

            include_btn.clicked.connect(_set_association)
            table.setCellWidget(row_num, 4, include_btn)
            table_items.append(row_items)
            row_num += 1

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(table)

        btns = QtWidgets.QHBoxLayout()
        ok_btn = QtWidgets.QPushButton("OK")
        ok_btn.clicked.connect(lambda: assoc_window.close())
        btns.addWidget(ok_btn)
        layout.addLayout(btns)

        assoc_window.setLayout(layout)
        assoc_window.show()
        assoc_window.finished.connect(lambda: self.experiment_view._refreshLayout())


class NuclearVerifier:
    def __init__(self, parent_window, nuclei, max_width_px, max_height_px, screen_dpi, image_loader, background):
        self.parent_window = parent_window
        self.nuclei = nuclei
        self.max_width_px = max_width_px
        self.max_height_px = max_height_px
        self.screen_dpi = screen_dpi
        self.image_loader = image_loader
        self.background = background

        self.window = QtWidgets.QDialog(self.parent_window)
        self.window.setModal(True)
        self.window.setWindowTitle("Verify nuclei")
        self.window.setGeometry(0, 60, 0.9 * self.max_width_px, 0.9 * self.max_height_px)
        self.window.show()

    def create_layout(self):
        self.main_layout = QtWidgets.QVBoxLayout()
        lineage_box = QtWidgets.QWidget()
        lineage_layout = QtWidgets.QVBoxLayout()

        for cell_id in self.nuclei.cell_id.unique():
            outlines = sorted(
                database.getOutlinesByCellId(cell_id),
                key=lambda x: x.frame_idx,
            )
            nuclei = pd.DataFrame(database.getNucleiByCellId(cell_id))

            width = self.max_width_px * 0.1
            cell_plot = Plotter(
                self.window,
                width=(width * len(outlines)) / self.screen_dpi,
                height=width / self.screen_dpi,
                dpi=self.screen_dpi,
                subplots=len(outlines),
            )
            cell_plot.mpl_connect("button_press_event", lambda x: print(x.button, x.inaxes))
            for ax, outline in zip(cell_plot.axes, outlines):
                c3 = self.image_loader.load_frame(outline.frame_idx, 2) - self.background
                im = c3[
                    outline.offset_left:outline.offset_left + 150,
                    outline.offset_top:outline.offset_top + 150,
                ]
                ax.imshow(im, cmap="binary")
                c = np.load(outline.coords_path)
                outline_poly = matplotlib.patches.Polygon(
                    np.array([c[:, 1], c[:, 0]]).T,
                    edgecolor="b",
                    fill=False,
                    lw=1,
                    linestyle="--",
                )
                ax.add_patch(outline_poly)

                outline_nuclei = nuclei[nuclei.outline_id == outline.outline_id]
                if len(outline_nuclei) > 2:
                    c = "r"
                else:
                    c = "y"

                for _, nucleus in outline_nuclei.iterrows():
                    n = np.load(nucleus.coords_path) - np.array([
                        outline.offset_left,
                        outline.offset_top,
                    ])
                    n_poly = matplotlib.patches.Polygon(
                        np.array([n[:, 1], n[:, 0]]).T,
                        edgecolor=c,
                        fill=False,
                        lw=1,
                    )
                    ax.add_patch(n_poly)

            cell_scroll_area = QtWidgets.QScrollArea()
            cell_scroll_area.verticalScrollBar().setEnabled(False)
            cell_scroll_area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            cell_scroll_area.setAlignment(QtCore.Qt.AlignLeft)
            cell_scroll_area.setWidget(cell_plot)

            lineage_layout.addWidget(cell_scroll_area)

        lineage_box.setLayout(lineage_layout)
        lineage_box.setMinimumWidth(0.9 * self.max_width_px - 40)
        lineage_box.setMaximumWidth(0.9 * self.max_width_px - 40)

        lineage_scroll_area = QtWidgets.QScrollArea()
        lineage_scroll_area.horizontalScrollBar().setEnabled(False)
        lineage_scroll_area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        lineage_scroll_area.setWidget(lineage_box)
        self.main_layout.addWidget(lineage_scroll_area)

        self.window.setLayout(self.main_layout)


class Plotter(FigureCanvas):
    def __init__(self, parent_window, width, height, dpi, subplots=3):
        fig = matplotlib.figure.Figure(figsize=(width, height), dpi=dpi)
        self.axes = []
        self.offsets = []
        for sp in range(subplots):
            ax = fig.add_subplot(1, subplots, sp + 1)
            # ax.axis("off")
            ax.set_xticks([], [])
            ax.set_yticks([], [])
            ax.set_aspect("equal")
            ax.autoscale("off")
            self.axes.append(ax)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent_window)

        FigureCanvas.setSizePolicy(
            self,
            QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding,
        )
        FigureCanvas.updateGeometry(self)
        fig.tight_layout()
