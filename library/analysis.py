#!/usr/bin/env python3

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
        self.main_layout = layout

        self.grid_layout = QtWidgets.QGridLayout()
        self.grid_layout.setHorizontalSpacing(100)

        cells = pd.DataFrame(database.getCellsByExperimentId(
            self._data.experiment_id,
            birth_observed=True,
            division_observed=True,
            is_wildtype=False,
        ))
        if len(cells) == 0:
            cells = pd.DataFrame(columns=[x[0] for x in database.CellRow.COLS])

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
            all_outlines.cell_id.isin(cells.cell_id.unique())
        ]
        self.wildtype_outlines = all_outlines[
            all_outlines.cell_id.isin(wildtype.cell_id.unique())
        ]

        self.grid_layout.addWidget(QtWidgets.QLabel("Cells observed from birth to division:"), 0, 0)
        self.grid_layout.addWidget(QtWidgets.QLabel("{0}".format(
            len(cells),
        )), 0, 1)

        self.grid_layout.addWidget(QtWidgets.QLabel("Outlines (total):"), 1, 0)
        self.grid_layout.addWidget(QtWidgets.QLabel("{0} ({1})".format(
            len(self.outlines),
            len(all_outlines) - len(self.wildtype_outlines),
        )), 1, 1)

        self.grid_layout.addWidget(QtWidgets.QLabel("Nuclei assigned:"), 2, 0)
        nuclei = pd.DataFrame(database.getNucleiByExperimentId(self._data.experiment_id))
        if len(nuclei) == 0:
            nuclei = pd.DataFrame(columns=[x[0] for x in database.NucleusRow.COLS])

        nuclei = nuclei[nuclei.outline_id.isin(self.outlines.outline_id.unique())]
        num_nuclei = len(nuclei.outline_id.unique())
        num_nuclei_label = QtWidgets.QLabel("{0}".format(
            num_nuclei,
        ))
        self.grid_layout.addWidget(num_nuclei_label, 2, 1)
        if num_nuclei < len(self.outlines):
            num_nuclei_label.setText("<font color='red'>{0}</font>".format(
                num_nuclei,
            ))
            nuclei_btn = QtWidgets.QPushButton("Assign")
            nuclei_btn.clicked.connect(lambda: self.assign_nuclei(nuclei))
            self.grid_layout.addWidget(nuclei_btn, 2, 2)
        else:
            btn_row = QtWidgets.QHBoxLayout()
            reassign = QtWidgets.QPushButton("Re-assign")
            reassign.clicked.connect(lambda: self.assign_nuclei())
            btn_row.addWidget(reassign)
            nuclei_btn = QtWidgets.QPushButton("Verify")
            nuclei_btn.clicked.connect(lambda: self.verify_nuclei())
            btn_row.addWidget(nuclei_btn)
            self.grid_layout.addLayout(btn_row, 2, 3)

        self.grid_layout.addWidget(QtWidgets.QLabel("Wildtype outlines (cells):"), 3, 0)
        self.grid_layout.addWidget(QtWidgets.QLabel("{0} ({1})".format(
            len(self.wildtype_outlines),
            len(wildtype),
        )))

        self.main_layout.addLayout(self.grid_layout)

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
        c3 = self.image_loader.load_frame(outline.frame_idx, 2)
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
            if not os.path.exists(os.path.dirname(coords_path)):
                os.makedirs(os.path.dirname(coords_path))

            np.save(coords_path, coords)
            database.insertNucleus(
                nucleus_id,
                outline.outline_id,
                outline.cell_id,
                outline.experiment_id,
                coords_path,
            )

    def verify_nuclei(self):
        print("Verify")
