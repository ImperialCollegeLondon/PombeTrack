#!/usr/bin/env python3

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.widgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import PyQt5.QtWidgets as QtWidgets
import PyQt5.QtGui as QtGui
import PyQt5.QtCore as QtCore
import PyQt5.Qt as Qt
import seaborn as sns
import os
import tifffile
import uuid

from . import database

sns.set_context("talk")
sns.set_style("white")

class Plotter(FigureCanvas):
    def __init__(self, parent_window, width, height, dpi, experiment_data, image_loader, subplots=3):
        fig = matplotlib.figure.Figure(figsize=(width, height), dpi=dpi)
        self.axes = []
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

        self.screen_dpi = dpi
        self._data = experiment_data
        self.image_loader = image_loader
        fig.tight_layout()


class Assigner:
    def __init__(self, experiment_data, image_loader):
        self.experiment_data = experiment_data
        self.image_loader = image_loader
        self.region_width = 75
        self.region_height = 75

    def set_screen_res(self, max_width_px, max_height_px, screen_dpi):
        self.max_width_px = max_width_px
        self.max_height_px = max_height_px
        self.screen_dpi = screen_dpi

    def _px_to_in(self, px):
        return px / self.screen_dpi

    def start_assigning(self, parent_window):
        if not hasattr(self, "max_width_px"):
            raise ValueError("Screen resolution has not been set")

        self.parent_window = parent_window
        self.window = QtWidgets.QDialog(self.parent_window)
        self.window.setModal(True)
        self.window.setGeometry(0, 0, self.max_width_px * 0.9, self.max_height_px * 0.9)
        self.window.setWindowTitle("Assign/Verify cell lineages")

        main_layout = QtWidgets.QVBoxLayout()

        menubar = QtWidgets.QMenuBar(self.window)
        file_menu = menubar.addMenu("&File")
        quit_action = QtWidgets.QAction("&Close", menubar)
        quit_action.triggered[bool].connect(lambda: self.window.close())
        file_menu.addAction(quit_action)
        main_layout.setMenuBar(menubar)

        outlines = database.getOutlinesByExperimentId(self.experiment_data.experiment_id)
        self.outlines = pd.DataFrame(outlines)
        unique_cells = self.outlines.cell_id.unique()

        lineage_layout = QtWidgets.QVBoxLayout()
        for cell_num, cell_id in enumerate(unique_cells):
            # cell_box = QtWidgets.QGroupBox("Cell #{0} ({1})".format(cell_num + 1, cell_id))
            cell_box = QtWidgets.QWidget()
            cell_layout = QtWidgets.QVBoxLayout()

            cell_title_layout = QtWidgets.QHBoxLayout()
            cell_title_layout.addWidget(QtWidgets.QLabel("Cell #{0} ({1})".format(cell_num + 1, cell_id)))
            verify_btn = QtWidgets.QPushButton("Assign Cell Lineage")
            verify_btn._cell_id = cell_id
            verify_btn.clicked.connect(self.assign_lineage)
            cell_title_layout.addWidget(verify_btn)
            cell_layout.addLayout(cell_title_layout)

            cell_outlines = self.outlines[self.outlines.cell_id == cell_id]
            width = len(cell_outlines) * self.max_width_px * 0.1
            cell_plot = Plotter(
                self.window,
                width=self._px_to_in(width),
                height=self._px_to_in(self.max_width_px * 0.1),
                dpi=self.screen_dpi,
                experiment_data=self.experiment_data,
                image_loader=self.image_loader,
                subplots=len(cell_outlines)
            )
            cell_plot.setMinimumWidth(width)
            cell_plot.setMaximumWidth(width)
            cell_plot.setMinimumHeight(self.max_width_px * 0.1)
            # i = 0 
            # for _, outline in cell_outlines.iterrows():
            #     roi = self.image_loader.load_frame(outline.frame_idx, 0)[
            #         outline.offset_left:outline.offset_left + (self.region_width * 2),
            #         outline.offset_top:outline.offset_top + (self.region_height * 2),
            #     ]
            #     cell_plot.axes[i].imshow(roi, cmap="gray")
            #     i += 1

            cell_scroll_area = QtWidgets.QScrollArea()
            cell_scroll_area.verticalScrollBar().setEnabled(False)
            cell_scroll_area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            cell_scroll_area.setAlignment(QtCore.Qt.AlignLeft)
            cell_scroll_area.setWidget(cell_plot)

            cell_layout.addWidget(cell_scroll_area)
            cell_box.setLayout(cell_layout)
            lineage_layout.addWidget(cell_box)

        lineage_widget = QtWidgets.QWidget()
        lineage_widget.setMinimumWidth(self.max_width_px * 0.9 - 50)
        lineage_scroll_area = QtWidgets.QScrollArea()
        lineage_scroll_area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        lineage_widget.setLayout(lineage_layout)
        lineage_scroll_area.setWidget(lineage_widget)

        main_layout.addWidget(lineage_scroll_area)

        self.plot = Plotter(
            self.window,
            width=self._px_to_in(self.max_width_px * 0.9),
            height=self._px_to_in(self.max_height_px * 0.2),
            dpi=self.screen_dpi,
            experiment_data=self.experiment_data,
            image_loader=self.image_loader,
        )
        self.plot.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.plot.setFocus()

        # tool_layout = QtWidgets.QVBoxLayout()
        # main_layout.addLayout(tool_layout)
        main_layout.addWidget(self.plot)

        self.window.setLayout(main_layout)
        self.window.show()

    def _clear_assignment_plot(self):
        for ax in self.plot.axes:
            ax.clear()
            ax.set_xticks([], [])
            ax.set_yticks([], [])
            ax.set_aspect("equal")
            ax.autoscale("off")

    def assign_lineage(self, *args):
        self._clear_assignment_plot()
        data = self.outlines[self.outlines.cell_id == self.window.sender()._cell_id]
        first_cell = data.iloc[0]
        im1 = self.image_loader.load_frame(first_cell.frame_idx, 0)
        self.plot.axes[1].imshow(im1, cmap="gray")
        self.plot.axes[1].set_xlim([
            first_cell.offset_top,
            first_cell.offset_top + (self.region_height * 2),
        ])
        self.plot.axes[1].set_ylim([
            first_cell.offset_left + (self.region_width * 2),
            first_cell.offset_left,
        ])
        c = np.load(first_cell.coords_path) + np.array([
            first_cell.offset_left, first_cell.offset_top
        ])
        p = matplotlib.patches.Polygon(np.array([c[:, 1], c[:, 0]]).T, edgecolor="r", fill=False, lw=1)
        self.plot.axes[1].add_patch(p)


        if len(data) > 1:
            cell_selected = data.iloc[1]
            fidx = data.iloc[1].frame_idx
            offt = data.iloc[1].offset_top
            offl = data.iloc[1].offset_left
            im2 = self.image_loader.load_frame(fidx, 0)
        elif first_cell.frame_idx + 1 < self.image_loader.num_frames:
            cell_selected = None
            fidx = first_cell.frame_idx + 1
            offt = first_cell.offset_top
            offl = first_cell.offset_left
            im2 = self.image_loader.load_frame(fidx, 0)
        else:
            cell_selected = None
            fidx = first_cell.frame_idx + 1
            im2 = np.zeros((self.region_width * 2, self.region_height * 2))
            offt = 0
            offl = 0

        self.plot.axes[2].imshow(im2, cmap="gray")
        self.plot.axes[2].set_xlim([
            offt, offt + (self.region_height * 2)
        ])
        self.plot.axes[2].set_ylim([
            offl + (self.region_width * 2), offl
        ])

        for outline in database.getOutlinesByFrameIdx(fidx, self.experiment_data.experiment_hash):
            c = np.load(outline.coords_path) + np.array([
                outline.offset_left, outline.offset_top
            ])
            if cell_selected is not None and outline.outline_id == cell_selected.outline_id:
                p = matplotlib.patches.Polygon(np.array([c[:, 1], c[:, 0]]).T, edgecolor="k", lw=1, facecolor="g", alpha=0.9)
            else:
                p = matplotlib.patches.Polygon(np.array([c[:, 1], c[:, 0]]).T, edgecolor="k", lw=1, facecolor="y", alpha=0.9)
            self.plot.axes[2].add_patch(p)

        if first_cell.parent_id:
            prev_cell = database.getOutlineById(first_cell.parent_id)
            im3 = self.image_loader.load_frame(prev_cell.frame_idx, 0)
            self.plot.axes[0].imshow(im3, cmap="gray")
            self.plot.axes[0].set_xlim([
                prev_cell.offset_top,
                prev_cell.offset_top + (self.region_height * 2),
            ])
            self.plot.axes[0].set_ylim([
                prev_cell.offset_left + (self.region_width * 2),
                prev_cell.offset_left,
            ])
            c = np.load(prev_cell.coords_path) + np.array([
                prev_cell.offset_left, prev_cell.offset_top
            ])
            p = matplotlib.patches.Polygon(np.array([c[:, 1], c[:, 0]]).T, edgecolor="y", linestyle="--", fill=False, lw=1)
            self.plot.axes[0].add_patch(p)

        elif first_cell.frame_idx - 1 >= 0:
            im3 = self.image_loader.load_frame(first_cell.frame_idx - 1, 0)
            self.plot.axes[0].imshow(im3, cmap="gray")
            self.plot.axes[0].set_xlim([
                first_cell.offset_top,
                first_cell.offset_top + (self.region_height * 2),
            ])
            self.plot.axes[0].set_ylim([
                first_cell.offset_left + (self.region_width * 2),
                first_cell.offset_left,
            ])

        else:
            prev_cell = None
            im3 = np.zeros((self.region_width * 2, self.region_height * 2))
            self.plot.axes[0].imshow(im3, cmap="gray")

        self.plot.draw()


