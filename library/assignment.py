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

        if not database.checkTable("cells"):
            database.createCellsTable()

    def set_screen_res(self, max_width_px, max_height_px, screen_dpi):
        self.max_width_px = max_width_px
        self.max_height_px = max_height_px
        self.screen_dpi = screen_dpi

    def _px_to_in(self, px):
        return px / self.screen_dpi

    def get_outlines(self):
        outlines = database.getOutlinesByExperimentId(self.experiment_data.experiment_id)
        self.outlines = pd.DataFrame(outlines)

    def start_assigning(self, parent_window):
        if not hasattr(self, "max_width_px"):
            raise ValueError("Screen resolution has not been set")

        self.parent_window = parent_window
        self.window = QtWidgets.QDialog(self.parent_window)
        self.window.setModal(True)
        self.window.setGeometry(0, 0, self.max_width_px * 0.5, self.max_height_px * 0.9)
        self.window.setWindowTitle("Assign/Verify cell lineages")
        self.assignment_queue = []

        self.main_layout = QtWidgets.QVBoxLayout()

        menubar = QtWidgets.QMenuBar(self.window)
        file_menu = menubar.addMenu("&File")
        quit_action = QtWidgets.QAction("&Close", menubar)
        quit_action.triggered[bool].connect(lambda: self.window.close())
        file_menu.addAction(quit_action)
        self.main_layout.setMenuBar(menubar)

        self.create_layout()

        self.plot = Plotter(
            self.window,
            width=self._px_to_in(self.max_width_px * 0.5),
            height=self._px_to_in((self.max_width_px * 0.5) / 3),
            dpi=self.screen_dpi,
            experiment_data=self.experiment_data,
            image_loader=self.image_loader,
        )
        self.plot.mpl_connect("key_press_event", self.key_press_event)
        self.plot.mpl_connect("pick_event", self.pick_event)

        # tool_layout = QtWidgets.QVBoxLayout()
        # self.main_layout.addLayout(tool_layout)
        self.main_layout.addWidget(self.plot)

        self.status_bar = QtWidgets.QStatusBar()
        self.main_layout.addWidget(self.status_bar)

        self.window.setLayout(self.main_layout)
        self.window.show()

        self.plot.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.plot.setFocus()

    def create_layout(self):
        self.get_outlines()
        unique_cells = self.outlines.cell_id.unique()
        lineage_layout = QtWidgets.QVBoxLayout()
        verification = (
            database.getCellById(cell_id)
            for cell_id in unique_cells
        )
        unique_cells = sorted(
            zip(unique_cells, verification),
            key=lambda x: bool(x[1]),
        )
        for cell_num, (cell_id, verified_cell) in enumerate(unique_cells):
            # cell_box = QtWidgets.QGroupBox("Cell #{0} ({1})".format(cell_num + 1, cell_id))
            cell_box = QtWidgets.QWidget()
            cell_layout = QtWidgets.QVBoxLayout()
            cell_title_layout = QtWidgets.QHBoxLayout()

            if verified_cell and verified_cell.is_wildtype:
                wildtype = True
                wildtype_btn = QtWidgets.QPushButton("Unset wildtype")
                # desc_str = "Wildtype cell #{0} ({1})".format(
                #     cell_num + 1, cell_id
                # )
                desc_str = "Wildtype cell #{0}".format(
                    cell_num + 1
                )
            else:
                wildtype = False
                wildtype_btn = QtWidgets.QPushButton("Set wildtype")
                # desc_str = "Cell #{0} ({1})".format(
                #     cell_num + 1, cell_id
                # )
                desc_str = "Cell #{0}".format(
                    cell_num + 1
                )

            if verified_cell:
                pixmap = QtGui.QPixmap("resources/tick.png")
            else:
                pixmap = QtGui.QPixmap("resources/cross.png")

            verification_sign = pixmap.scaledToWidth(20)
            verification_label = QtWidgets.QLabel()
            verification_label.setPixmap(verification_sign)
            cell_title_layout.addWidget(verification_label)

            desc_label = QtWidgets.QLabel(desc_str)
            desc_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
            cell_title_layout.addWidget(desc_label)

            verify_btn = QtWidgets.QPushButton("Assign Cell Lineage")
            verify_btn._cell_id = cell_id
            verify_btn.clicked.connect(self.assign_lineage)
            cell_title_layout.addWidget(verify_btn)
            cell_layout.addLayout(cell_title_layout)

            # create proper lineage
            cell_outlines = self.outlines[
                self.outlines.cell_id == cell_id
            ].sort_values("frame_idx")
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
            i = 0
            for _, outline in cell_outlines.iterrows():
                roi = self.image_loader.load_frame(outline.frame_idx, 0)[
                    outline.offset_left:outline.offset_left + (self.region_width * 2),
                    outline.offset_top:outline.offset_top + (self.region_height * 2),
                ]
                cell_plot.axes[i].imshow(roi, cmap="gray")
                cell_plot.axes[i].text(
                    0, 0,
                    "F{0}".format(outline.frame_idx + 1),
                    horizontalalignment="left",
                    verticalalignment="top",
                    color="y",
                    fontweight="bold",
                    fontsize=10,
                )
                i += 1

            cell_scroll_area = QtWidgets.QScrollArea()
            cell_scroll_area.verticalScrollBar().setEnabled(False)
            cell_scroll_area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            cell_scroll_area.setAlignment(QtCore.Qt.AlignLeft)
            cell_scroll_area.setWidget(cell_plot)

            cell_layout.addWidget(cell_scroll_area)
            cell_box.setLayout(cell_layout)
            lineage_layout.addWidget(cell_box)

        lineage_widget = QtWidgets.QWidget()
        lineage_widget.setMinimumWidth(self.max_width_px * 0.5 - 50)
        lineage_widget.setLayout(lineage_layout)
        replace = hasattr(self, "lineage_scroll_area")
        if not replace:
            self.lineage_scroll_area = QtWidgets.QScrollArea()
            self.lineage_scroll_area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)

        self.lineage_scroll_area.setWidget(lineage_widget)

        if not replace:
            self.main_layout.addWidget(self.lineage_scroll_area)

    def _clear_assignment_plot(self):
        for ax in self.plot.axes:
            ax.clear()
            ax.set_xticks([], [])
            ax.set_yticks([], [])
            ax.set_aspect("equal")
            ax.autoscale("off")

    def assign_lineage(self, click=False, outline_id=None):
        if not outline_id:
            potential_outlines = self.outlines[self.outlines.cell_id == self.window.sender()._cell_id]
        else:
            potential_outlines = self.outlines[self.outlines.outline_id == outline_id]

        self.lineage = [potential_outlines[
            potential_outlines.frame_idx == potential_outlines.frame_idx.min()
        ].iloc[0]]
        self.selected_outlines = []
        self.display_frame()
        self.plot.setFocus()

    def display_frame(self):
        self.selected_outlines = []
        self._clear_assignment_plot()
        first_outline = self.lineage[-1]

        self.plot.axes[0].set_title("Frame {0}".format(first_outline.frame_idx))
        self.plot.axes[1].set_title("Frame {0}".format(first_outline.frame_idx + 1))
        self.plot.axes[2].set_title("Frame {0}".format(first_outline.frame_idx + 2))

        im1 = self.image_loader.load_frame(first_outline.frame_idx, 0)
        self.plot.axes[1].imshow(im1, cmap="gray")
        self.plot.axes[1].set_xlim([
            first_outline.offset_top,
            first_outline.offset_top + (self.region_height * 2),
        ])
        self.plot.axes[1].set_ylim([
            first_outline.offset_left + (self.region_width * 2),
            first_outline.offset_left,
        ])
        c = np.load(first_outline.coords_path) + np.array([
            first_outline.offset_left, first_outline.offset_top
        ])
        p = matplotlib.patches.Polygon(
            np.array([c[:, 1], c[:, 0]]).T,
            edgecolor="r",
            fill=False,
            lw=1
        )
        self.plot.axes[1].add_patch(p)

        if first_outline.child_id1 is not None:
            selected_outline = self.outlines[
                self.outlines.outline_id == first_outline.child_id1
            ].iloc[0]
            self.selected_outlines.append(selected_outline)
            fidx = selected_outline.frame_idx
            offt = selected_outline.offset_top
            offl = selected_outline.offset_left

            if first_outline.child_id2 is not None:
                selected_outline = self.outlines[
                    self.outlines.outline_id == first_outline.child_id2
                ].iloc[0]
                self.selected_outlines.append(selected_outline)
                offt2 = selected_outline.offset_top
                offl2 = selected_outline.offset_left
                offt = (offt + offt2) / 2
                offl = (offl + offl2) / 2

            im2 = self.image_loader.load_frame(fidx, 0)

        elif first_outline.frame_idx + 1 < self.image_loader.num_frames:
            fidx = first_outline.frame_idx + 1
            offt = first_outline.offset_top
            offl = first_outline.offset_left
            im2 = self.image_loader.load_frame(fidx, 0)
        else:
            fidx = first_outline.frame_idx + 1
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

        for outline in database.getOutlinesByFrameIdx(fidx, self.experiment_data.experiment_id):
            c = np.load(outline.coords_path) + np.array([
                outline.offset_left, outline.offset_top
            ])
            kws = dict(
                edgecolor="k",
                lw=1,
                alpha=0.9,
                picker=True,
            )
            if (len(self.selected_outlines) > 0 and
                outline.outline_id in [x.outline_id
                                       for x in self.selected_outlines]):
                p = matplotlib.patches.Polygon(
                    np.array([c[:, 1], c[:, 0]]).T,
                    facecolor="g",
                    **kws
                )
                p.selected = True
            else:
                p = matplotlib.patches.Polygon(
                    np.array([c[:, 1], c[:, 0]]).T,
                    facecolor="y",
                    **kws
                )
                p.selected = False

            p.outline_id = outline.outline_id
            self.plot.axes[2].add_patch(p)

        if first_outline.parent_id:
            prev_outline = database.getOutlineById(first_outline.parent_id)
            im3 = self.image_loader.load_frame(prev_outline.frame_idx, 0)
            self.plot.axes[0].imshow(im3, cmap="gray")
            self.plot.axes[0].set_xlim([
                prev_outline.offset_top,
                prev_outline.offset_top + (self.region_height * 2),
            ])
            self.plot.axes[0].set_ylim([
                prev_outline.offset_left + (self.region_width * 2),
                prev_outline.offset_left,
            ])
            c = np.load(prev_outline.coords_path) + np.array([
                prev_outline.offset_left, prev_outline.offset_top
            ])
            p = matplotlib.patches.Polygon(
                np.array([c[:, 1], c[:, 0]]).T,
                edgecolor="y",
                linestyle="--",
                fill=False,
                lw=1
            )
            self.plot.axes[0].add_patch(p)

        elif first_outline.frame_idx - 1 >= 0:
            im3 = self.image_loader.load_frame(first_outline.frame_idx - 1, 0)
            self.plot.axes[0].imshow(im3, cmap="gray")
            self.plot.axes[0].set_xlim([
                first_outline.offset_top,
                first_outline.offset_top + (self.region_height * 2),
            ])
            self.plot.axes[0].set_ylim([
                first_outline.offset_left + (self.region_width * 2),
                first_outline.offset_left,
            ])

        else:
            prev_outline = None
            im3 = np.zeros((self.region_width * 2, self.region_height * 2))
            self.plot.axes[0].imshow(im3, cmap="gray")


        status_message = "Defining cell lineage {0}: frame {1} ".format(
            self.lineage[0].cell_id,
            first_outline.frame_idx + 1,
        )
        if len(self.selected_outlines) == 1:
            status_message += "[Press ENTER to move to the next frame]"
        elif len(self.selected_outlines) == 2:
            status_message += "[Press ENTER to denote division]"

        self.status_bar.showMessage(status_message)
        self.plot.draw()

    def key_press_event(self, evt):
        if evt.key == "enter":
            if (hasattr(self, "lineage") and len(self.lineage) == 0) or not hasattr(self, "selected_outlines"):
                return

            if len(self.selected_outlines) == 1:
                # next frame
                self.lineage.append(self.selected_outlines[0])
                self.display_frame()
            elif len(self.selected_outlines) > 2:
                print("> 2")
                self.status_bar.showMessage(
                    "!!! Too many outlines selected !!!"
                )
            else:
                write_success = self.write_lineage()
                if not write_success:
                    return

                self.create_layout()
                for outline in self.selected_outlines:
                    self.assignment_queue.append(outline.outline_id)

                if self.assignment_queue:
                    self.assign_lineage(outline_id=self.assignment_queue.pop(0))
                else:
                    self._clear_assignment_plot()
                    self.lineage = []
                    self.selected_outlines = []
                    self.plot.draw()

                self.plot.setFocus()

    def pick_event(self, evt):
        if not evt.artist.selected:
            evt.artist.selected = True
            evt.artist.set_facecolor("g")
            selected_outline = self.outlines[
                self.outlines.outline_id == evt.artist.outline_id
            ].iloc[0]
            self.selected_outlines.append(selected_outline)
        else:
            for i, x in enumerate(self.selected_outlines):
                if x.outline_id == evt.artist.outline_id:
                    self.selected_outlines.pop(i)
                    break

            evt.artist.selected = False
            evt.artist.set_facecolor("y")

        status_message = "Defining cell lineage {0}: frame {1} ".format(
            self.lineage[0].cell_id,
            self.lineage[-1].frame_idx + 1,
        )
        if len(self.selected_outlines) == 1:
            status_message += "[Press ENTER to move to the next frame]"
        elif len(self.selected_outlines) == 2:
            status_message += "[Press ENTER to denote division]"
        elif len(self.selected_outlines) > 2:
            status_message += "[!!! Too many outlines selected !!!]"

        self.status_bar.showMessage(status_message)
        self.plot.draw()

    def write_lineage(self):
        self.status_bar.showMessage("Writing lineage, please wait...")
        if len(self.selected_outlines) == 0:
            death_confirm = QtWidgets.QMessageBox().question(
                self.window,
                "Confirm end of cell lineage",
                "Are you sure this lineage ends here?"
            )
            if death_confirm == QtWidgets.QMessageBox.No:
                status_message = "Defining cell lineage {0}: frame {1} ".format(
                    self.lineage[0].cell_id,
                    self.lineage[-1].frame_idx + 1,
                )
                self.status_bar.showMessage(status_message)
                return False

        elif len(self.selected_outlines) != 2:
            print("Trying to write with {0} selected outlines".format(
                len(self.selected_outlines)
            ))
            self.status_bar.showMessage("Something went wrong, more than 2 outlines were selected")
            return False

        self.temp_window = QtWidgets.QDialog(self.window)
        self.temp_window.setGeometry(0, 0, self.max_width_px * 0.5, self.max_height_px * 0.9)
        self.temp_window.setWindowTitle("Please wait patiently")
        temp_layout = QtWidgets.QVBoxLayout()
        temp_layout.addWidget(QtWidgets.QLabel("Writing lineage, please wait..."))
        self.temp_window.setLayout(temp_layout)
        self.temp_window.show()

        self.status_bar.showMessage("Removing cell_id from extra outlines")
        cell_id = self.lineage[0].cell_id
        for outline in self.lineage:
            database.updateOutlineById(
                outline.outline_id,
                cell_id=None,
            )

        existing_outlines = database.getOutlinesByCellId(cell_id)
        if len(existing_outlines) > 0:
            new_cell_id = str(uuid.uuid4())

        for outline in existing_outlines:
            database.updateOutlineById(
                outline.outline_id,
                cell_id=new_cell_id,
            )

        self.status_bar.showMessage("Assigning parents and children")
        parent_replacements = []
        child_replacements = []
        for i, outline in enumerate(self.lineage):
            if i == 0:
                parent_id = self.lineage[0].parent_id
            else:
                parent_id = self.lineage[i - 1].outline_id

            if i < len(self.lineage) - 1:
                child_id1 = self.lineage[i + 1].outline_id
                child_id2 = ""
            else:
                if len(self.selected_outlines) == 0:
                    child_id1 = ""
                    child_id2 = ""
                else:
                    child_id1 = self.selected_outlines[0].outline_id
                    child_id2 = self.selected_outlines[1].outline_id
                    database.updateOutlineById(
                        child_id1,
                        parent_id=outline.outline_id,
                    )
                    database.updateOutlineById(
                        child_id2,
                        parent_id=outline.outline_id,
                    )

            database.updateOutlineById(
                outline.outline_id,
                cell_id=cell_id,
                parent_id=parent_id,
                child_id1=child_id1,
                child_id2=child_id2,
            )

        # add cell to cells table
        wildtype = False
        database.insertCell(
            cell_id,
            self.lineage[0].experiment_id,
            int(self.lineage[0].frame_idx),
            int(self.lineage[-1].frame_idx),
            self.lineage[0].parent_id is not None,
            self.lineage[-1].child_id2 is not None,
            wildtype,
            self.lineage[0].outline_id,
            self.lineage[-1].outline_id,
        )
        self.status_bar.showMessage("Finished writing lineage {0}".format(cell_id))
        self.temp_window.close()
        self.temp_window.deleteLater()
        self.temp_window = None
        self.plot.setFocus()
        return True
