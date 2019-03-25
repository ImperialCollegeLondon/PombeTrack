#!/usr/bin/env python3

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.widgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.figure
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
from . import movie_generator

sns.set_context("talk")
sns.set_style("white")

class Toolbar(NavigationToolbar):
    def __init__(self, figure_canvas, parent=None):
        self.toolitems = [
            ("Home", "Home", "home_large", "home_event"),
            # (None, None, None, None),
            ("Pan", "Pan", "move_large", "pan"),
            ("Zoom", "Zoom", "zoom_to_rect_large", "zoom"),
            (None, None, None, None),
            ("Previous", "Previous frame", "previous", "previous_frame"),
            ("Next", "Next frame", "next", "next_frame"),
            (None, None, None, None),
            ("Accept", "Accept assignment", "accept", "accept"),
            ("Cancel", "Cancel assignment", "cancel", "cancel"),
        ]
        NavigationToolbar.__init__(self, figure_canvas, parent=parent)

    def _icon(self, name):
        path = os.path.join("resources", name)
        if not os.path.exists(path):
            path = os.path.join(self.basedir, name)

        pm = QtGui.QPixmap(path)
        if hasattr(pm, "setDevicePixelRatio"):
            pm.setDevicePixelRatio(self.canvas._dpi_ratio)

        return QtGui.QIcon(pm)

    def home_event(self, *args, **kwargs):
        for ax, lims in zip(self.canvas.axes, self.canvas.offsets):
            ax.set_xlim(lims[3], lims[1])
            ax.set_ylim(lims[0], lims[2])
        self.canvas.draw()

    def previous_frame(self, *args, **kwargs):
        self.parent.trigger_previous_frame()

    def next_frame(self, *args, **kwargs):
        self.parent.trigger_next_frame()

    def accept(self, *args, **kwargs):
        self.parent.trigger_accept_all()

    def cancel(self, *args, **kwargs):
        self.parent.trigger_cancel()


class Plotter(FigureCanvas):
    def __init__(self, parent_window, width, height, dpi, experiment_data, image_loader, subplots=3):
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
            self.offsets.append((0, 0, 0, 0))  # top, right, bottom, left (like CSS)

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
        self.window.setGeometry(0, 60, self.max_width_px * 0.5, self.max_height_px * 0.9)
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

        self.window.trigger_previous_frame = self.previous_frame
        self.window.trigger_next_frame = self.next_frame
        self.window.trigger_accept_all = self.accept_all
        self.window.trigger_cancel = self.cancel_assignment

        self.window.toolbar = Toolbar(self.plot, self.window)
        tool_layout = QtWidgets.QVBoxLayout()
        tool_layout.addWidget(self.window.toolbar)

        self.main_layout.addLayout(tool_layout)
        self.main_layout.addWidget(self.plot)

        self.status_bar = QtWidgets.QStatusBar()
        self.main_layout.addWidget(self.status_bar)

        self.window.setLayout(self.main_layout)
        self.window.show()

        self.plot.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.window.toolbar.hide()
        self.plot.hide()

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
            cell_box = QtWidgets.QWidget()
            # create proper lineage
            cell_outlines = self.outlines[
                self.outlines.cell_id == cell_id
            ].sort_values("frame_idx")

            plot_layout = QtWidgets.QHBoxLayout()

            spacer_l = QtWidgets.QHBoxLayout()
            spacer = QtWidgets.QWidget()
            if verified_cell:
                spacer.setStyleSheet("background-color:green")
            else:
                spacer.setStyleSheet("background-color:red")
            spacer.setMaximumWidth(5)
            spacer_l.addWidget(spacer)
            plot_layout.addLayout(spacer_l)

            # only take the first and last frames
            cell_plots = []
            for outline_num in [0, -1]:
                outline = cell_outlines.iloc[outline_num]
                width = self.max_width_px * 0.1
                cell_plot = Plotter(
                    self.window,
                    width=self._px_to_in(width),
                    height=self._px_to_in(width),
                    dpi=self.screen_dpi,
                    experiment_data=self.experiment_data,
                    image_loader=self.image_loader,
                    subplots=1
                )
                cell_plot.setMinimumWidth(width)
                cell_plot.setMaximumWidth(width)
                cell_plot.setMinimumHeight(width)
                roi = self.image_loader.load_frame(outline.frame_idx, 0)[
                    outline.offset_left:outline.offset_left + (self.region_width * 2),
                    outline.offset_top:outline.offset_top + (self.region_height * 2),
                ]
                cell_plot.axes[0].imshow(roi, cmap="gray")
                cell_plot.axes[0].set_title("F{0}".format(outline.frame_idx + 1))
                cell_plot.current_channel = 0
                cell_plots.append(cell_plot)
                plot_layout.addWidget(cell_plot)

            control_layout = QtWidgets.QVBoxLayout()
            control_layout.setAlignment(QtCore.Qt.AlignTop)
            info_layout = QtWidgets.QHBoxLayout()
            if verified_cell:
                pixmap = QtGui.QPixmap("resources/tick.png")
                if verified_cell.is_wildtype:
                    wildtype = True
                    wildtype_btn = QtWidgets.QPushButton("Unset wildtype")
                    desc_str = "Wildtype cell {0}".format(cell_id)
                else:
                    wildtype = False
                    wildtype_btn = QtWidgets.QPushButton("Set wildtype")
                    desc_str = "Cell {0}".format(cell_id)
            else:
                pixmap = QtGui.QPixmap("resources/cross.png")
                wildtype = False
                wildtype_btn = None
                desc_str = "Cell {0}".format(cell_id)

            verification_sign = pixmap.scaledToWidth(20)
            verification_label = QtWidgets.QLabel()
            verification_label.setPixmap(verification_sign)
            info_layout.addWidget(verification_label)

            desc_label = QtWidgets.QLabel(desc_str)
            info_layout.addWidget(desc_label)
            control_layout.addLayout(info_layout)

            row1 = QtWidgets.QHBoxLayout()
            verify_btn = QtWidgets.QPushButton("Assign Cell Lineage")
            verify_btn._cell_id = cell_id
            verify_btn.clicked.connect(self.assign_lineage)
            row1.addWidget(verify_btn)
            export_btn = QtWidgets.QPushButton("Export movie")
            export_btn._cell_id = cell_id
            export_btn.clicked.connect(self.export_movie)
            row1.addWidget(export_btn)
            control_layout.addLayout(row1)

            if wildtype_btn:
                row2 = QtWidgets.QHBoxLayout()
                switch_channel = QtWidgets.QPushButton("Change channel")
                switch_channel._cell_id = cell_id
                switch_channel._cell_plots = cell_plots
                switch_channel.clicked.connect(self.switch_channel)
                row2.addWidget(switch_channel)
                wildtype_btn._cell_id = cell_id
                wildtype_btn.clicked.connect(self.toggle_wildtype)
                row2.addWidget(wildtype_btn)
                control_layout.addLayout(row2)

            details_label = QtWidgets.QLabel(
                "{0} cell with {1} frames (F{2} - F{3}), ending in {4}".format(
                    verified_cell and "Verified" or "Unverified",
                    cell_outlines.iloc[-1].frame_idx + 1 - cell_outlines.iloc[0].frame_idx,
                    cell_outlines.iloc[0].frame_idx + 1,
                    cell_outlines.iloc[-1].frame_idx + 1,
                    cell_outlines.iloc[-1].child_id1 and "division" or "loss",
                )
            )
            control_layout.addWidget(details_label)

            plot_layout.addLayout(control_layout)
            cell_box.setLayout(plot_layout)
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

    def exit_assignment(self):
        self.create_layout()
        self.lineage_scroll_area.show()
        self._clear_assignment_plot()
        self.window.toolbar.hide()
        self.plot.hide()
        self.lineage = []
        self.selected_outlines = []

    def get_child_ids(self, cell_id):
        this_cell = database.getCellById(cell_id)
        cell_ids = []
        if this_cell.child_cell_id1:
            cell_ids.append(this_cell.child_cell_id1)
            cell_ids.extend(self.get_child_ids(this_cell.child_cell_id1))
            cell_ids.append(this_cell.child_cell_id2)
            cell_ids.extend(self.get_child_ids(this_cell.child_cell_id2))
        return cell_ids

    def toggle_wildtype(self, click=False):
        cell_id = self.window.sender()._cell_id
        this_cell = database.getCellById(cell_id)
        if this_cell:
            setting = not this_cell.is_wildtype
            database.updateCellById(cell_id, is_wildtype=setting)
            while this_cell.parent_cell_id:
                this_cell = database.getCellById(this_cell.parent_cell_id)
                database.updateCellById(this_cell.cell_id, is_wildtype=setting)

            for child_id in self.get_child_ids(cell_id):
                database.updateCellById(child_id, is_wildtype=setting)
            self.create_layout()

    def switch_channel(self, click=False):
        cell_outlines = self.outlines[
            self.outlines.cell_id == self.window.sender()._cell_id
        ].sort_values("frame_idx")
        cell_plots = self.window.sender()._cell_plots
        for plot, outline in zip(
            cell_plots,
            [cell_outlines.iloc[0], cell_outlines.iloc[-1]]
        ):
            plot.current_channel += 1
            if plot.current_channel == self.image_loader.num_channels:
                plot.current_channel = 0

            roi = self.image_loader.load_frame(outline.frame_idx, plot.current_channel)[
                outline.offset_left:outline.offset_left + (self.region_width * 2),
                outline.offset_top:outline.offset_top + (self.region_height * 2),
            ]

            if plot.current_channel != 0:
                plot.axes[0].imshow(roi, cmap="binary")
                plot.axes[0].set_title("F{0} (C{1})".format(
                    outline.frame_idx + 1,
                    plot.current_channel + 1,
                ))
            else:
                plot.axes[0].imshow(roi, cmap="gray")
                plot.axes[0].set_title("F{0}".format(
                    outline.frame_idx + 1,
                ))

            plot.draw()

    def export_movie(self, click=False):
        cell_id = self.window.sender()._cell_id
        settings = QtWidgets.QDialog(self.window)
        settings.setModal(True)
        settings.setWindowTitle("Cell lineage exporter")
        layout = QtWidgets.QVBoxLayout()

        checkbox_layout = QtWidgets.QGridLayout()
        descendants = QtWidgets.QCheckBox()
        descendants.setCheckState(QtCore.Qt.Checked)
        checkbox_layout.addWidget(descendants, 0, 0)
        checkbox_layout.addWidget(QtWidgets.QLabel("Include descendants"), 0, 1)

        outlines = QtWidgets.QCheckBox()
        outlines.setCheckState(QtCore.Qt.Checked)
        checkbox_layout.addWidget(outlines, 1, 0)
        checkbox_layout.addWidget(QtWidgets.QLabel("Include outlines"), 1, 1)

        frames = QtWidgets.QCheckBox()
        frames.setCheckState(QtCore.Qt.Checked)
        checkbox_layout.addWidget(frames, 2, 0)
        checkbox_layout.addWidget(QtWidgets.QLabel("Include frame numbers"), 2, 1)

        scalebar = QtWidgets.QCheckBox()
        scalebar.setCheckState(QtCore.Qt.Checked)
        checkbox_layout.addWidget(scalebar, 3, 0)
        checkbox_layout.addWidget(QtWidgets.QLabel("Include scale bar"), 3, 1)

        channels = QtWidgets.QCheckBox()
        checkbox_layout.addWidget(channels, 4, 0)
        checkbox_layout.addWidget(QtWidgets.QLabel("Include all channels"), 4, 1)

        def generate_movie():
            settings = {
                "experiment_id": self.experiment_data.experiment_id,
                "cell_id": cell_id,
                "descendants": descendants.checkState() == 2,
                "outlines": outlines.checkState() == 2,
                "frames": frames.checkState() == 2,
                "scale": scalebar.checkState() == 2,
                "channels": channels.checkState() == 2,
                "image_loader": self.image_loader,
                "screen_dpi": self.screen_dpi,
            }
            m = movie_generator.MovieMaker(**settings)

        button_layout = QtWidgets.QHBoxLayout()
        submit = QtWidgets.QPushButton("Generate")
        submit.clicked.connect(generate_movie)
        button_layout.addWidget(submit)

        layout.addLayout(checkbox_layout)
        layout.addLayout(button_layout)
        settings.setLayout(layout)
        settings.show()

    def assign_lineage(self, click=False, outline_id=None):
        self.lineage_scroll_area.hide()
        self.get_outlines()
        if not outline_id:
            potential_outlines = self.outlines[self.outlines.cell_id == self.window.sender()._cell_id]
        else:
            potential_outlines = self.outlines[self.outlines.outline_id == outline_id]

        self.lineage = [potential_outlines[
            potential_outlines.frame_idx == potential_outlines.frame_idx.min()
        ].iloc[0]]
        self.selected_outlines = []
        self.display_frame()
        self.window.toolbar.show()
        self.plot.show()
        self.plot.setFocus()

    def display_frame(self):
        self.selected_outlines = []
        self._clear_assignment_plot()
        first_outline = self.lineage[-1]

        self.plot.axes[0].set_title("Frame {0}".format(first_outline.frame_idx))
        self.plot.axes[1].set_title("Frame {0}".format(first_outline.frame_idx + 1))
        self.plot.axes[2].set_title("Frame {0}".format(first_outline.frame_idx + 2))

        im1 = self.image_loader.load_frame(first_outline.frame_idx, 0)
        offset_left, offset_right = (first_outline.offset_top,
                                     first_outline.offset_top + (self.region_height * 2))
        offset_top, offset_bottom = (first_outline.offset_left + (self.region_width * 2),
                                     first_outline.offset_left)
        self.plot.offsets[1] = (offset_top, offset_right, offset_bottom, offset_left)
        self.plot.axes[1].imshow(im1, cmap="gray")
        self.plot.axes[1].set_xlim([offset_left, offset_right])
        self.plot.axes[1].set_ylim([offset_top, offset_bottom])
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

            if first_outline.child_id2 is not None and first_outline.child_id1 != first_outline.child_id2:
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
        self.plot.offsets[2] = (
            offl + (self.region_width * 2),
            offt + (self.region_height * 2),
            offl,
            offt
        )
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
            self.plot.offsets[0] = (
                prev_outline.offset_left + (self.region_width * 2),
                prev_outline.offset_top + (self.region_height * 2),
                prev_outline.offset_left,
                prev_outline.offset_top,
            )
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
            self.plot.offsets[0] = (
                first_outline.offset_left + (self.region_width * 2),
                first_outline.offset_top + (self.region_height * 2),
                first_outline.offset_left,
                first_outline.offset_top,
            )
        else:
            prev_outline = None
            im3 = np.zeros((self.region_width * 2, self.region_height * 2))
            self.plot.axes[0].imshow(im3, cmap="gray")

        self.plot.axes[0].set_xlim([
            self.plot.offsets[0][3],
            self.plot.offsets[0][1]
        ])
        self.plot.axes[0].set_ylim([
            self.plot.offsets[0][0],
            self.plot.offsets[0][2]
        ])


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
        if evt.key == "enter" or evt.key == "right":
            self.next_frame()
        elif evt.key == "left":
            self.previous_frame()

    def accept_all(self):
        confirm_accept = QtWidgets.QMessageBox().question(
            self.window,
            "Confirm automatic assignment",
            "Please confirm that you wish to accept this cell lineage",
            QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes,
            QtWidgets.QMessageBox.No,
        )
        if confirm_accept == QtWidgets.QMessageBox.No:
            return

        outline = self.lineage[-1]
        while outline.child_id1 and not outline.child_id2:
            outline = database.getOutlineById(outline.child_id1)
            self.lineage.append(outline)

        if outline.child_id1:
            selection = [
                database.getOutlineById(outline.child_id1),
                database.getOutlineById(outline.child_id2),
            ]
            self.assignment_queue.append(outline.child_id1)
            self.assignment_queue.append(outline.child_id2)
        else:
            selection = []

        self.write_lineage(selected_outlines=selection)

        if self.assignment_queue:
            next_id = self.assignment_queue.pop()
            self.assign_lineage(outline_id=next_id)
        else:
            self.exit_assignment()

        self.plot.setFocus()

    def cancel_assignment(self):
        confirm_cancel = QtWidgets.QMessageBox().question(
            self.window,
            "Confirm cancellation",
            "Please confirm that you wish to cancel assigning this cell lineage",
            QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes,
            QtWidgets.QMessageBox.No,
        )
        if confirm_cancel == QtWidgets.QMessageBox.No:
            return

        self.exit_assignment()

    def previous_frame(self):
        if (hasattr(self, "lineage") and len(self.lineage) <= 1):
            return

        self.lineage.pop()
        self.display_frame()
        self.plot.setFocus()

    def next_frame(self):
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

            for outline in self.selected_outlines:
                self.assignment_queue.append(outline.outline_id)

            if self.assignment_queue:
                next_id = self.assignment_queue.pop()
                self.assign_lineage(outline_id=next_id)
            else:
                self.exit_assignment()

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

    def write_lineage(self, selected_outlines=None):
        self.status_bar.showMessage("Writing lineage, please wait...")
        if selected_outlines is None:
            if len(self.selected_outlines) == 0:
                death_confirm = QtWidgets.QMessageBox().question(
                    self.window,
                    "Confirm end of cell lineage",
                    "Are you sure this lineage ends here?",
                    QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes,
                    QtWidgets.QMessageBox.No,
                )
                if death_confirm == QtWidgets.QMessageBox.No:
                    status_message = "Defining cell lineage {0}: frame {1} ".format(
                        self.lineage[0].cell_id,
                        self.lineage[-1].frame_idx + 1,
                    )
                    self.status_bar.showMessage(status_message)
                    return False

            elif len(self.selected_outlines) == 2:
                div_confirm = QtWidgets.QMessageBox().question(
                    self.window,
                    "Confirm cell division",
                    "Are you sure this cell divides?\nThere is no going back if you say yes...",
                    QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes,
                    QtWidgets.QMessageBox.No,
                )
                if div_confirm == QtWidgets.QMessageBox.No:
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

        else:
            self.selected_outlines = selected_outlines

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

        wildtype = False
        parent_cell_id = None
        child_cell_id1 = None
        child_cell_id2 = None

        if self.lineage[0].parent_id:
            parent_cell_id = database.getOutlineById(
                self.lineage[0].parent_id
            ).cell_id
            parent_cell = database.getCellById(parent_cell_id)
            wildtype = parent_cell and parent_cell.is_wildtype or False

        if self.lineage[-1].child_id2:
            child_cell_id1 = database.getOutlineById(
                self.lineage[-1].child_id1
            ).cell_id
            child_cell_id2 = database.getOutlineById(
                self.lineage[-1].child_id2
            ).cell_id

        existing = database.getCellById(cell_id)
        if existing:
            database.deleteCellById(cell_id)

        kwargs = {
            "cell_id": cell_id,
            "experiment_id": self.lineage[0].experiment_id,
            "start_frame_idx": int(self.lineage[0].frame_idx),
            "end_frame_idx": int(self.lineage[-1].frame_idx),
            "birth_observed": self.lineage[0].parent_id is not None,
            "division_observed": self.lineage[-1].child_id2 is not None,
            "is_wildtype": wildtype,
            "first_outline_id": self.lineage[0].outline_id,
            "last_outline_id": self.lineage[-1].outline_id,
            "parent_cell_id": parent_cell_id,
            "child_cell_id1": child_cell_id1,
            "child_cell_id2": child_cell_id2,
        }
        database.insertCell(**kwargs)
        self.status_bar.showMessage("Finished writing lineage {0}".format(cell_id))
        self.temp_window.close()
        self.temp_window.deleteLater()
        self.temp_window = None
        self.plot.setFocus()
        return True
