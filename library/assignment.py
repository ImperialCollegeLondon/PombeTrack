#!/usr/bin/env python3

""" Interface for assembling cells and assigning lineages.  """

import os
import uuid

import matplotlib
import matplotlib.figure
import matplotlib.widgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

import numpy as np
import pandas as pd
import PyQt5.QtWidgets as QtWidgets
import PyQt5.QtGui as QtGui
import PyQt5.QtCore as QtCore
import seaborn as sns

from . import database
from . import movie_generator

matplotlib.use('Qt5Agg')

sns.set_context("talk")
sns.set_style("white")

class Toolbar(NavigationToolbar):
    """Custom toolbar for controlling cell assembly and lineage assignment."""
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

        pixmap = QtGui.QPixmap(path)
        if hasattr(pixmap, "setDevicePixelRatio"):
            pixmap.setDevicePixelRatio(self.canvas._dpi_ratio)

        return QtGui.QIcon(pixmap)

    def home_event(self):
        """Set axis limits to the default limits.

        Triggered when the home button is pressed in the toolbar."""
        for axis, lims in zip(self.canvas.axes, self.canvas.offsets):
            axis.set_xlim(lims[0], lims[2])
            axis.set_ylim(lims[1], lims[3])

        self.canvas.draw()

    def previous_frame(self):
        """Change frame to previous frame on toolbar click.

        Actually calls `Assigner.previous_frame` function."""
        self.parent.trigger_previous_frame()

    def next_frame(self):
        """Change frame to next frame on toolbar click.

        Actually calls `Assigner.next_frame` function."""
        self.parent.trigger_next_frame()

    def accept(self):
        """Accepts the current lineage on toolbar click.

        Includes all relatives of the current cell.
        Actually calls `Assigner.accept_all` function."""
        self.parent.trigger_accept_all()

    def cancel(self):
        """Cancels assignment of lineage, discarding any changes.

        Actually calls `Assigner.cancel_assignment` function."""
        self.parent.trigger_cancel()


class Plotter(FigureCanvas):
    """Interface for assigning outlines to cells and constructing cell lineages.

    Interface creates three subplots representing consecutive frames of the
    timelapse with a single outline focus.

    Actual interactivity functionality is encoded in the Analysis class, but is
    also described here for convenience.

    ---- Left subplot ----
    This shows the frame prior to the current outline, with the assigned
    "parent" outline (if division has just occurred this will be a true parent,
    otherwise it will display the _same_ outline just in the previous frame)
    outlined with a yellow dotted line.

    If there is no prior frame, and black image will be displayed.
    If no parent outline is known, no cell will be outlined.

    --- Middle subplot ---
    This shows the outline in its current frame, outlined in red.

    --- Right subplot ---
    This shows the subsequent frame compared to the current frame.
    All known cells in the vicinity of the current outline are displayed as
    shaded polygons.

    Outlines that are connected to the current outline are coloured in green,
    else in yellow.
    If a single outline is coloured in green, a "cell growth" event will be
    assigned, i.e. no division occurs.
    If two outlines are coloured in green, a division event will be assigned.
    If no outlines are coloured in green, a loss/death event will be assigned.

    The initial state of the right-hand subplot denotes the currently known
    state of the data.

    These assignment options can be altered simply by clicking on the outlines
    to toggle their colour between green and yellow.

    Assignments can be accepted by pressing the ENTER key, the RIGHT arrow key,
    or the > arrow button in the toolbar, advancing the timelapse by one frame.
    Past assignments can be be altered by pressing the LEFT arrow key, or the <
    arrow button in the toolbar, which moves the timelapse back one frame,
    permitting reassignment.
    This alteration method however only works within a single cell (i.e. not
    across division events).

    The tick toolbar button results in the automatic acceptance of the cell
    assignments, including all "growth" events and all division events for the
    cell and all its ancestors.

    The cross button results in the discarding of all assignments and the
    closing of the assignment interface."""
    def __init__(self, parent_window, fig_dimensions, experiment_data,
                 image_loader):
        self.current_channel = None
        width, height, dpi, subplots = fig_dimensions
        fig = matplotlib.figure.Figure(figsize=(width, height), dpi=dpi)
        self.axes = []
        self.offsets = []
        for sp_idx in range(subplots):
            axis = fig.add_subplot(1, subplots, sp_idx + 1)
            # axis.axis("off")
            axis.set_xticks([], [])
            axis.set_yticks([], [])
            axis.set_aspect("equal")
            axis.autoscale("off")
            self.axes.append(axis)
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

    def set_channel(self, channel):
        """Assign channel number to internal variable."""
        self.current_channel = channel


class Assigner:
    """Provides the assigner interface

    The interface will display all potention "cells" and allow them to be
    defined, and to be collected into lineages.

    Each outline has the parameters `parent_id`, `child_id1`, and `child_id2`.
    `parent_id` is not necessarily a "parent", since we are dealing with
    outlines in single frames at a time.
    Instead, the `parent_id` is defined as the `outline_id` of the outline in
    the previous frame which is associated with each outline.
    If division has occurred between the two frames, then the `parent_id` will
    reflect true divison, otherwise it will reflect growth.

    Similarly, `child_id1` does not imply that division occurs, it is populated
    with the `outline_id` of the associated outline in the subsequent frame.
    If `child_id1` is not defined, the lineage is terminated at that point.
    If `child_id2` is not defined (but `child_id1` is), a "growth event" has
    occurred.
    If both parameters are defined, then true division has occurred.

    Outlines from birth until division are defined as a "cell", which are
    comprised of collections of outlines over multiple frames.

    This class provides an interface to assemble/verify which outlines belong
    to the same cell, and a method to define division events.

    The interface is arranged as a series of "cards" vertically.
    Each card represents a single cell.
    From left to right, each card contains:
    - A red/green bar signalling whether the cell has been verified.
    - Two images with the prospective cell outlined at birth (or the first
      frame in which it is observed) and the cell outlined at division or
      loss.
    - A control/information panel:
        - A cross or tick image signalling the cell verification state
        - The `cell_id` and whether it has been set to be wildtype
        - Information about the cell e.g.:
            "Unverified cell with 4 frames (F2 - F5), ending in loss"
            "Verified cell with 40 frames (F10 - F49), ending in division"
          I believe these are self-explanatory.
        - The buttons <Assign Cell Lineage>, <Export movie>, <Change channel>,
          <Set/Unset wildtype>; see below.

    ___Control buttons___
    __Assign Cell Lineage__
    This button launches the assignment interface in a dialog window.
    See the `Plotter` class for information about the interface.

    __Export movie__
    This button allows movies to be generated.
    A popup interface allows parameters to be selected:
        - <Include descendants>: a movie rooted from the selected cell will be
            generated, including all its descendants.
        - <Include outlines>: cell outlines will be drawn in the movie.
        - <Include frame numbers>: frame numbers will be annotated.
        - <Include scale bar>: a 10 µm scale bar will be annotated.
        - <Include all channels>: if selected, the movie will display up to
            three panels with brightfield, and two fluorescent channels
            displayed. If not selected, only brightfield will be exported.

    Data will be exported as a series of PNG files for each frame, into the
    directory data/movies/<experiment_id>/<cell_id>.

    True movies can be subsequently generated from these files with third-party
    software, e.g. ImageJ.

    __Change channel__
    This button allows cycling through each channel in the left-hand images,
    for example to determine quickly whether a cell is wildtype.

    __Set/Unset wildtype__
    Toggles the wildtype status of a cell (and all its descendants).
    This information may be useful for subsequent analysis.
    This button is only displayed when a cell has been verified (green bar and
    tick mark).
    """
    # pylint: disable=too-many-instance-attributes
    # not sure how to get fewer attributes for functionality really
    # pylint: disable=no-member
    # apparently database.py can't be introspected properly
    def __init__(self, experiment_data, image_loader, screen_dimensions):
        self.experiment_data = experiment_data
        self.image_loader = image_loader
        self.region_halfwidth = 75
        self.region_halfheight = 75

        if not database.checkTable("cells"):
            database.createCellsTable()

        self.max_width_px = screen_dimensions[0]
        self.max_height_px = screen_dimensions[1]
        self.screen_dpi = screen_dimensions[2]

        self.outlines = pd.DataFrame()
        self.parent_window = None
        self.window = None
        self.temp_window = None
        self.main_layout = None
        self.plot = None
        self.status_bar = None
        self.lineage_scroll_area = None
        self.lineage = []
        self.selected_outlines = []
        self.assignment_queue = []

    def _px_to_in(self, num_pixels):
        return num_pixels / self.screen_dpi

    def get_outlines(self):
        outlines = database.getOutlinesByExperimentId(self.experiment_data.experiment_id)
        self.outlines = pd.DataFrame(outlines)

    def start_assigning(self, parent_window):
        self.parent_window = parent_window
        self.window = QtWidgets.QDialog(self.parent_window)
        self.window.setModal(True)
        self.window.setGeometry(0, 60, self.max_width_px * 0.5, self.max_height_px * 0.9)
        self.window.setWindowTitle("Assign/Verify cell lineages")

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
            fig_dimensions=(
                self._px_to_in(self.max_width_px * 0.5),
                self._px_to_in((self.max_width_px * 0.5) / 3),
                self.screen_dpi,
                3,
            ),
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

    def get_offsets(self, centre):
        offset_left = int(round(centre[0] - self.region_halfwidth))
        offset_top = int(round(centre[1] - self.region_halfheight))
        img = self.image_loader.load_frame(0)
        if offset_left < 0:
            offset_left = 0
        elif offset_left >= img.shape[0] - (self.region_halfwidth * 2):
            offset_left = img.shape[0] - (self.region_halfwidth * 2)

        if offset_top < 0:
            offset_top = 0
        elif offset_top >= img.shape[1] - (self.region_halfheight * 2):
            offset_top = img.shape[1] - (self.region_halfheight * 2)
        del img

        return offset_left, offset_top

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
        for cell_id, verified_cell in unique_cells:
            lineage_layout.addWidget(
                self.create_cell_box(cell_id, verified_cell)
            )

        lineage_widget = QtWidgets.QWidget()
        lineage_widget.setMinimumWidth(self.max_width_px * 0.5 - 50)
        lineage_widget.setLayout(lineage_layout)

        if self.lineage_scroll_area is None:
            self.lineage_scroll_area = QtWidgets.QScrollArea()
            self.lineage_scroll_area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)

        self.lineage_scroll_area.setWidget(lineage_widget)
        self.main_layout.addWidget(self.lineage_scroll_area)

    def create_cell_box(self, cell_id, verified_cell):
        cell_box = QtWidgets.QWidget()
        # create proper lineage
        cell_outlines = self.outlines[
            self.outlines.cell_id == cell_id
        ].sort_values("frame_idx")

        plot_layout = QtWidgets.QHBoxLayout()
        spacer = QtWidgets.QWidget()
        if verified_cell:
            spacer.setStyleSheet("background-color:green")
        else:
            spacer.setStyleSheet("background-color:red")
        spacer.setMaximumWidth(5)
        plot_layout.addWidget(spacer)

        # only take the first and last frames
        cell_plot1 = self.create_cell_plot(cell_outlines.iloc[0])
        cell_plot2 = self.create_cell_plot(cell_outlines.iloc[-1])
        plot_layout.addWidget(cell_plot1)
        plot_layout.addWidget(cell_plot2)

        control_layout = self.create_control_layout(
            cell_outlines,
            cell_id,
            verified_cell,
            [cell_plot1, cell_plot2],
        )
        plot_layout.addLayout(control_layout)
        cell_box.setLayout(plot_layout)
        return cell_box

    def create_cell_plot(self, outline):
        width = self.max_width_px * 0.1
        cell_plot = Plotter(
            self.window,
            fig_dimensions=(
                self._px_to_in(width),
                self._px_to_in(width),
                self.screen_dpi,
                1,
            ),
            experiment_data=self.experiment_data,
            image_loader=self.image_loader,
        )
        cell_plot.setMinimumWidth(width)
        cell_plot.setMaximumWidth(width)
        cell_plot.setMinimumHeight(width)
        centre = outline.centre_y, outline.centre_x
        outline.offset_left, outline.offset_top = self.get_offsets(centre)
        roi = self.image_loader.load_frame(outline.frame_idx, 0)[
            outline.offset_left:outline.offset_left + (self.region_halfwidth * 2),
            outline.offset_top:outline.offset_top + (self.region_halfheight * 2),
        ]
        cell_plot.axes[0].imshow(roi, cmap="gray")
        cell_plot.axes[0].set_title("F{0}".format(outline.frame_idx + 1))
        coords = np.load(outline.coords_path) - np.array(
            [outline.offset_left, outline.offset_top]
        )
        outline_poly = matplotlib.patches.Polygon(
            np.array([coords[:, 1], coords[:, 0]]).T,
            edgecolor="r",
            fill=False,
            lw=1
        )
        cell_plot.axes[0].add_patch(outline_poly)
        cell_plot.set_channel(0)  # current_channel = 0
        return cell_plot

    def create_control_layout(self, cell_outlines, cell_id,
                              verified_cell, cell_plots):
        control_layout = QtWidgets.QVBoxLayout()
        control_layout.setAlignment(QtCore.Qt.AlignTop)
        info_layout = self.create_info_layout(cell_id, verified_cell)
        control_layout.addLayout(info_layout)

        button_layouts = self.create_control_button_layouts(
            cell_id, verified_cell, cell_plots
        )
        for layout in button_layouts:
            control_layout.addLayout(layout)

        control_layout.addWidget(QtWidgets.QLabel(
            "{0} cell with {1} frames (F{2} - F{3}), ending in {4}".format(
                verified_cell and "Verified" or "Unverified",
                cell_outlines.iloc[-1].frame_idx + 1 - cell_outlines.iloc[0].frame_idx,
                cell_outlines.iloc[0].frame_idx + 1,
                cell_outlines.iloc[-1].frame_idx + 1,
                cell_outlines.iloc[-1].child_id1 and "division" or "loss",
            )
        ))
        return control_layout

    @staticmethod
    def create_info_layout(cell_id, verified_cell):
        info_layout = QtWidgets.QHBoxLayout()
        if verified_cell:
            pixmap = QtGui.QPixmap("resources/tick.png").scaledToWidth(20)
            if verified_cell.is_wildtype:
                desc_str = "Wildtype cell {0}".format(cell_id)
            else:
                desc_str = "Cell {0}".format(cell_id)
        else:
            pixmap = QtGui.QPixmap("resources/cross.png").scaledToWidth(20)
            desc_str = "Cell {0}".format(cell_id)

        verification_label = QtWidgets.QLabel()
        verification_label.setPixmap(pixmap)
        info_layout.addWidget(verification_label)

        info_layout.addWidget(QtWidgets.QLabel(desc_str))
        return info_layout

    def create_control_button_layouts(self, cell_id, verified_cell,
                                      cell_plots):
        row1 = QtWidgets.QHBoxLayout()
        verify_btn = QtWidgets.QPushButton("Assign Cell Lineage")
        verify_btn.cell_id = cell_id
        verify_btn.clicked.connect(self.assign_lineage)
        row1.addWidget(verify_btn)
        export_btn = QtWidgets.QPushButton("Export movie")
        export_btn.cell_id = cell_id
        export_btn.clicked.connect(self.export_movie)
        row1.addWidget(export_btn)

        if verified_cell:
            row2 = QtWidgets.QHBoxLayout()
            switch_channel = QtWidgets.QPushButton("Change channel")
            switch_channel.cell_id = cell_id
            switch_channel.cell_plots = cell_plots
            switch_channel.clicked.connect(self.switch_channel)
            row2.addWidget(switch_channel)

            if verified_cell.is_wildtype:
                wildtype_btn = QtWidgets.QPushButton("Unset wildtype")
            else:
                wildtype_btn = QtWidgets.QPushButton("Set wildtype")

            wildtype_btn.cell_id = cell_id
            wildtype_btn.clicked.connect(self.toggle_wildtype)
            row2.addWidget(wildtype_btn)
            return [row1, row2]

        return [row1]


    def _clear_assignment_plot(self):
        for axis in self.plot.axes:
            axis.clear()
            axis.set_xticks([], [])
            axis.set_yticks([], [])
            axis.set_aspect("equal")
            axis.autoscale("off")

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
        if this_cell and this_cell.child_cell_id1:
            cell_ids.append(this_cell.child_cell_id1)
            cell_ids.extend(self.get_child_ids(this_cell.child_cell_id1))
            cell_ids.append(this_cell.child_cell_id2)
            cell_ids.extend(self.get_child_ids(this_cell.child_cell_id2))
        return cell_ids

    def toggle_wildtype(self):
        cell_id = self.window.sender().cell_id
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

    def switch_channel(self):
        cell_outlines = self.outlines[
            self.outlines.cell_id == self.window.sender().cell_id
        ].sort_values("frame_idx")
        cell_plots = self.window.sender().cell_plots
        for plot, outline in zip(
            cell_plots,
            [cell_outlines.iloc[0], cell_outlines.iloc[-1]]
        ):
            plot.set_channel(plot.current_channel + 1)
            if plot.current_channel == self.image_loader.num_channels:
                plot.set_channel(0)

            roi = self.image_loader.load_frame(outline.frame_idx, plot.current_channel)[
                outline.offset_left:outline.offset_left + (self.region_halfwidth * 2),
                outline.offset_top:outline.offset_top + (self.region_halfheight * 2),
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

    def export_movie(self):
        cell_id = self.window.sender().cell_id
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
            movie_generator.MovieMaker(**settings)

        button_layout = QtWidgets.QHBoxLayout()
        submit = QtWidgets.QPushButton("Generate")
        submit.clicked.connect(generate_movie)
        button_layout.addWidget(submit)

        layout.addLayout(checkbox_layout)
        layout.addLayout(button_layout)
        settings.setLayout(layout)
        settings.show()

    def assign_lineage(self, outline_id=None):
        self.lineage_scroll_area.hide()
        self.get_outlines()
        if not outline_id:
            potential_outlines = self.outlines[
                (self.outlines.cell_id == self.window.sender().cell_id)
            ]
        else:
            potential_outlines = self.outlines[
                (self.outlines.outline_id == outline_id)
            ]

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

        self.display_frame_centre(first_outline)
        self.display_frame_right(first_outline)
        self.display_frame_left(first_outline)
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

    def display_frame_centre(self, outline):
        # y0, x1, y1, x0
        img = self.image_loader.load_frame(outline.frame_idx, 0)
        offset_left, offset_top = self.get_offsets((
            outline.centre_y,
            outline.centre_x,
        ))
        self.plot.offsets[1] = (
            offset_top,
            offset_left + self.region_halfwidth * 2,
            offset_top + self.region_halfheight * 2,
            offset_left,
        )
        self.plot.axes[1].imshow(img, cmap="gray")
        self.plot.axes[1].set_xlim([
            self.plot.offsets[1][0],
            self.plot.offsets[1][2],
        ])
        self.plot.axes[1].set_ylim([
            self.plot.offsets[1][1],
            self.plot.offsets[1][3],
        ])
        coords = np.load(outline.coords_path)
        outline_poly = matplotlib.patches.Polygon(
            np.array([coords[:, 1], coords[:, 0]]).T,
            edgecolor="r",
            fill=False,
            lw=1
        )
        self.plot.axes[1].add_patch(outline_poly)

    def display_frame_right(self, current_outline):
        if current_outline.child_id1 is not None:
            selected_outline = self.outlines[
                self.outlines.outline_id == current_outline.child_id1
            ].iloc[0]
            self.selected_outlines.append(selected_outline)
            fidx = selected_outline.frame_idx
            offl, offt = self.get_offsets((
                selected_outline.centre_y,
                selected_outline.centre_x,
            ))

            if (current_outline.child_id2 is not None and
                    current_outline.child_id1 != current_outline.child_id2):
                selected_outline = self.outlines[
                    self.outlines.outline_id == current_outline.child_id2
                ].iloc[0]
                self.selected_outlines.append(selected_outline)
                offl2, offt2 = self.get_offsets((
                    selected_outline.centre_y,
                    selected_outline.centre_x,
                ))
                offt = (offt + offt2) / 2
                offl = (offl + offl2) / 2

            img = self.image_loader.load_frame(fidx, 0)

        elif current_outline.frame_idx + 1 < self.image_loader.num_frames:
            fidx = current_outline.frame_idx + 1
            offl, offt = self.plot.offsets[1][3], self.plot.offsets[1][0]
            img = self.image_loader.load_frame(fidx, 0)
        else:
            fidx = current_outline.frame_idx + 1
            img = np.zeros((self.region_halfwidth * 2, self.region_halfheight * 2))
            offt = 0
            offl = 0

        self.plot.axes[2].imshow(img, cmap="gray")
        self.plot.offsets[2] = (
            offt,
            offl + (self.region_halfwidth * 2),
            offt + (self.region_halfheight * 2),
            offl,
        )
        self.plot.axes[2].set_xlim([
            offt, offt + (self.region_halfheight * 2)
        ])
        self.plot.axes[2].set_ylim([
            offl + (self.region_halfwidth * 2), offl
        ])

        for outline in database.getOutlinesByFrameIdx(fidx, self.experiment_data.experiment_id):
            coords = np.load(outline.coords_path)
            kws = dict(
                edgecolor="k",
                lw=1,
                alpha=0.9,
                picker=True,
            )
            if (self.selected_outlines and
                    outline.outline_id in [
                        x.outline_id
                        for x in self.selected_outlines
                    ]):
                outline_poly = matplotlib.patches.Polygon(
                    np.array([coords[:, 1], coords[:, 0]]).T,
                    facecolor="g",
                    **kws
                )
                outline_poly.selected = True
            else:
                outline_poly = matplotlib.patches.Polygon(
                    np.array([coords[:, 1], coords[:, 0]]).T,
                    facecolor="y",
                    **kws
                )
                outline_poly.selected = False

            outline_poly.outline_id = outline.outline_id
            self.plot.axes[2].add_patch(outline_poly)

    def display_frame_left(self, current_outline):
        if current_outline.parent_id:
            outline = database.getOutlineById(current_outline.parent_id)
            offset_left, offset_top = self.get_offsets((
                outline.centre_y,
                outline.centre_x,
            ))
            img = self.image_loader.load_frame(outline.frame_idx, 0)
            self.plot.axes[0].imshow(img, cmap="gray")
            self.plot.offsets[0] = (
                offset_top,
                offset_left + self.region_halfwidth * 2,
                offset_top + self.region_halfheight * 2,
                offset_left,
            )
            coords = np.load(outline.coords_path)
            outline_poly = matplotlib.patches.Polygon(
                np.array([coords[:, 1], coords[:, 0]]).T,
                edgecolor="y",
                linestyle="--",
                fill=False,
                lw=1
            )
            self.plot.axes[0].add_patch(outline_poly)

        elif current_outline.frame_idx - 1 >= 0:
            img = self.image_loader.load_frame(current_outline.frame_idx - 1, 0)
            self.plot.axes[0].imshow(img, cmap="gray")
            self.plot.offsets[0] = self.plot.offsets[1]
        else:
            outline = None
            img = np.zeros((
                self.region_halfwidth * 2 + 2,
                self.region_halfheight * 2 + 2
            ))
            self.plot.axes[0].imshow(img, cmap="gray")
            self.plot.offsets[0] = (
                1, self.region_halfwidth * 2 + 1,
                self.region_halfwidth * 2 + 1,
                1,
            )

        self.plot.axes[0].set_xlim([
            self.plot.offsets[0][0],
            self.plot.offsets[0][2]
        ])
        self.plot.axes[0].set_ylim([
            self.plot.offsets[0][1],
            self.plot.offsets[0][3]
        ])

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

    @staticmethod
    def get_cell_area(coords):
        area = np.dot(
            coords[:, 0],
            np.roll(coords[:, 1], 1)
        ) - np.dot(
            coords[:, 1],
            np.roll(coords[:, 0], 1)
        )
        return area

    def previous_frame(self):
        if not self.lineage:
            return

        self.lineage.pop()
        self.display_frame()
        self.plot.setFocus()

    def next_frame(self):
        if not self.lineage or not self.selected_outlines:
            return

        if len(self.selected_outlines) == 1:
            # next frame
            outline = self.selected_outlines[0]
            coords_current = np.load(outline.coords_path)
            area_current = self.get_cell_area(coords_current)
            coords_prev = np.load(self.lineage[-1].coords_path)
            area_prev = self.get_cell_area(coords_prev)
            if (area_current / area_prev < 0.7 or
                    area_prev / area_current < 0.7):
                confirm = QtWidgets.QMessageBox().question(
                    self.window,
                    "Confirm assignment",
                    ("The outline you have selected causes a greater than 70% "
                     "change in cell area.\n"
                     "Are you sure you wish to assign it to the same cell?"),
                    QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes,
                    QtWidgets.QMessageBox.No,
                )
                if confirm == QtWidgets.QMessageBox.No:
                    return

            self.lineage.append(outline)
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
            for i, outline in enumerate(self.selected_outlines):
                if outline.outline_id == evt.artist.outline_id:
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

    def determine_state(self, selected_outlines):
        if selected_outlines is not None:
            self.selected_outlines = selected_outlines
            return True

        if not self.selected_outlines:
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
            self.status_bar.showMessage(
                "Something went wrong, more than  2 outlines were selected"
            )
            return False

        return True

    def write_lineage(self, selected_outlines=None):
        self.status_bar.showMessage("Writing lineage, please wait...")
        stop_condition = self.determine_state(selected_outlines)
        if not stop_condition:
            return False

        self.temp_window = QtWidgets.QDialog(self.window)
        self.temp_window.setGeometry(0, 0, self.max_width_px * 0.5, self.max_height_px * 0.9)
        self.temp_window.setWindowTitle("Please wait patiently")
        temp_layout = QtWidgets.QVBoxLayout()
        temp_layout.addWidget(QtWidgets.QLabel("Writing lineage, please wait..."))
        self.temp_window.setLayout(temp_layout)
        self.temp_window.show()

        cell_id = self.lineage[0].cell_id
        self.clear_extra_outlines(cell_id)
        rel_params = self.assign_relationships(cell_id)

        kwargs = {
            "cell_id": cell_id,
            "experiment_id": self.lineage[0].experiment_id,
            "start_frame_idx": int(self.lineage[0].frame_idx),
            "end_frame_idx": int(self.lineage[-1].frame_idx),
            "birth_observed": self.lineage[0].parent_id is not None,
            "division_observed": len(self.selected_outlines) == 2,
            "is_wildtype": rel_params["wildtype"],
            "first_outline_id": self.lineage[0].outline_id,
            "last_outline_id": self.lineage[-1].outline_id,
            "parent_cell_id": rel_params["parent_cell_id"],
            "child_cell_id1": rel_params["child_cell_id1"],
            "child_cell_id2": rel_params["child_cell_id2"],
        }
        database.insertCell(**kwargs)
        self.status_bar.showMessage("Finished writing lineage {0}".format(cell_id))
        self.temp_window.close()
        self.temp_window.deleteLater()
        self.temp_window = None
        self.plot.setFocus()
        return True

    def clear_extra_outlines(self, cell_id):
        self.status_bar.showMessage("Removing cell_id from extra outlines")
        for outline in self.lineage:
            database.updateOutlineById(
                outline.outline_id,
                cell_id=None,
            )

        existing_outlines = database.getOutlinesByCellId(cell_id)
        if existing_outlines:
            random_cell_id = str(uuid.uuid4())
            for outline in existing_outlines:
                database.updateOutlineById(
                    outline.outline_id,
                    cell_id=random_cell_id,
                )

    def assign_relationships(self, cell_id):
        self.status_bar.showMessage("Assigning parents and children")
        for outline_idx, outline in enumerate(self.lineage):
            if outline_idx == 0:
                parent_id = self.lineage[0].parent_id
            else:
                parent_id = self.lineage[outline_idx - 1].outline_id

            child_id1, child_id2 = self.get_children(outline_idx, outline)
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
            wildtype = parent_cell.is_wildtype if parent_cell else False

        if len(self.selected_outlines) == 2:
            child_cell_id1 = database.getOutlineById(
                self.selected_outlines[0].outline_id
            ).cell_id
            child_cell_id2 = database.getOutlineById(
                self.selected_outlines[1].outline_id
            ).cell_id

        existing = database.getCellById(cell_id)
        if existing:
            database.deleteCellById(cell_id)

        return {
            "wildtype": wildtype,
            "parent_cell_id": parent_cell_id,
            "child_cell_id1": child_cell_id1,
            "child_cell_id2": child_cell_id2,
        }

    def get_children(self, outline_idx, outline):
        if outline_idx < len(self.lineage) - 1:
            child_id1 = self.lineage[outline_idx + 1].outline_id
            child_id2 = ""
            database.updateOutlineById(
                child_id1,
                parent_id=outline.outline_id,
            )
        elif not self.selected_outlines:
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

        return child_id1, child_id2
