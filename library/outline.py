#!/usr/bin/env python3

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.widgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import PyQt5.QtWidgets as QtWidgets
import PyQt5.QtGui as QtGui
import PyQt5.QtCore as QtCore
import PyQt5.Qt as Qt
import seaborn as sns
import os
import tifffile
import uuid

from . import balloon
from . import database

sns.set_context("talk")
sns.set_style("white")

class Plotter(FigureCanvas):
    def __init__(self, parent_window, width, height, dpi, experiment_data, image_loader):
        fig = matplotlib.figure.Figure(figsize=(width, height), dpi=dpi)
        self.main_ax = fig.add_subplot(121)
        self.sub_ax = fig.add_subplot(122)

        self.decorate_axis(self.main_ax)
        self.decorate_axis(self.sub_ax)

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
        self.image_percentile = 1.0
        self.outline_store = os.path.join(
            "data", "outlines", self._data.experiment_hash
        )
        if not os.path.exists(self.outline_store):
            os.makedirs(self.outline_store)

        if not database.checkTable("outlines"):
            database.createOutlinesTable()

        self.image_loader = image_loader
        self.load_metadata()
        self.region_width, self.region_height = 75, 75
        # self.region_width, self.region_height = 100, 100

        self.cell_outlines = []
        self.cell_outline_text = []
        self.subfigure_patches = []
        self.dragging = False
        self.previous_id = None
        self.current_frame_idx = 0
        self.current_channel = 0

        self.mpl_connect("key_press_event", self._key_press_event)
        self.mpl_connect("button_press_event", self._button_press_event)
        self.mpl_connect("button_release_event", self._button_release_event)
        self.mpl_connect("motion_notify_event", self._motion_notify_event)

        self.main_frame = self.main_ax.imshow(self.load_frame(), cmap="gray")
        self.plot_existing_outlines()

    def load_metadata(self):
        self.num_frames = self.image_loader.num_frames
        self.num_channels = self.image_loader.num_channels

    def load_frame(self, frame_idx=None):
        if not frame_idx:
            frame_idx = self.current_frame_idx

        if frame_idx < 0 or frame_idx > (self.num_frames - 1):
            return

        return self.image_loader.load_frame(frame_idx, self.current_channel)

    def refresh(self):
        self.draw()

    def decorate_axis(self, ax):
        ax.axis("off")
        ax.set_aspect("equal")
        ax.autoscale("off")

    def plot_existing_outlines(self):
        self.main_ax.set_title("Frame = {0}".format(self.current_frame_idx + 1))
        while True:
            try:
                self.main_ax.lines.pop()
            except IndexError:
                break

        while True:
            try:
                t = self.cell_outlines.pop()
                t.remove()
            except IndexError:
                break

        while True:
            try:
                t = self.cell_outline_text.pop()
                t.remove()
            except IndexError:
                break

        for t in self.cell_outline_text:
            t.remove()

        outline_data = database.getOutlinesByFrameIdx(self.current_frame_idx, self._data.experiment_hash)
        for i, outline in enumerate(outline_data):
            if not os.path.exists(outline.coords_path):
                database.deleteOutlineById(outline.outline_id)
                continue

            c = np.load(outline.coords_path) + np.array([outline.offset_left, outline.offset_top])
            p = matplotlib.patches.Polygon(np.array([c[:, 1], c[:, 0]]).T, edgecolor="r", fill=False, lw=1)
            p._outline_id = outline.outline_id
            self.main_ax.add_patch(p)
            self.cell_outlines.append(p)
            centre = c.mean(axis=0)
            t = self.main_ax.text(
                centre[1], centre[0],
                "{0}".format(i + 1),
                verticalalignment="center",
                horizontalalignment="center",
                color="w",
            )
            self.cell_outline_text.append(t)

    def save_outline(self):
        coords_path = os.path.join(
            self.outline_store,
            "{0}.npy".format(self.outline_id)
        )

        data = {
            "outline_id": self.outline_id,
            "cell_id": self.cell_id,
            "experiment_id": self._data.experiment_id,
            "experiment_hash": self._data.experiment_hash,
            "image_path": self._data.image_path,
            "frame_idx": self.current_frame_idx,
            "coords_path": coords_path,
            "offset_left": self.offset_left,
            "offset_top": self.offset_top,
            "parent_id": self.previous_id or "",
        }
        if os.path.exists(coords_path):
            os.remove(coords_path)
        else:
            database.insertOutline(**data)

        coords = np.array([(n.x, n.y) for n in self.balloon_obj.nodes])
        np.save(coords_path, coords)

        if self.previous_id:
            database.addOutlineChild(self.previous_id, child1=self.outline_id)

        self.previous_id = str(self.outline_id)

    def fit_outline(self, roi, init_nodes=None):
        centre = [self.region_height, self.region_width]
        if init_nodes is None:
            self.outline_id = str(uuid.uuid4())
            radius = 5
            num_nodes = 10
            init_nodes = balloon.initial_nodes(centre, radius, num_nodes)

        self.balloon_obj = balloon.Balloon(init_nodes, roi)
        self.sub_ax.imshow(roi, cmap="gray")
        self.sub_ax.set_xlim([0, self.region_height * 2])
        self.sub_ax.set_ylim([self.region_width * 2, 0])
        self._plot_nodes()
        self.draw()

    def _plot_nodes(self, line=True, patches=True):
        nodes = self.balloon_obj.nodes
        if line:
            try:
                self.sub_ax.lines.pop(0)
            except IndexError:
                pass
            coords = np.array([(n.x, n.y) for n in nodes])
            self.sub_ax.plot(
                coords[:, 1],
                coords[:, 0],
                color="y",
                lw=1,
            )

        if patches:
            for i in range(len(self.subfigure_patches)):
                p = self.subfigure_patches.pop()
                p.remove()

            for n in nodes:
                patch = matplotlib.patches.Circle(
                    (n.y, n.x),
                    1,
                    fc="y",
                )
                patch.this_node = n
                self.subfigure_patches.append(patch)
                self.sub_ax.add_artist(patch)

        self.draw()

    def _key_press_event(self, evt):
        if evt.key == "left":
            if self.current_channel <= 0:
                return
            self.current_channel -= 1
            new_im = self.load_frame()
            self.main_frame.set_data(new_im)
            self.plot_existing_outlines()
            self.main_frame.set_clim([new_im.min(), new_im.max()])
            self.draw()

        elif evt.key == "right":
            if self.current_channel >= self.num_channels - 1:
                return
            self.current_channel += 1
            new_im = self.load_frame()
            self.main_frame.set_data(new_im)
            self.plot_existing_outlines()
            self.main_frame.set_clim([new_im.min(), new_im.max()])
            self.draw()

        elif evt.key == "up":
            if self.current_frame_idx >= self.num_frames - 1:
                return
            self.current_frame_idx += 1
            new_im = self.load_frame()
            self.main_frame.set_data(new_im)
            self.plot_existing_outlines()
            self.draw()

        elif evt.key == "down":
            if self.current_frame_idx <= 0:
                return
            self.current_frame_idx -= 1
            new_im = self.load_frame()
            self.main_frame.set_data(new_im)
            self.plot_existing_outlines()
            self.draw()

        elif evt.key == "r" and self.subfigure_patches:
            self._refine_event()

        elif evt.key == "d" and self.subfigure_patches:
            self._delete_event()

        elif evt.key == "." and self.subfigure_patches:
            self.balloon_obj.evolve(image_percentile=self.image_percentile)
            self._plot_nodes()

        elif evt.key == "enter" and self.subfigure_patches:
            # save outline
            self.save_outline()

            offset_centre = self.balloon_obj.get_centre()

            # clear plot
            self.balloon_obj = None
            self.dragging = False
            self.subfigure_patches = []
            self.sub_ax.clear()
            self.decorate_axis(self.sub_ax)
            self.draw()

            # fit next
            centre = [offset_centre[0] + self.offset_left,
                      offset_centre[1] + self.offset_top]
            self.offset_left = int(round(centre[0] - self.region_width))
            self.offset_top = int(round(centre[1] - self.region_height))
            self.current_frame_idx += 1
            bf_frame = self.load_frame()
            self.main_frame.set_data(bf_frame)
            self.plot_existing_outlines()
            roi = bf_frame[
                self.offset_left:self.offset_left + (self.region_width * 2),
                self.offset_top:self.offset_top + (self.region_height * 2)
            ]
            self.fit_outline(roi)

        else:
            print("Unknown key:", evt.key)

    def _button_press_event(self, evt):
        if self.parent().toolbar.mode:
            return

        if evt.inaxes == self.main_ax:
            # check is not in an existing outline
            hit = False
            for outline in self.cell_outlines:
                hit, _ = outline.contains(evt)
                if hit:
                    break

            if hit:
                outline_info = database.getOutlineById(outline._outline_id)
                self.previous_id = outline_info.parent_id
                self.cell_id = outline_info.cell_id
                self.offset_left = outline_info.offset_left
                self.offset_top = outline_info.offset_top
                roi = self.load_frame()[
                    self.offset_left:self.offset_left + (self.region_width * 2),
                    self.offset_top:self.offset_top + (self.region_height * 2)
                ]
                self.outline_id = outline_info.outline_id
                current_nodes = np.load(outline_info.coords_path)
                self.fit_outline(roi, init_nodes=current_nodes)

            else:
                self.previous_id = None
                self.cell_id = str(uuid.uuid4())
                centre = [evt.ydata, evt.xdata]
                self.offset_left = int(round(centre[0] - self.region_width))
                self.offset_top = int(round(centre[1] - self.region_height))
                roi = self.load_frame()[
                    self.offset_left:self.offset_left + (self.region_width * 2),
                    self.offset_top:self.offset_top + (self.region_height * 2)
                ]

                self.fit_outline(roi)

        elif evt.inaxes == self.sub_ax:
            if evt.button == 1:
                for p in self.subfigure_patches:
                    if p.contains_point((evt.x, evt.y)):
                        p.set_facecolor("r")
                        self.dragging = p
                        self.draw()
                        break
            elif evt.button == 3:
                for p in self.subfigure_patches:
                    if p.contains_point((evt.x, evt.y)):
                        self.balloon_obj.remove_node(p.this_node)
                        self._plot_nodes()
                        self.draw()

    def _button_release_event(self, evt):
        if not self.dragging:
            return

        self.dragging.set_facecolor("y")
        self.dragging.this_node.set_position((evt.ydata, evt.xdata))
        self.dragging.this_node.apply_changes()
        self._plot_nodes()
        self.dragging = False

    def _motion_notify_event(self, evt):
        if not self.dragging:
            return

        self.dragging.center = evt.xdata, evt.ydata
        self.draw()

    def _home_event(self):
        f = self.load_frame()
        self.main_ax.set_xlim([0, f.shape[0]])
        self.main_ax.set_ylim([f.shape[1], 0])
        self.sub_ax.set_xlim([0, self.region_width * 2])
        self.sub_ax.set_ylim([self.region_height * 2, 0])
        self.draw()


    def _delete_event(self):
        if not hasattr(self, "outline_id") or not self.subfigure_patches or self.outline_id is None:
            return

        alert = QtWidgets.QMessageBox()
        delete_confirm = alert.question(
            self.parent(),
            "Delete outline?",
            "Are you really sure you want to delete this outline permanently?"
        )
        if delete_confirm == QtWidgets.QMessageBox.Yes:
            database.deleteOutlineById(self.outline_id)
            self.sub_ax.clear()
            self.decorate_axis(self.sub_ax)
            self.outline_id = None
            self.balloon_obj = None
            self.dragging = False
            self.subfigure_patches = []
            self.plot_existing_outlines()
            self.draw()

    def _refine_event(self):
        if not hasattr(self, "outline_id") or not self.subfigure_patches or self.outline_id is None:
            return

        for i in range(10):
            self.balloon_obj.evolve(image_percentile=self.image_percentile)
        self._plot_nodes()


class Toolbar(NavigationToolbar):
    def __init__(self, figure_canvas, parent=None):
        self.toolitems = [
            ("Home", "Home", "home_large", "home_event"),
            (None, None, None, None),
            ("Pan", "Pan", "move_large", "pan"),
            ("Zoom", "Zoom", "zoom_to_rect_large", "zoom"),
            (None, None, None, None),
            ("Save", "Save view", "filesave_large", "save_figure"),
            (None, None, None, None),
            ("Delete", "Delete outline", "delete", "delete"),
            ("Refine", "Refine outline", "recycle", "refine"),
        ]
        NavigationToolbar.__init__(self, figure_canvas, parent=None)

    def _icon(self, name):
        path = os.path.join("resources", name)
        if not os.path.exists(path):
            path = os.path.join(self.basedir, name)

        pm = QtGui.QPixmap(path)
        if hasattr(pm, "setDevicePixelRatio"):
            pm.setDevicePixelRatio(self.canvas._dpi_ratio)

        return QtGui.QIcon(pm)

    def home_event(self, *args, **kwargs):
        self.canvas._home_event()

    def delete(self):
        self.canvas._delete_event()

    def refine(self):
        self.canvas._refine_event()

class Outliner:
    def __init__(self, experiment_data, image_loader):
        self.experiment_data = experiment_data
        self.image_loader = image_loader

    def set_screen_res(self, max_width_px, max_height_px, screen_dpi):
        self.max_width_px = max_width_px
        self.max_height_px = max_height_px
        self.screen_dpi = screen_dpi

    def _px_to_in(self, px):
        return px / self.screen_dpi

    def start_outlining(self, parent_window):
        if not hasattr(self, "max_width_px"):
            raise ValueError("Screen resolution has not been set")

        dim = (self._px_to_in(self.max_width_px * 0.5),
               self._px_to_in(self.max_height_px * 0.5))

        self.parent_window = parent_window
        self.window = QtWidgets.QDialog(self.parent_window)
        self.window.setModal(True)
        self.window.setGeometry(self.max_width_px // 10, self.max_height_px // 10, self.max_width_px * 0.5, self.max_height_px * 0.5)
        self.window.setWindowTitle("Outline cells")

        main_layout = QtWidgets.QVBoxLayout()

        menubar = QtWidgets.QMenuBar(self.window)
        file_menu = menubar.addMenu("&File")
        quit_action = QtWidgets.QAction("&Close", menubar)
        quit_action.triggered[bool].connect(lambda: self.window.close())
        file_menu.addAction(quit_action)
        main_layout.setMenuBar(menubar)

        self.plot = Plotter(
            self.window,
            dim[0],
            dim[1],
            dpi=self.screen_dpi,
            experiment_data=self.experiment_data,
            image_loader=self.image_loader,
        )
        self.plot.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.plot.setFocus()

        self.window.toolbar = Toolbar(self.plot, self.window)

        tool_layout = QtWidgets.QVBoxLayout()
        tool_layout.addWidget(self.window.toolbar)

        # add tolerance box thing
        label = QtWidgets.QLabel("Tolerance:")
        self.tolerance_widget = QtWidgets.QLineEdit()
        self.tolerance_widget.setText(str(self.plot.image_percentile))
        self.tolerance_widget.textChanged[str].connect(lambda text: self._submit_tolerance(text))
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(label)
        layout.addWidget(self.tolerance_widget)
        tool_layout.addLayout(layout)

        main_layout.addLayout(tool_layout)
        main_layout.addWidget(self.plot)

        self.window.setLayout(main_layout)
        self.window.show()

    def _submit_tolerance(self, text):
        try:
            if text:
                self.plot.image_percentile = float(text)
            else:
                self.plot.image_percentile = 0.0
                self.tolerance_widget.setText(str(self.plot.image_percentile))
        except ValueError:
            self.tolerance_widget.setText(str(self.plot.image_percentile))
