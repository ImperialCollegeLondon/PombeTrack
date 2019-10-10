#!/usr/bin/env python3

import os
import time
import uuid

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

import matplotlib
import matplotlib.figure
import matplotlib.path
import matplotlib.widgets
import numpy as np
import PyQt5.QtCore as QtCore
import PyQt5.QtGui as QtGui
import PyQt5.QtWidgets as QtWidgets
import seaborn as sns

from . import balloon
from . import database
from . import segmentation

matplotlib.use('Qt5Agg')
sns.set_context("talk")
sns.set_style("white")

class Plotter(FigureCanvas):
    def __init__(self, parent_window, width, height, dpi, experiment_data, image_loader, status_bar):
        self.fig = matplotlib.figure.Figure(figsize=(width, height), dpi=dpi)
        self.main_ax = self.fig.add_subplot(121)
        self.sub_ax = self.fig.add_subplot(122)

        self.decorate_axis(self.main_ax)
        self.decorate_axis(self.sub_ax)

        FigureCanvas.__init__(self, self.fig)
        self.parent_window = parent_window
        self.setParent(self.parent_window)
        self.status_bar = status_bar

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
            "data", "outlines", self._data.experiment_id
        )
        if not os.path.exists(self.outline_store):
            os.makedirs(self.outline_store)

        if not database.checkTable("outlines"):
            database.createOutlinesTable()

        if not database.checkTable("cells"):
            database.createCellsTable()

        self.image_loader = image_loader
        self.load_metadata()
        self.region_halfwidth, self.region_halfheight = 75, 75
        # self.region_halfwidth, self.region_halfheight = 100, 100

        self.cell_outlines = []
        self.sub_outlines = []
        self.subfigure_patches = []
        self.main_dragging = False
        self.main_dragging_rect = None
        self.main_dragging_background = None
        self.selected_outlines = []
        self.dragging = False
        self.sub_dragging_background = None

        self.previous_id = None
        self.current_frame_idx = 0
        self.current_channel = 0
        self.current_slice = 0
        self.cell_id = None
        self.outline_id = None
        self.full_coords = None
        self.sub_outline = None
        self.balloon_obj = None
        self.centre_x, self.centre_y = None, None
        self.offset_left, self.offset_top = None, None

        self.mpl_connect("key_press_event", self._key_press_event)
        self.mpl_connect("button_press_event", self._button_press_event)
        self.mpl_connect("button_release_event", self._button_release_event)
        self.mpl_connect("motion_notify_event", self._motion_notify_event)

        self.current_status = None
        self.main_frame = self.main_ax.imshow(self.load_frame(), cmap="gray")
        self.plot_existing_outlines()

    def set_status(self, text=None, clear=False, status=None):
        if not clear and text is not None:
            self.status_bar.showMessage(text)
        else:
            self.status_bar.clearMessage()

        if not status and self.current_status:
            QtGui.QGuiApplication.restoreOverrideCursor()
            # QtGui.QGuiApplication.setOverrideCursor(QtGui.QCursor(
            #     QtCore.Qt.ArrowCursor
            # ))
        elif status == "working" and self.current_status != status:
            QtGui.QGuiApplication.setOverrideCursor(QtGui.QCursor(
                QtCore.Qt.WaitCursor
            ))

        self.current_status = status
        self.status_bar.repaint()

    def automatic_segmentation(self, display=True):
        self.set_status(
            "Preprocessing",
            status="working",
        )
        # load_frame: frame, z-slice, channel
        img_mid = self.load_frame(
            self.current_frame_idx,
            int(np.floor(self.num_slices / 2)),
            0
        )
        img_up = self.load_frame(
            self.current_frame_idx,
            int(np.floor(self.num_slices / 2) - 1),
            0
        )
        img = np.maximum(img_mid, img_up)


        img_pp = segmentation.preprocessing(img)
        img_i = segmentation.find_cellinterior(img_pp)
        img_wat = segmentation.find_watershed(img_i)
        #  bd = segmentation.find_bd(im_wat)

        if display:
            background = self.fig.canvas.copy_from_bbox(self.main_ax.bbox)

        cell_count = 0

        for index in range(1, img_wat.max()+1):
            self.set_status(
                "Processing cell {0} of {1}".format(index, img_wat.max()),
                status="working",
            )
            img_ii = img_wat == index
            if (np.any(np.asarray(img_ii.nonzero()) == 0) or
                    np.any(np.asarray(img_ii.nonzero()) == 2047)):
                continue

            img_ii_bd = segmentation.find_boundaries(img_ii, mode="inner")
            # Sort in radial
            bd_ii_sorted = segmentation.sort_in_order(img_ii_bd)

            # Define the balloon object
            balloon_obj, origin_y, origin_x, halfwidth = segmentation.find_balloon_obj(
                bd_ii_sorted.astype(int)[::5],
                img,
            )
            #  Test if the cell exists
            overlap = False
            outline_data = database.getOutlinesByFrameIdx(
                self.current_frame_idx,
                self._data.experiment_id
            )
            for outline in outline_data:
                if not os.path.exists(outline.coords_path):
                    continue
                if (outline.centre_x not in range(origin_x, origin_x + 2 * halfwidth) or
                        outline.centre_y not in range(origin_y, origin_y + 2 * halfwidth)):
                    continue

                polygonpath = matplotlib.path.Path(
                    np.append(balloon_obj.get_coordinates(),
                              balloon_obj.get_coordinates()[1, :].reshape(1, 2),
                              axis=0),
                    closed=True
                )
                if polygonpath.contains_point([outline.centre_y-origin_y,
                                               outline.centre_x - origin_x]):
                    overlap = True

            if overlap:
                continue

            # Evolve the contour
            try:
                sensitivity = self.image_percentile
                area_init = balloon_obj.get_area()
                for _ in range(20):
                    balloon_obj.evolve(image_percentile=sensitivity)
                    if (balloon_obj.get_area() > 1.5 * area_init or
                            balloon_obj.get_area() < 0.5 * area_init):
                        raise ValueError

                self.full_coords = balloon_obj.get_coordinates() + [origin_y, origin_x]
                self.outline_id = str(uuid.uuid4())
                self.cell_id = str(uuid.uuid4())
                self.centre_y, self.centre_x = self.full_coords.mean(axis=0).astype(int)
                self.save_outline(auto=True)

                # Draw the cell
                if display:
                    self.fig.canvas.restore_region(background)
                    outline_poly = matplotlib.patches.Polygon(
                        np.array([self.full_coords[:, 1],
                                  self.full_coords[:, 0]]).T,
                        edgecolor="r",
                        fill=False,
                        lw=1
                    )
                    outline_poly.outline_id = self.outline_id
                    self.main_ax.add_patch(outline_poly)
                    self.cell_outlines.append(outline_poly)
                    self.main_ax.draw_artist(outline_poly)
                    self.fig.canvas.blit(self.main_ax.bbox)
                    background = self.fig.canvas.copy_from_bbox(self.main_ax.bbox)
            except ValueError:
                continue

            cell_count += 1
            self.draw()
        self.set_status(text="Auto segmentation finished with {0} cells".format(cell_count))



    def load_metadata(self):
        self.num_frames = self.image_loader.num_frames
        self.num_channels = self.image_loader.num_channels
        self.num_slices = self.image_loader.num_slices

    def load_frame(self, frame_idx=None, slice_idx=None, channel_idx=None):
        if frame_idx is None:
            frame_idx = self.current_frame_idx

        if channel_idx is None:
            channel_idx = self.current_channel

        if slice_idx is None:
            slice_idx = self.current_slice

        if frame_idx < 0 or frame_idx > (self.num_frames - 1):
            return np.zeros((100, 100))

        return self.image_loader.load_frame(frame_idx, slice_idx, channel_idx)

    def refresh(self):
        self.draw()

    @staticmethod
    def decorate_axis(axis):
        axis.axis("off")
        axis.set_aspect("equal")
        axis.autoscale("off")

    def clear_sub_ax(self):
        self.sub_ax.clear()
        self.decorate_axis(self.sub_ax)
        self.outline_id = None
        self.balloon_obj = None
        self.main_dragging = False
        self.dragging = False
        self.subfigure_patches = []
        self.sub_outline = []

    def clear_sub_outlines(self):
        while True:
            try:
                outline_poly = self.sub_outlines.pop()
                outline_poly.remove()
            except ValueError:
                pass
            except IndexError:
                break

    def plot_existing_outlines(self):
        if self._data.image_mode == "static":
            title = "Image #{0}".format(self.current_frame_idx + 1)
        else:
            title = "F={0}".format(self.current_frame_idx + 1)

        if self._data.num_slices > 1:
            title += " Z={0}".format(self.current_slice + 1)

        if self._data.num_channels > 1:
            title += " C={0}".format(self.current_channel + 1)

        self.main_ax.set_title(title)
        while True:
            try:
                self.main_ax.lines.pop()
            except IndexError:
                break

        while True:
            try:
                outline_poly = self.cell_outlines.pop()
                outline_poly.remove()
            except IndexError:
                break

        selected_ids = [x.outline_id for x in self.selected_outlines]
        fresh_selections = []
        outline_data = database.getOutlinesByFrameIdx(
            self.current_frame_idx,
            self._data.experiment_id,
        )
        for outline in outline_data:
            if not os.path.exists(outline.coords_path):
                # database.deleteOutlineById(outline.outline_id)
                continue

            coords = np.load(outline.coords_path)
            outline_poly = matplotlib.patches.Polygon(
                np.array([coords[:, 1], coords[:, 0]]).T,
                edgecolor="r",
                fill=False,
                lw=1
            )
            outline_poly.outline_id = outline.outline_id
            outline_poly.cell_id = outline.cell_id

            if outline.outline_id in selected_ids:
                prev_outline = self.selected_outlines[
                    selected_ids.index(outline.outline_id)
                ]
                if hasattr(prev_outline, "is_modified") and prev_outline.is_modified:
                    outline_poly.set_xy(prev_outline.get_xy())
                    outline_poly.is_modified = True
                outline_poly.set_edgecolor("yellow")
                outline_poly.is_selected = True
                fresh_selections.append(outline_poly)

            self.main_ax.add_patch(outline_poly)
            self.cell_outlines.append(outline_poly)

        self.selected_outlines = fresh_selections

    def save_outline(self, auto=False, explicit=None):
        if auto:
            coords = self.full_coords
        elif explicit:
            self.outline_id = explicit.outline_id
            outline_info = database.getOutlineById(self.outline_id)
            xy_coords = explicit.get_xy()
            coords = np.array([xy_coords[:, 1], xy_coords[:, 0]]).T
            self.cell_id = explicit.cell_id
            self.previous_id = outline_info.parent_id  # pylint: disable=no-member
            self.centre_y, self.centre_x = coords.mean(axis=0)
        else:
            coords_offset = np.array([(n.x, n.y) for n in self.balloon_obj.nodes])
            # Determine cell centre
            coords = coords_offset + np.array([self.offset_left, self.offset_top])
            self.centre_y, self.centre_x = coords.mean(axis=0)

        coords_path = os.path.join(
            self.outline_store,
            "{0}.npy".format(self.outline_id)
        )
        data = {
            "outline_id": self.outline_id,
            "cell_id": self.cell_id,
            "experiment_num": self._data.experiment_num,
            "experiment_id": self._data.experiment_id,
            "image_path": self._data.image_path,
            "frame_idx": self.current_frame_idx,
            "coords_path": coords_path,
            "parent_id": self.previous_id or "",
            "centre_x":int(self.centre_x),
            "centre_y":int(self.centre_y),
        }
        if self._data.image_mode == "static":
            data["parent_id"] = ""

        if os.path.exists(coords_path):
            os.remove(coords_path)
            database.updateOutlineById(
                self.outline_id,
                centre_x=int(self.centre_x),
                centre_y=int(self.centre_y),
            )
        else:
            database.insertOutline(**data)

        np.save(coords_path, coords)

        if self.previous_id:
            database.addOutlineChild(self.previous_id, child1=self.outline_id)

        database.updateExperimentById(
            self._data.experiment_id,
            verified=False,
        )
        database.deleteCellById(self.cell_id)

    def fit_outline(self, roi, init_nodes=None, centre_offset_left=0, centre_offset_top=0):
        if init_nodes is None:
            centre = [
                self.region_halfwidth - centre_offset_left,
                self.region_halfheight - centre_offset_top,
            ]
            self.outline_id = str(uuid.uuid4())
            radius = 5
            num_nodes = 10
            init_nodes = balloon.initial_nodes(centre, radius, num_nodes)

        self.balloon_obj = balloon.Balloon(init_nodes, roi)
        self.sub_ax.imshow(self.load_frame(), cmap="gray")
        self.sub_ax.set_ylim([self.offset_left + self.region_halfwidth * 2, self.offset_left])
        self.sub_ax.set_xlim([self.offset_top, self.region_halfheight * 2 + self.offset_top])
        self.sub_ax.set_title("Frame = {0}".format(self.current_frame_idx + 1))
        self._plot_nodes()

        self.clear_sub_outlines()
        outline_data = database.getOutlinesByFrameIdx(
            self.current_frame_idx,
            self._data.experiment_id,
        )
        for outline in outline_data:
            if not os.path.exists(outline.coords_path) or outline.outline_id == self.outline_id:
                continue

            coords = np.load(outline.coords_path)
            poly = matplotlib.patches.Polygon(
                np.array([coords[:, 1], coords[:, 0]]).T,
                edgecolor="r",
                fill=False,
                lw=1,
            )
            poly.outline_id = outline.outline_id
            self.sub_ax.add_patch(poly)
            self.sub_outlines.append(poly)

        self.draw()

    def _plot_nodes(self, line=True, patches=True):
        nodes = self.balloon_obj.nodes
        if line:
            try:
                self.sub_ax.lines.pop(0)
            except IndexError:
                pass
            coords = np.array([
                (n.x + self.offset_left,
                 n.y + self.offset_top)
                for n in nodes
            ])
            self.sub_ax.plot(
                coords[:, 1],
                coords[:, 0],
                color="y",
                lw=1,
            )

        if patches:
            for _ in range(len(self.subfigure_patches)):
                patch = self.subfigure_patches.pop()
                patch.remove()

            for node_idx, node in enumerate(nodes):
                patch = matplotlib.patches.Circle(
                    (node.y + self.offset_top,
                     node.x + self.offset_left),
                    1.5,
                    fc="y",
                )
                patch.this_node = node
                patch.node_idx = node_idx
                self.subfigure_patches.append(patch)
                self.sub_ax.add_artist(patch)

        self.draw()

    def slice_change(self, delta):
        if delta < 0 and self.current_slice <= 0:
            return

        if delta > 0 and self.current_slice >= self.num_slices - 1:
            return

        self.current_slice += delta
        new_im = self.load_frame()
        self.main_frame.set_data(new_im)
        self.plot_existing_outlines()
        self.main_frame.set_clim([new_im.min(), new_im.max()])
        self.draw()

    def channel_change(self, delta):
        if delta < 0 and self.current_channel <= 0:
            return

        if delta > 0 and self.current_channel >= self.num_channels - 1:
            return

        self.current_channel += delta
        new_im = self.load_frame()
        self.main_frame.set_data(new_im)
        self.plot_existing_outlines()
        self.main_frame.set_clim([new_im.min(), new_im.max()])
        self.draw()

    def frame_change(self, delta):
        if delta < 0 and self.current_frame_idx <= 0:
            return

        if delta > 0 and self.current_frame_idx >= self.num_frames - 1:
            return

        if not self.check_selected_outlines():
            return

        self.current_frame_idx += delta
        new_im = self.load_frame()
        self.main_frame.set_data(new_im)
        self.deselect_outlines()
        self.plot_existing_outlines()
        self.main_frame.set_clim([new_im.min(), new_im.max()])
        self.draw()

    def _key_press_event(self, evt):
        if evt.key == "left" or evt.key == "a":
            self.channel_change(-1)

        elif evt.key == "right" or evt.key == "d":
            self.channel_change(1)

        elif evt.key == "up" or evt.key == "w":
            self.frame_change(1)

        elif evt.key == "down" or evt.key == "s":
            self.frame_change(-1)

        elif evt.key == "q":
            self.slice_change(-1)

        elif evt.key == "e":
            self.slice_change(1)

        elif evt.key == "r" and len(self.selected_outlines) > 1:
            self.refine_multi()

        elif evt.key == "r" and self.subfigure_patches:
            self.refine_event()

        elif evt.key == "delete" and len(self.selected_outlines) > 1:
            self.delete_multi()

        elif evt.key == "delete" and self.subfigure_patches:
            self.delete_event()

        elif evt.key == "." and len(self.selected_outlines) > 1:
            self.refine_multi(1)

        elif evt.key == "." and self.subfigure_patches:
            self.refine_event(1)

        elif evt.key == "enter" and len(self.selected_outlines) > 1:
            self.accept_multi()

        elif evt.key == "enter" and self.subfigure_patches:
            self.accept_event()

        elif evt.key == "shift" or evt.key == "control":
            pass

        else:
            print("Unknown key:", evt.key)

    def get_offsets(self, centre):
        offset_left = int(round(centre[0] - self.region_halfwidth))
        offset_top = int(round(centre[1] - self.region_halfheight))
        centre_offset_left, centre_offset_top = 0, 0
        img = self.load_frame()
        if offset_left < 0:
            centre_offset_left = -offset_left
            offset_left = 0
        elif offset_left >= img.shape[0] - (self.region_halfwidth * 2):
            centre_offset_left = img.shape[0] - (self.region_halfwidth * 2) - offset_left
            offset_left = img.shape[0] - (self.region_halfwidth * 2)

        if offset_top < 0:
            centre_offset_top = -offset_top
            offset_top = 0
        elif offset_top >= img.shape[1] - (self.region_halfheight * 2):
            centre_offset_top = img.shape[1] - (self.region_halfheight * 2) - offset_top
            offset_top = img.shape[1] - (self.region_halfheight * 2)
        del img

        return offset_left, offset_top, centre_offset_left, centre_offset_top

    def _button_press_event(self, evt):
        if self.parent().toolbar.mode:
            return

        if evt.inaxes == self.main_ax:
            self.main_dragging = [evt.xdata, evt.ydata]
            self.main_dragging_rect = matplotlib.patches.Rectangle(
                self.main_dragging,
                0, 0,
                edgecolor="yellow",
                facecolor=(0.9, 0.9, 0.7, 0.3),
            )
            self.main_ax.add_patch(self.main_dragging_rect)
            self.main_dragging_background = self.fig.canvas.copy_from_bbox(
                self.main_ax.bbox
            )
            self.draw()

        elif evt.inaxes == self.sub_ax:
            if evt.button == 1:
                for patch in self.subfigure_patches:
                    if patch.contains_point((evt.x, evt.y)):
                        patch.set_visible(False)
                        self.draw()
                        self.sub_dragging_background = self.fig.canvas.copy_from_bbox(
                            self.sub_ax.bbox
                        )
                        self.sub_ax.lines[0].set_color("red")
                        patch.set_visible(True)
                        self.dragging = patch
                        patch.set_facecolor("r")
                        # self.draw()
                        break

            elif evt.button == 3:
                for patch in self.subfigure_patches:
                    if patch.contains_point((evt.x, evt.y)):
                        self.balloon_obj.remove_node(patch.this_node)
                        self._plot_nodes()
                        self.draw()

    def check_selected_outlines(self):
        if sum([hasattr(x, "is_modified") and x.is_modified or False
                for x in self.selected_outlines]) > 0:
            alert = QtWidgets.QMessageBox()
            message = ("The selected outlines have been modified"
                       "\nWould you like to save them?")
            add_confirm = alert.question(
                self.parent(),
                "Save changes?",
                message,
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Cancel,
                QtWidgets.QMessageBox.Cancel,
            )
            if add_confirm == QtWidgets.QMessageBox.Yes:
                for outline_num, outline in enumerate(self.selected_outlines):
                    self.set_status(
                        "Saving {0} of {1} outlines".format(
                            outline_num + 1,
                            len(self.selected_outlines),
                        )
                    )
                    self.save_outline(auto=False, explicit=outline)
                    outline.is_modified = False
                self.set_status(
                    "Modified outlines saved ({0} outlines)".format(
                        len(self.selected_outlines),
                    )
                )
                self.plot_existing_outlines()
                self.draw()
                return True

            if add_confirm == QtWidgets.QMessageBox.No:
                self.set_status(
                    "Modified outlines discarded ({0} outlines)".format(
                        len(self.selected_outlines),
                    )
                )

                for outline in self.selected_outlines:
                    outline.is_modified = False

                self.plot_existing_outlines()
                self.draw()
                return True

            if add_confirm == QtWidgets.QMessageBox.Cancel:
                self.set_status(clear=True)
                return False

        return True

    def deselect_outlines(self):
        for outline in self.selected_outlines:
            outline.set_edgecolor("red")
            outline.is_selected = False
            outline.is_modified = False

        self.selected_outlines = []

    def _button_release_main(self, evt):
        if not evt.xdata or not evt.ydata or evt.inaxes != self.main_ax:
            rect_width = self.main_dragging_rect.get_width()
            rect_height = self.main_dragging_rect.get_height()
        else:
            rect_width = evt.xdata - self.main_dragging[0]
            rect_height = evt.ydata - self.main_dragging[1]

        keyboard_mods = QtGui.QGuiApplication.queryKeyboardModifiers()
        # if keyboard_mods & (QtCore.Qt.ShiftModifier | QtCore.Qt.ControlModifier):
        if keyboard_mods & QtCore.Qt.ShiftModifier:
            # additive mode
            additive_flag = True
            reductive_flag = False

        elif keyboard_mods & QtCore.Qt.ControlModifier:
            additive_flag = False
            reductive_flag = True
            if not self.check_selected_outlines():
                self.main_dragging = False
                self.main_dragging_rect.remove()
                del self.main_dragging_rect
                self.draw()
                return
        else:
            additive_flag = False
            reductive_flag = False
            if not self.check_selected_outlines():
                self.main_dragging = False
                self.main_dragging_rect.remove()
                del self.main_dragging_rect
                self.draw()
                return

            self.deselect_outlines()

        # determine whether any outline within the bounds
        new_outlines = self._get_outlines_in_rect(evt)

        if not reductive_flag:
            for outline in new_outlines:
                outline.set_edgecolor("yellow")
                outline.is_selected = True
                self.selected_outlines.append(outline)

        elif reductive_flag:
            for outline in new_outlines:
                if hasattr(outline, "is_selected") and outline.is_selected:
                    outline.set_edgecolor("red")
                    outline.is_selected = False
                    self.selected_outlines.pop(
                        self.selected_outlines.index(outline)
                    )

        self.main_dragging = False
        self.main_dragging_rect.remove()
        del self.main_dragging_rect

        if len(self.selected_outlines) > 1:
            self.clear_sub_ax()

        elif len(self.selected_outlines) == 1:
            self._select_hit(self.selected_outlines[0])

        elif not self.selected_outlines:
            self.clear_sub_ax()
            if (not additive_flag and not reductive_flag and
                    rect_width < 10 and rect_height < 10):
                self.previous_id = None
                self.cell_id = str(uuid.uuid4())
                centre = [evt.ydata, evt.xdata]
                (self.offset_left, self.offset_top,
                 centre_offset_left, centre_offset_top) = self.get_offsets(centre)

                roi = self.load_frame(channel_idx=0)[
                    self.offset_left:self.offset_left + (self.region_halfwidth * 2),
                    self.offset_top:self.offset_top + (self.region_halfheight * 2)
                ]

                self.fit_outline(
                    roi,
                    centre_offset_left=centre_offset_left,
                    centre_offset_top=centre_offset_top,
                )

        self.draw()

    def _get_outlines_in_rect(self, evt):
        rect_path = self.main_dragging_rect.get_path()
        rect_transform = self.main_dragging_rect.get_transform()
        true_rect = rect_transform.transform_path(rect_path)
        new_outlines = []
        for outline in self.cell_outlines:
            out_transform = outline.get_transform()
            out_path = outline.get_path()
            true_outline = out_transform.transform_path(out_path)
            contains = true_rect.contains_points(
                true_outline.vertices
            )
            if sum(contains) > 0:
                new_outlines.append(outline)

        if not new_outlines:
            hit = self._check_hit(evt)
            if hit:
                new_outlines.append(hit)

        return new_outlines

    def _button_release_event(self, evt):
        if self.parent().toolbar.mode:
            return

        if self.main_dragging:
            self._button_release_main(evt)
            return

        if evt.inaxes == self.sub_ax:
            if not self.dragging:
                return

            self.dragging.set_facecolor("y")
            self.dragging.this_node.set_position((
                evt.ydata - self.offset_left,
                evt.xdata - self.offset_top,
            ))
            self.dragging.this_node.apply_changes()
            self._plot_nodes()
            self.dragging = False
            return

    def _check_hit(self, evt):
        for outline in self.cell_outlines:
            hit, _ = outline.contains(evt)
            if hit:
                return outline
        return False

    def _select_hit(self, outline):
        outline.set_edgecolor("yellow")
        outline.is_selected = True
        outline_info = database.getOutlineById(outline.outline_id)
        self.previous_id = outline_info.parent_id  # pylint: disable=no-member
        self.cell_id = outline_info.cell_id  # pylint: disable=no-member
        centre = outline_info.centre_y, outline_info.centre_x  # pylint: disable=no-member
        (self.offset_left, self.offset_top,
         centre_offset_left, centre_offset_top) = self.get_offsets(centre)
        roi = self.load_frame(channel_idx=0)[
            self.offset_left:self.offset_left + (self.region_halfwidth * 2),
            self.offset_top:self.offset_top + (self.region_halfheight * 2),
        ]
        self.outline_id = outline_info.outline_id  # pylint: disable=no-member
        current_nodes = np.load(outline_info.coords_path)  # pylint: disable=no-member
        self.fit_outline(
            roi,
            init_nodes=current_nodes - np.array([self.offset_left, self.offset_top]),
            centre_offset_left=centre_offset_left,
            centre_offset_top=centre_offset_top,
        )
        self.balloon_obj.refining_cycles = 1

    def _motion_notify_event(self, evt):
        if self.parent().toolbar.mode:
            return

        if self.main_dragging and evt.inaxes == self.main_ax:
            self.fig.canvas.restore_region(self.main_dragging_background)
            new_width = evt.xdata - self.main_dragging[0]
            new_height = evt.ydata - self.main_dragging[1]
            self.main_dragging_rect.set_width(new_width)
            self.main_dragging_rect.set_height(new_height)
            self.main_ax.draw_artist(self.main_dragging_rect)
            self.fig.canvas.blit(self.main_ax.bbox)

        elif self.dragging and evt.inaxes == self.sub_ax:
            self.fig.canvas.restore_region(self.sub_dragging_background)
            self.dragging.center = evt.xdata, evt.ydata
            xdata, ydata = self.sub_ax.lines[0].get_data()
            xdata[self.dragging.node_idx] = evt.xdata
            ydata[self.dragging.node_idx] = evt.ydata
            self.sub_ax.lines[0].set_data(xdata, ydata)
            self.sub_ax.draw_artist(self.sub_ax.lines[0])
            self.sub_ax.draw_artist(self.dragging)
            self.fig.canvas.blit(self.sub_ax.bbox)

    def home_event(self):
        frame = self.load_frame()
        self.main_ax.set_xlim([0, frame.shape[0]])
        self.main_ax.set_ylim([frame.shape[1], 0])
        self.sub_ax.set_xlim([0, self.region_halfwidth * 2])
        self.sub_ax.set_ylim([self.region_halfheight * 2, 0])
        self.draw()

    @staticmethod
    def get_area(coords):
        area = np.dot(
            coords[:, 0],
            np.roll(coords[:, 1], 1)
        ) - np.dot(
            coords[:, 1],
            np.roll(coords[:, 0], 1)
        )
        return area * 0.5

    def accept_event(self):
        if not hasattr(self, "outline_id") or not self.subfigure_patches or self.outline_id is None:
            return

        self.set_status("Saving outline", status="working")
        # check balloon object has been refined at least once
        if self.balloon_obj.refining_cycles == 0:
            alert = QtWidgets.QMessageBox()
            message = ("The outline you are about to add has not been refined."
                       "\nAre you sure you want to add it?")
            add_confirm = alert.question(
                self.parent(),
                "Add outline?",
                message,
            )
            if add_confirm != QtWidgets.QMessageBox.Yes:
                self.set_status("Save aborted")
                return

        # check overlap with other outlines?
        # objects already exist (red lines) so can compare perhaps?
        x_points, y_points = self.sub_ax.lines[0].get_data()
        existing_points = np.row_stack([x_points, y_points]).T
        for outline in self.sub_outlines:
            contains_new = sum(outline.get_path().contains_points(
                existing_points,
            ))
            if contains_new > 0:
                alert = QtWidgets.QMessageBox()
                message = ("The outline you are about to add may overlap with "
                           "existing outlines.\n"
                           "Are you sure you would like to add it?")
                add_confirm = alert.question(
                    self.parent(),
                    "Add outline?",
                    message,
                )
                if add_confirm != QtWidgets.QMessageBox.Yes:
                    self.set_status("Save aborted")
                    return

        if self._data.image_mode == "movie":
            area_diff = self._area_diff()
            if area_diff > 0.3:
                alert = QtWidgets.QMessageBox()
                message = ("The outline you are about to add has a large "
                           "difference in cell area compared to its "
                           "predecessor in the previous frame ({0:.0f}%).\n"
                           "Are you sure you want to add it?")
                add_confirm = alert.question(
                    self.parent(),
                    "Add outline?",
                    message.format(area_diff * 100),
                )
                if add_confirm != QtWidgets.QMessageBox.Yes:
                    self.set_status("Save aborted")
                    return

            self.save_outline()
            self.deselect_outlines()

            offset_centre = self.balloon_obj.get_centre()

            # clear plot
            self.previous_id = str(self.outline_id)
            self.clear_sub_ax()

            # fit next
            centre = [offset_centre[0] + self.offset_left,
                      offset_centre[1] + self.offset_top]

            (self.offset_left, self.offset_top,
             centre_offset_left, centre_offset_top) = self.get_offsets(centre)

            if self.current_frame_idx == self.num_frames - 1:
                bf_frame = self.load_frame()
                self.main_frame.set_data(bf_frame)
                self.previous_id = None
                self.outline_id = None
                self.plot_existing_outlines()
                self.clear_sub_outlines()
                self.draw()
                return

            self.current_frame_idx += 1
            bf_frame = self.load_frame()
            self.main_frame.set_data(bf_frame)
            self.plot_existing_outlines()
            self.clear_sub_outlines()

            roi = self.load_frame(channel_idx=0)[
                self.offset_left:self.offset_left + (self.region_halfwidth * 2),
                self.offset_top:self.offset_top + (self.region_halfheight * 2)
            ]
            self.fit_outline(
                roi,
                centre_offset_left=centre_offset_left,
                centre_offset_top=centre_offset_top,
            )
            self.set_status(clear=True)

        elif self._data.image_mode == "static":
            self.save_outline()
            self.clear_sub_ax()

            # update existing outlines
            self.deselect_outlines()
            self.plot_existing_outlines()
            self.draw()
            self.set_status("Outline saved")

    def accept_multi(self):
        for outline_num, outline in enumerate(self.selected_outlines):
            self.set_status(
                "Saving {0: 2d} of {1} outlines".format(
                    outline_num + 1,
                    len(self.selected_outlines),
                ),
                status="working"
            )
            self.save_outline(auto=False, explicit=outline)

        self.set_status(
            "Outlines saved ({0} outlines)".format(
                len(self.selected_outlines)
            )
        )

        self.deselect_outlines()
        self.draw()

    def _area_diff(self):
        # check difference in total area is small
        if self.previous_id:
            previous_outline_entry = database.getOutlineById(self.previous_id)
            previous_outline = np.load(previous_outline_entry.coords_path)  # pylint: disable=no-member
            if previous_outline_entry.cell_id == self.cell_id:  # pylint: disable=no-member
                previous_area = self.get_area(previous_outline)
                current_area = self.get_area(np.array([(n.x, n.y) for n in self.balloon_obj.nodes]))
                area_diff = abs(previous_area - current_area) / previous_area
                return area_diff
        return 0

    def delete_event(self):
        if not hasattr(self, "outline_id") or not self.subfigure_patches or self.outline_id is None:
            return

        alert = QtWidgets.QMessageBox()
        delete_confirm = alert.question(
            self.parent(),
            "Delete outline?",
            "Are you really sure you want to delete this outline permanently?"
        )
        if delete_confirm == QtWidgets.QMessageBox.Yes:
            self.set_status("Deleting outline", status="working")
            database.deleteOutlineById(self.outline_id)
            database.updateExperimentById(
                self._data.experiment_id,
                verified=False,
            )
            database.deleteCellById(self.cell_id)

            self.clear_sub_ax()
            self.plot_existing_outlines()
            self.draw()
            self.set_status("Outline deleted")

    def delete_multi(self):
        if len(self.selected_outlines) < 2:
            return

        alert = QtWidgets.QMessageBox()
        delete_confirm = alert.question(
            self.parent(),
            "Delete selected {0} outlines?".format(len(self.selected_outlines)),
            "Are you really sure you want to delete these outlines permanently?"
        )
        if delete_confirm == QtWidgets.QMessageBox.Yes:
            for outline_num, outline in enumerate(self.selected_outlines):
                self.set_status(
                    "Deleting {0: 2d} of {1} outlines".format(
                        outline_num + 1,
                        len(self.selected_outlines),
                    ),
                    status="working",
                )
                database.deleteOutlineById(outline.outline_id)
                database.updateExperimentById(
                    self._data.experiment_id,
                    verified=False,
                )
                database.deleteCellById(outline.cell_id)

            self.set_status("Outlines deleted ({0} outlines)".format(
                len(self.selected_outlines)
            ))
            self.clear_sub_ax()
            self.plot_existing_outlines()
            self.draw()

    def refine_event(self, num_ref=10):
        if not hasattr(self, "outline_id") or not self.subfigure_patches or self.outline_id is None:
            return

        self.set_status("Refining outline", status="working")
        for i in range(num_ref):
            try:
                self.balloon_obj.evolve(image_percentile=self.image_percentile)
            except ValueError:
                QtWidgets.QMessageBox().warning(
                    self.parent(),
                    "I'm sorry, I'm afraid I can't do that",
                    "The outline has shrunk too much, try reducing the tolerance",
                )
                self.set_status(
                    "Refinement failed after {0} cycles".format(i + 1)
                )
                self._plot_nodes()
                self.draw()
                return

        if self.selected_outlines:
            self.selected_outlines[0].is_modified = True

        self.set_status("Refinement complete ({0} cycles)".format(num_ref))
        self._plot_nodes()
        self.draw()

    def refine_multi(self, num_ref=10):
        if len(self.selected_outlines) < 2:
            return

        self.outline_id = None
        for outline_num, outline in enumerate(self.selected_outlines):
            self.set_status(
                "Refining {0: 2d} of {1} outlines".format(
                    outline_num + 1,
                    len(self.selected_outlines),
                ),
                status="working"
            )

            outline_info = database.getOutlineById(outline.outline_id)
            centre = outline_info.centre_y, outline_info.centre_x  # pylint: disable=no-member
            offset_left, offset_top, _, _ = self.get_offsets(centre)
            roi = self.load_frame(channel_idx=0)[
                offset_left:offset_left + (self.region_halfwidth * 2),
                offset_top:offset_top + (self.region_halfheight * 2),
            ]
            if hasattr(outline, "is_modified") and outline.is_modified:
                xy_coords = outline.get_xy()
                current_nodes = np.array([xy_coords[:, 1], xy_coords[:, 0]]).T
            else:
                current_nodes = np.load(outline_info.coords_path)  # pylint: disable=no-member

            balloon_obj = balloon.Balloon(
                current_nodes - np.array([offset_left, offset_top]),
                roi
            )
            balloon_obj.refining_cycles = 1
            for _i in range(num_ref):
                try:
                    balloon_obj.evolve(image_percentile=self.image_percentile)
                except ValueError:
                    break

            coords = np.array([(n.y, n.x) for n in balloon_obj.nodes]) + np.array([
                offset_top,
                offset_left,
            ])
            outline.set_xy(coords)
            outline.is_modified = True
            self.main_ax.draw_artist(outline)

        self.set_status(
            "Refinement complete ({0} outlines)".format(
                len(self.selected_outlines)
            )
        )
        self.draw()


class Toolbar(NavigationToolbar):
    def __init__(self, figure_canvas, experiment_data, parent=None):
        self.toolitems = [
            ("Home", "Home", "home_large", "home_event"),
            (None, None, None, None),
            ("Auto", "Automatic segmentation", "auto_segmentation_1", "auto_segmentation"),
            (None, None, None, None),
            ("Pan", "Pan", "move_large", "pan"),
            ("Zoom", "Zoom", "zoom_to_rect_large", "zoom"),
            (None, None, None, None),
            ("Save", "Save view", "filesave_large", "save_figure"),
            (None, None, None, None),
            ("Accept", "Accept outline (enter)", "accept", "accept"),
            ("Delete", "Delete outline (delete)", "delete", "delete"),
            ("Refine", "Refine outline (r)", "recycle", "refine"),
            ("Refine1", "Refine one step (.)", "recycle_single", "refine_single"),
            (None, None, None, None),
        ]
        if experiment_data.num_channels > 1:
            self.toolitems.extend([
                ("ChannelLeft", "Previous channel (left or a)", "channel_prev", "channel_prev"),
                ("ChannelRight", "Next channel (right or d)", "channel_next", "channel_next"),
            ])

        if experiment_data.num_slices > 1:
            self.toolitems.extend([
                ("SlicePrev", "Previous slice (q)", "slice_prev", "slice_prev"),
                ("SliceNext", "Next slice (e)", "slice_next", "slice_next"),
            ])

        if experiment_data.image_mode == "movie":
            self.toolitems.extend([
                ("FrameUp", "Next frame (up or w)", "frame_next", "frame_next"),
                ("FrameDown", "Previous frame (down or s)", "frame_prev", "frame_prev"),
            ])
        elif experiment_data.image_mode == "static":
            if experiment_data.num_frames > 1:
                self.toolitems.extend([
                    ("ImageUp", "Next image (up or w)", "image_next", "frame_next"),
                    ("ImageDown", "Previous image (down or s)", "image_prev", "frame_prev"),
                ])

        NavigationToolbar.__init__(self, figure_canvas, parent=None)

    def _icon(self, name):
        """Note used by matplotlib toolbar internally"""
        path = os.path.join("resources", name)
        if not os.path.exists(path):
            path = os.path.join(self.basedir, name)

        pixmap = QtGui.QPixmap(path)
        if hasattr(pixmap, "setDevicePixelRatio"):
            # no other way to avoid using _dpi_ratio that I know of
            pixmap.setDevicePixelRatio(self.canvas._dpi_ratio)  # pylint: disable=protected-access

        return QtGui.QIcon(pixmap)

    def home_event(self):
        self.canvas.home_event()

    def auto_segmentation(self):
        startt = time.time()
        self.canvas.automatic_segmentation()
        #  self.canvas.plot_existing_outlines()
        #  self.canvas.refresh()

        print(time.time() - startt)


    def delete(self):
        if len(self.canvas.selected_outlines) > 1:
            self.canvas.delete_multi()
        else:
            self.canvas.delete_event()

    def refine(self):
        if len(self.canvas.selected_outlines) > 1:
            self.canvas.refine_multi()
        else:
            self.canvas.refine_event()

    def refine_single(self):
        if len(self.canvas.selected_outlines) > 1:
            self.canvas.refine_multi(1)
        else:
            self.canvas.refine_event(1)

    def accept(self):
        if len(self.canvas.selected_outlines) > 1:
            self.canvas.accept_multi()
        else:
            self.canvas.accept_event()

    def channel_prev(self):
        self.canvas.channel_change(-1)

    def channel_next(self):
        self.canvas.channel_change(1)

    def frame_next(self):
        self.canvas.frame_change(1)

    def frame_prev(self):
        self.canvas.frame_change(-1)

    def slice_next(self):
        self.canvas.slice_change(1)

    def slice_prev(self):
        self.canvas.slice_change(-1)


class Outliner:
    def __init__(self, experiment_data, image_loader, parent_window):
        self.experiment_data = experiment_data
        self.image_loader = image_loader

        self.plot = None
        self.tolerance_widget = QtWidgets.QLineEdit()
        self.parent_window = parent_window
        self.window = QtWidgets.QDialog(self.parent_window)

        # must be set later
        self.screen_dpi = None
        self.max_width_px = None
        self.max_height_px = None

    def set_screen_res(self, max_width_px, max_height_px, screen_dpi):
        self.max_width_px = max_width_px
        self.max_height_px = max_height_px
        self.screen_dpi = screen_dpi

    def _px_to_in(self, pixels):
        return pixels / self.screen_dpi

    def start_outlining(self):
        if not hasattr(self, "max_width_px"):
            raise ValueError("Screen resolution has not been set")

        dim = (self._px_to_in(self.max_width_px * 0.5),
               self._px_to_in(self.max_height_px * 0.5))

        self.window.setModal(True)
        self.window.setGeometry(
            0, 60,
            self.max_width_px * 0.5,
            self.max_height_px - 120
        )
        self.window.setWindowTitle("Outline cells")

        main_layout = QtWidgets.QVBoxLayout()

        menubar = QtWidgets.QMenuBar(self.window)
        file_menu = menubar.addMenu("&File")
        quit_action = QtWidgets.QAction("&Close", menubar)
        quit_action.triggered.connect(self.window.close)
        file_menu.addAction(quit_action)
        main_layout.setMenuBar(menubar)

        # add tolerance box thing
        label = QtWidgets.QLabel("Tolerance:")
        self.tolerance_widget.setText("1.0")
        # pylint suppression, this lambda is absolutely required for some reason
        self.tolerance_widget.textChanged.connect(lambda text: self._submit_tolerance(text))  # pylint: disable=unnecessary-lambda
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(label)
        layout.addWidget(self.tolerance_widget)

        status_label = QtWidgets.QLabel("Status:")
        status_bar = QtWidgets.QStatusBar()
        status_layout = QtWidgets.QHBoxLayout()
        status_layout.setAlignment(QtCore.Qt.AlignLeft)
        status_layout.addWidget(status_label)
        status_layout.addWidget(status_bar)

        self.plot = Plotter(
            self.window,
            dim[0],
            dim[1],
            dpi=self.screen_dpi,
            experiment_data=self.experiment_data,
            image_loader=self.image_loader,
            status_bar=status_bar,
        )
        self.plot.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.plot.setFocus()
        self.window.toolbar = Toolbar(self.plot, self.experiment_data, self.window)

        tool_layout = QtWidgets.QVBoxLayout()
        tool_layout.addWidget(self.window.toolbar)
        tool_layout.addLayout(layout)
        tool_layout.addLayout(status_layout)

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
