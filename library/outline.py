#!/usr/bin/env python3

""" Module that handles outlining (automated and manually supervised)

The interface opens with a toolbar containing standard matplotlib icons: home,
pan, zoom, save; and custom controls for PombeTrack: automatic segmentation,
accept, delete, refine, change channel, change frame, and change slice.

Underneath the toolbar is an input box for the "tolerance" parameter, which
refers to the expandability of the balloon (see library.balloon).

Below this is a status indicator which will contain text when certain
operations are performed.

The main interface consists of two matplotlib panels, a "main" panel on the
left which displays the whole field-of-view of a particular frame, and a right
"sub" panel, which is used to display and control the outlining procedure.
See the Plotter class below for how these panels are controlled.
"""

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

class DragRect:
    """ Convenience class for controlling a draggable select box on the main axis.

    Arguments:
        bound_fig (matplotlib.figure.Figure): the figure containing the axis.
        target_ax (matplotlib.axes.Axes): the axis in question.
    """
    INACTIVE = 0
    ACTIVE = 1
    def __init__(self, bound_fig, target_ax):
        self.fig = bound_fig
        self.axis = target_ax
        self.rect = self.init_rect()
        self.status = self.INACTIVE
        self.background = self.fig.canvas.copy_from_bbox(
            self.axis.bbox
        )

    @staticmethod
    def init_rect():
        """ Create a placeholder rectangle object.

        Arguments: N/A
        Returns:
            rect (matplotlib.patches.Rectangle):
                Placeholder with zero width and height for the draggable
                rectangle.

        The rectangle is initialised with zero width and height, with a yellow
        border colour and a partially transparent gray fill colour.
        """
        rect = matplotlib.patches.Rectangle(
            (0, 0),
            0, 0,
            edgecolor="yellow",
            facecolor=(0.9, 0.9, 0.7, 0.3),
        )
        return rect

    def place(self, pos_x, pos_y):
        """ Set the bottom left corner of the rectangle, and draw it.

        Arguments:
            pos_x (float): x coordinate.
            pos_y (float): y coordinate.

        Returns:
            None

        If the rectangle has not yet been drawn to the axis, it will be added,
        and the current state of the axis (pre-drawing) saved.
        This saving of the "background" permits smooth animations when the
        rectangle is modified.
        """
        self.rect.set_xy((pos_x, pos_y))
        if self.status == self.INACTIVE:
            self.axis.add_patch(self.rect)
            self.background = self.fig.canvas.copy_from_bbox(
                self.axis.bbox
            )

        self.status = self.ACTIVE

        self.fig.canvas.draw()

    def resize(self, evt):
        """ Change the size of the rectangle according to a mouse event.

        Arguments:
            evt (matplotlib.backend_bases.MouseEvent):
                Mouse event triggering the resizing.

        Returns:
            None

        Determines what size the rectangle should be according to the position
        of the mouse cursor, then animates the new rectangle smoothly (using
        blitting).
        """
        self.fig.canvas.restore_region(self.background)
        width, height = self.get_dimensions(evt)
        self.rect.set_width(width)
        self.rect.set_height(height)
        self.axis.draw_artist(self.rect)
        self.fig.canvas.blit(self.axis.bbox)

    def get_dimensions(self, evt=None):
        """ Calculate the size of the rectangle.

        Arguments:
            evt (matplotlib.backend_bases.MouseEvent):
                Mouse event that triggered a resizing (see DragRect.resize).
                Optional argument, will return the current size if omitted.

        Returns:
            rect_width (float): Width in pixels of the rectangle.
            rect_height (float): Height in pixels of the rectangle.

        If the evt argument is passed, the width/height will be calculated
        according to the distance from the mouse cursor to the rectangle bottom
        left corner (negative if necessary).
        """

        if not evt or not evt.xdata or not evt.ydata or evt.inaxes != self.axis:
            rect_width = self.rect.get_width()
            rect_height = self.rect.get_height()
        else:
            rect_width = evt.xdata - self.rect.get_x()
            rect_height = evt.ydata - self.rect.get_y()

        return rect_width, rect_height

    def get_path(self, transform=True):
        """ Get path for the rectangle.

        Arguments:
            transform (bool):
                Whether to transform the path into axis coordinates.
                Defaults to True.

        Returns:
            rect_path (matplotlib.path.Path):
                The Path describing the rectangle.
        """
        rect_path = self.rect.get_path()
        if transform:
            rect_transform = self.rect.get_transform()
            true_rect = rect_transform.transform_path(rect_path)
            return true_rect

        return rect_path

    def reset(self):
        """ Remove the rectangle from the axis.

        Arguments: N/A
        Returns: None

        Only removes the rectangle if one exists.
        """
        if self.is_active():
            self.status = self.INACTIVE
            self.rect.remove()
            self.rect = self.init_rect()
            self.fig.canvas.draw()

    def is_active(self):
        """ Get rectangle status.

        Arguments: N/A
        Returns:
            active (bool): Rectangle status (True if drawn, False if not).
        """
        return self.status == self.ACTIVE


class Plotter(FigureCanvas):
    def __init__(self, parent_window, dims, experiment_data, status_bar):
        self.fig = matplotlib.figure.Figure(
            figsize=(dims[0], dims[1]),
            dpi=dims[2],
        )
        self.fig.frame_idx = 0
        self.fig.channel_idx = 0
        self.fig.slice_idx = 0

        self.main_ax = self.fig.add_subplot(121)
        self.sub_ax = self.fig.add_subplot(122)

        self.decorate_axis(self.main_ax)
        self.decorate_axis(self.sub_ax)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent_window)

        FigureCanvas.setSizePolicy(
            self,
            QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding,
        )
        FigureCanvas.updateGeometry(self)

        self.experiment = experiment_data
        self.image_percentile = 1.0

        if not database.checkTable("outlines"):
            database.createOutlinesTable()

        if not database.checkTable("cells"):
            database.createCellsTable()

        self.load_metadata()

        self.main_ax.drag = DragRect(self.fig, self.main_ax)
        self.main_ax.image = self.main_ax.imshow(self.load_frame(), cmap="gray")
        self.main_ax.outlines = []
        self.main_ax.selected_outlines = []

        self.sub_ax.outlines = []
        self.sub_ax.nodes = []
        self.sub_ax.is_dragging = False
        # dummy patch for drag node (is replaced when a node is dragged)
        self.sub_ax.drag_node = matplotlib.patches.Circle((0, 0), 1.5)
        self.sub_ax.background = self.fig.canvas.copy_from_bbox(
            self.sub_ax.bbox
        )
        self.sub_ax.previous_id = None
        self.sub_ax.cell_id = None
        self.sub_ax.outline_id = None
        self.sub_ax.balloon = None
        self.sub_ax.centre_x = None
        self.sub_ax.centre_y = None
        self.sub_ax.offset_left = None
        self.sub_ax.offset_top = None
        self.sub_ax.region_halfwidth = 75
        self.sub_ax.region_halfheight = 75

        self.mpl_connect("key_press_event", self._key_press_event)
        self.mpl_connect("button_press_event", self._button_press_event)
        self.mpl_connect("button_release_event", self._button_release_event)
        self.mpl_connect("motion_notify_event", self._motion_notify_event)

        self.status_bar = status_bar
        self.current_status = None

        self.plot_existing_outlines()

    def count_selected(self):
        return len(self.main_ax.selected_outlines)

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

    def preprocess(self):
        self.set_status(
            "Preprocessing",
            status="working",
        )
        # load_frame: frame, z-slice, channel
        img_mid = self.load_frame(
            self.fig.frame_idx,
            int(np.floor(self.num_slices / 2)),
            0
        )
        try:
            img_up = self.load_frame(
                self.fig.frame_idx,
                int(np.floor(self.num_slices / 2) - 1),
                0
            )
        except:
            self.set_status(text="Preprocessing failed since not enough z-stacks are provided")
            return
        img = np.maximum(img_mid, img_up)

        img_pp = segmentation.preprocessing(img)
        img_i = segmentation.find_cellinterior(img_pp)
        img_wat = segmentation.find_watershed(img_i)
        #  bd = segmentation.find_bd(img_wat)
        return img, img_wat

    def automatic_segmentation(self, display=True):
        img, img_wat = self.preprocess()
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
            balloon_results = segmentation.find_balloon_obj(
                bd_ii_sorted.astype(int)[::5],
                img,
            )
            self.sub_ax.balloon = balloon_results[0]
            self.sub_ax.offset_left = balloon_results[1]
            self.sub_ax.offset_top = balloon_results[2]
            halfwidth = balloon_results[3]

            #  Test if the cell exists
            overlap = False
            outline_data = database.getOutlinesByFrameIdx(
                self.fig.frame_idx,
                self.experiment.experiment_id
            )
            for outline in outline_data:
                if not os.path.exists(outline.coords_path):
                    continue

                test_range_x = range(self.sub_ax.offset_top,
                                     self.sub_ax.offset_top + 2 * halfwidth)
                test_range_y = range(self.sub_ax.offset_left,
                                     self.sub_ax.offset_left + 2 * halfwidth)
                if (outline.centre_x not in test_range_x or
                        outline.centre_y not in test_range_y):
                    continue

                polygonpath = matplotlib.path.Path(
                    np.append(self.sub_ax.balloon.get_coordinates(),
                              self.sub_ax.balloon.get_coordinates()[1, :].reshape(1, 2),
                              axis=0),
                    closed=True
                )
                if polygonpath.contains_point([outline.centre_y - self.sub_ax.offset_left,
                                               outline.centre_x - self.sub_ax.offset_top]):
                    overlap = True

            if overlap:
                continue

            # Evolve the contour
            try:
                sensitivity = self.image_percentile
                area_init = self.sub_ax.balloon.get_area()
                for _ in range(20):
                    self.sub_ax.balloon.evolve(image_percentile=sensitivity)
                    if (self.sub_ax.balloon.get_area() > 1.5 * area_init or
                            self.sub_ax.balloon.get_area() < 0.5 * area_init):
                        raise ValueError

                full_coords = self.sub_ax.balloon.get_coordinates() + [
                    self.sub_ax.offset_left, self.sub_ax.offset_top
                ]
                self.sub_ax.outline_id = str(uuid.uuid4())
                self.sub_ax.cell_id = str(uuid.uuid4())
                self.sub_ax.centre_y, self.sub_ax.centre_x = full_coords.mean(axis=0).astype(int)
                self.save_outline()

                # Draw the cell
                if display:
                    background = self.draw_outline(background, full_coords)
            except ValueError:
                continue

            cell_count += 1
            self.draw()

        self.set_status(text="Auto segmentation finished with {0} cells".format(cell_count))

    def draw_outline(self, background, full_coords):
        self.fig.canvas.restore_region(background)
        outline_poly = matplotlib.patches.Polygon(
            np.array([full_coords[:, 1],
                      full_coords[:, 0]]).T,
            edgecolor="r",
            fill=False,
            lw=1
        )
        outline_poly.outline_id = self.sub_ax.outline_id
        self.main_ax.add_patch(outline_poly)
        self.main_ax.outlines.append(outline_poly)
        self.main_ax.draw_artist(outline_poly)
        self.fig.canvas.blit(self.main_ax.bbox)
        background = self.fig.canvas.copy_from_bbox(self.main_ax.bbox)
        return background

    def load_metadata(self):
        self.num_frames = self.experiment.image_loader.num_frames
        self.num_channels = self.experiment.image_loader.num_channels
        self.num_slices = self.experiment.image_loader.num_slices

    def load_frame(self, frame_idx=None, slice_idx=None, channel_idx=None):
        if frame_idx is None:
            frame_idx = self.fig.frame_idx

        if channel_idx is None:
            channel_idx = self.fig.channel_idx

        if slice_idx is None:
            slice_idx = self.fig.slice_idx

        if frame_idx < 0 or frame_idx > (self.num_frames - 1):
            return np.zeros((100, 100))

        return self.experiment.image_loader.load_frame(frame_idx, slice_idx, channel_idx)

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
        self.sub_ax.outline_id = None
        self.sub_ax.balloon = None
        self.main_ax.drag.reset()
        self.sub_ax.is_dragging = False
        self.sub_ax.nodes = []

    def clear_sub_outlines(self):
        while True:
            try:
                outline_poly = self.sub_ax.outlines.pop()
                outline_poly.remove()
            except ValueError:
                pass
            except IndexError:
                break

    def plot_existing_outlines(self):
        if self.experiment.image_mode == "static":
            title = "Image #{0}".format(self.fig.frame_idx + 1)
        else:
            title = "F={0}".format(self.fig.frame_idx + 1)

        if self.experiment.num_slices > 1:
            title += " Z={0}".format(self.fig.slice_idx + 1)

        if self.experiment.num_channels > 1:
            title += " C={0}".format(self.fig.channel_idx + 1)

        self.main_ax.set_title(title)
        while True:
            try:
                self.main_ax.lines.pop()
            except IndexError:
                break

        while True:
            try:
                outline_poly = self.main_ax.outlines.pop()
                outline_poly.remove()
            except IndexError:
                break

        selected_ids = [x.outline_id for x in self.main_ax.selected_outlines]
        fresh_selections = []
        outline_data = database.getOutlinesByFrameIdx(
            self.fig.frame_idx,
            self.experiment.experiment_id,
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
                prev_outline = self.main_ax.selected_outlines[
                    selected_ids.index(outline.outline_id)
                ]
                if hasattr(prev_outline, "is_modified") and prev_outline.is_modified:
                    outline_poly.set_xy(prev_outline.get_xy())
                    outline_poly.is_modified = True
                outline_poly.set_edgecolor("yellow")
                outline_poly.is_selected = True
                fresh_selections.append(outline_poly)

            self.main_ax.add_patch(outline_poly)
            self.main_ax.outlines.append(outline_poly)

        self.main_ax.selected_outlines = fresh_selections

    def save_outline(self, explicit=None):
        if explicit:
            self.sub_ax.outline_id = explicit.outline_id
            outline_info = database.getOutlineById(self.sub_ax.outline_id)
            xy_coords = explicit.get_xy()
            coords = np.array([xy_coords[:, 1], xy_coords[:, 0]]).T
            self.sub_ax.cell_id = explicit.cell_id
            self.sub_ax.previous_id = outline_info.parent_id  # pylint: disable=no-member
            self.sub_ax.centre_y, self.sub_ax.centre_x = coords.mean(axis=0)
        else:
            coords_offset = np.array([(n.x, n.y) for n in self.sub_ax.balloon.nodes])
            # Determine cell centre
            coords = coords_offset + np.array([self.sub_ax.offset_left, self.sub_ax.offset_top])
            self.sub_ax.centre_y, self.sub_ax.centre_x = coords.mean(axis=0)

        outline_store = os.path.join(
            "data", "outlines", self.experiment.experiment_id
        )
        if not os.path.exists(outline_store):
            os.makedirs(outline_store)

        coords_path = os.path.join(
            outline_store,
            "{0}.npy".format(self.sub_ax.outline_id)
        )
        data = {
            "outline_id": self.sub_ax.outline_id,
            "cell_id": self.sub_ax.cell_id,
            "experiment_num": self.experiment.experiment_num,
            "experiment_id": self.experiment.experiment_id,
            "image_path": self.experiment.image_path,
            "frame_idx": self.fig.frame_idx,
            "coords_path": coords_path,
            "parent_id": self.sub_ax.previous_id or "",
            "centre_x":int(self.sub_ax.centre_x),
            "centre_y":int(self.sub_ax.centre_y),
        }
        if self.experiment.image_mode == "static":
            data["parent_id"] = ""

        if os.path.exists(coords_path):
            os.remove(coords_path)
            database.updateOutlineById(
                self.sub_ax.outline_id,
                centre_x=int(self.sub_ax.centre_x),
                centre_y=int(self.sub_ax.centre_y),
            )
        else:
            database.insertOutline(**data)

        np.save(coords_path, coords)

        if self.sub_ax.previous_id:
            database.addOutlineChild(self.sub_ax.previous_id, child1=self.sub_ax.outline_id)

        database.updateExperimentById(
            self.experiment.experiment_id,
            verified=False,
        )
        database.deleteCellById(self.sub_ax.cell_id)

    def fit_outline(self, roi, init_nodes=None, centre_offset_left=0, centre_offset_top=0):
        if init_nodes is None:
            centre = [
                self.sub_ax.region_halfwidth - centre_offset_left,
                self.sub_ax.region_halfheight - centre_offset_top,
            ]
            self.sub_ax.outline_id = str(uuid.uuid4())
            radius = 5
            num_nodes = 10
            init_nodes = balloon.initial_nodes(centre, radius, num_nodes)

        self.sub_ax.balloon = balloon.Balloon(init_nodes, roi)
        self.sub_ax.imshow(self.load_frame(), cmap="gray")
        self.sub_ax.set_ylim([
            self.sub_ax.offset_left + self.sub_ax.region_halfwidth * 2,
            self.sub_ax.offset_left
        ])
        self.sub_ax.set_xlim([
            self.sub_ax.offset_top,
            self.sub_ax.region_halfheight * 2 + self.sub_ax.offset_top
        ])
        self.sub_ax.set_title("Frame = {0}".format(self.fig.frame_idx + 1))
        self._plot_nodes()

        self.clear_sub_outlines()
        outline_data = database.getOutlinesByFrameIdx(
            self.fig.frame_idx,
            self.experiment.experiment_id,
        )
        for outline in outline_data:
            if (not os.path.exists(outline.coords_path) or
                    outline.outline_id == self.sub_ax.outline_id):
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
            self.sub_ax.outlines.append(poly)

        self.draw()

    def _plot_nodes(self, line=True, patches=True):
        nodes = self.sub_ax.balloon.nodes
        if line:
            try:
                self.sub_ax.lines.pop(0)
            except IndexError:
                pass
            coords = np.array([
                (n.x + self.sub_ax.offset_left,
                 n.y + self.sub_ax.offset_top)
                for n in nodes
            ])
            self.sub_ax.plot(
                coords[:, 1],
                coords[:, 0],
                color="y",
                lw=1,
            )

        if patches:
            for _ in range(len(self.sub_ax.nodes)):
                patch = self.sub_ax.nodes.pop()
                patch.remove()

            for node_idx, node in enumerate(nodes):
                patch = matplotlib.patches.Circle(
                    (node.y + self.sub_ax.offset_top,
                     node.x + self.sub_ax.offset_left),
                    1.5,
                    fc="y",
                )
                patch.this_node = node
                patch.node_idx = node_idx
                self.sub_ax.nodes.append(patch)
                self.sub_ax.add_artist(patch)

        self.draw()

    def slice_change(self, delta):
        if delta < 0 and self.fig.slice_idx <= 0:
            return

        if delta > 0 and self.fig.slice_idx >= self.num_slices - 1:
            return

        self.fig.slice_idx += delta
        new_im = self.load_frame()
        self.main_ax.image.set_data(new_im)
        self.plot_existing_outlines()
        self.main_ax.image.set_clim([new_im.min(), new_im.max()])
        self.draw()

    def channel_change(self, delta):
        if delta < 0 and self.fig.channel_idx <= 0:
            return

        if delta > 0 and self.fig.channel_idx >= self.num_channels - 1:
            return

        self.fig.channel_idx += delta
        new_im = self.load_frame()
        self.main_ax.image.set_data(new_im)
        self.plot_existing_outlines()
        self.main_ax.image.set_clim([new_im.min(), new_im.max()])
        self.draw()

    def frame_change(self, delta):
        if delta < 0 and self.fig.frame_idx <= 0:
            return

        if delta > 0 and self.fig.frame_idx >= self.num_frames - 1:
            return

        if not self.check_selected_outlines():
            return

        self.fig.frame_idx += delta
        new_im = self.load_frame()
        self.main_ax.image.set_data(new_im)
        self.deselect_outlines()
        self.plot_existing_outlines()
        self.main_ax.image.set_clim([new_im.min(), new_im.max()])
        self.draw()

    def _key_press_event(self, evt):
        key_events = {
            "left": (self.channel_change, -1),
            "a": (self.channel_change, -1),
            "right": (self.channel_change, 1),
            "d": (self.channel_change, 1),
            "up": (self.frame_change, 1),
            "w": (self.frame_change, 1),
            "down": (self.frame_change, -1),
            "s": (self.frame_change, -1),
            "q": (self.slice_change, -1),
            "e": (self.slice_change, 1),
            "r": (self.refine,),
            ".": (self.refine, 1),
            "delete": (self.delete,),
            "enter": (self.accept,),
            "shift": (lambda: None,),
            "control": (lambda: None,),
        }

        if evt.key in key_events:
            func, *args = key_events[evt.key]
            func(*args)
        else:
            print("Unknown key:", evt.key)

    def get_offsets(self, centre):
        offset_left = int(round(centre[0] - self.sub_ax.region_halfwidth))
        offset_top = int(round(centre[1] - self.sub_ax.region_halfheight))
        centre_offset_left, centre_offset_top = 0, 0
        img = self.load_frame()
        if offset_left < 0:
            centre_offset_left = -offset_left
            offset_left = 0
        elif offset_left >= img.shape[0] - (self.sub_ax.region_halfwidth * 2):
            centre_offset_left = img.shape[0] - (self.sub_ax.region_halfwidth * 2) - offset_left
            offset_left = img.shape[0] - (self.sub_ax.region_halfwidth * 2)

        if offset_top < 0:
            centre_offset_top = -offset_top
            offset_top = 0
        elif offset_top >= img.shape[1] - (self.sub_ax.region_halfheight * 2):
            centre_offset_top = img.shape[1] - (self.sub_ax.region_halfheight * 2) - offset_top
            offset_top = img.shape[1] - (self.sub_ax.region_halfheight * 2)
        del img

        return offset_left, offset_top, centre_offset_left, centre_offset_top

    def _button_press_event(self, evt):
        if self.parent().toolbar.mode:
            return

        if evt.inaxes == self.main_ax:
            self.main_ax.drag.place(evt.xdata, evt.ydata)

        elif evt.inaxes == self.sub_ax:
            if evt.button == 1:
                for patch in self.sub_ax.nodes:
                    if patch.contains_point((evt.x, evt.y)):
                        patch.set_visible(False)
                        self.draw()
                        self.sub_ax.background = self.fig.canvas.copy_from_bbox(
                            self.sub_ax.bbox
                        )
                        self.sub_ax.lines[0].set_color("red")
                        patch.set_visible(True)
                        self.sub_ax.is_dragging = True
                        self.sub_ax.drag_node = patch
                        patch.set_facecolor("r")
                        # self.draw()
                        break

            elif evt.button == 3:
                for patch in self.sub_ax.nodes:
                    if patch.contains_point((evt.x, evt.y)):
                        self.sub_ax.balloon.remove_node(patch.this_node)
                        self._plot_nodes()
                        self.draw()

    def check_selected_outlines(self):
        if sum([hasattr(x, "is_modified") and x.is_modified or False
                for x in self.main_ax.selected_outlines]) > 0:
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
                for outline_num, outline in enumerate(self.main_ax.selected_outlines):
                    self.set_status(
                        "Saving {0} of {1} outlines".format(
                            outline_num + 1,
                            len(self.main_ax.selected_outlines),
                        )
                    )
                    self.save_outline(explicit=outline)
                    outline.is_modified = False
                self.set_status(
                    "Modified outlines saved ({0} outlines)".format(
                        len(self.main_ax.selected_outlines),
                    )
                )
                self.plot_existing_outlines()
                self.draw()
                return True

            if add_confirm == QtWidgets.QMessageBox.No:
                self.set_status(
                    "Modified outlines discarded ({0} outlines)".format(
                        len(self.main_ax.selected_outlines),
                    )
                )

                for outline in self.main_ax.selected_outlines:
                    outline.is_modified = False

                self.plot_existing_outlines()
                self.draw()
                return True

            if add_confirm == QtWidgets.QMessageBox.Cancel:
                self.set_status(clear=True)
                return False

        return True

    def deselect_outlines(self):
        for outline in self.main_ax.selected_outlines:
            outline.set_edgecolor("red")
            outline.is_selected = False
            outline.is_modified = False

        self.main_ax.selected_outlines = []

    def _button_release_main(self, evt):
        rect_width, rect_height = self.main_ax.drag.get_dimensions(evt)
        keyboard_mods = QtGui.QGuiApplication.queryKeyboardModifiers()
        # if keyboard_mods & (QtCore.Qt.ShiftModifier | QtCore.Qt.ControlModifier):
        if keyboard_mods & QtCore.Qt.ShiftModifier:
            # additive mode
            additive_flag = True
            reductive_flag = False
            self.select_new_outlines(check=False)

        elif keyboard_mods & QtCore.Qt.ControlModifier:
            additive_flag = False
            reductive_flag = True
            check = self.deselect_new_outlines()
            if not check:
                return

        else:
            additive_flag = False
            reductive_flag = False
            check = self.select_new_outlines(check=True, deselect=True)
            if not check:
                return

        self.main_ax.drag.reset()

        if len(self.main_ax.selected_outlines) > 1:
            self.clear_sub_ax()

        elif len(self.main_ax.selected_outlines) == 1:
            self._select_hit(self.main_ax.selected_outlines[0])

        elif not self.main_ax.selected_outlines:
            self.clear_sub_ax()
            if (not additive_flag and not reductive_flag and
                    rect_width < 10 and rect_height < 10):
                self.sub_ax.previous_id = None
                self.sub_ax.cell_id = str(uuid.uuid4())
                centre = [evt.ydata, evt.xdata]
                (self.sub_ax.offset_left, self.sub_ax.offset_top,
                 centre_offset_left, centre_offset_top) = self.get_offsets(centre)

                right = self.sub_ax.offset_left + (self.sub_ax.region_halfwidth * 2)
                bottom = self.sub_ax.offset_top + (self.sub_ax.region_halfheight * 2)
                roi = self.load_frame(channel_idx=0)[
                    self.sub_ax.offset_left:right,
                    self.sub_ax.offset_top:bottom,
                ]

                self.fit_outline(
                    roi,
                    centre_offset_left=centre_offset_left,
                    centre_offset_top=centre_offset_top,
                )

        self.draw()


    def select_new_outlines(self, check=True, deselect=False):
        if check:
            response = self.check_selected_outlines()
            if not response:
                self.main_ax.drag.reset()
                return False

        if deselect:
            self.deselect_outlines()

        new_outlines = self.get_outlines_in_rect()
        for outline in new_outlines:
            outline.set_edgecolor("yellow")
            outline.is_selected = True
            self.main_ax.selected_outlines.append(outline)

        return True

    def deselect_new_outlines(self, check=True):
        if check:
            response = self.check_selected_outlines()
            if not response:
                self.main_ax.drag.reset()
                return False

        new_outlines = self.get_outlines_in_rect()
        for outline in new_outlines:
            if hasattr(outline, "is_selected") and outline.is_selected:
                outline.set_edgecolor("red")
                outline.is_selected = False
                self.main_ax.selected_outlines.pop(
                    self.main_ax.selected_outlines.index(outline)
                )

        return True


    def get_outlines_in_rect(self, evt=None):
        rect_path = self.main_ax.drag.get_path()
        new_outlines = []
        for outline in self.main_ax.outlines:
            out_transform = outline.get_transform()
            out_path = outline.get_path()
            true_outline = out_transform.transform_path(out_path)
            contains = rect_path.contains_points(
                true_outline.vertices
            )
            if sum(contains) > 0:
                new_outlines.append(outline)

        if not new_outlines:
            if not evt:
                evt = rect_path.vertices.mean(axis=0)
                hit = self._check_hit(evt, mouse=False)
            else:
                hit = self._check_hit(evt)

            if hit:
                new_outlines.append(hit)

        return new_outlines

    def _button_release_event(self, evt):
        if self.parent().toolbar.mode:
            return

        if self.main_ax.drag.is_active():
            self._button_release_main(evt)
            return

        if evt.inaxes == self.sub_ax:
            if not self.sub_ax.is_dragging:
                return

            self.sub_ax.drag_node.set_facecolor("y")
            self.sub_ax.drag_node.this_node.set_position((
                evt.ydata - self.sub_ax.offset_left,
                evt.xdata - self.sub_ax.offset_top,
            ))
            self.sub_ax.drag_node.this_node.apply_changes()
            self._plot_nodes()
            self.sub_ax.is_dragging = False
            return

    def _check_hit(self, evt, mouse=True):
        for outline in self.main_ax.outlines:
            if mouse:
                hit, _ = outline.contains(evt)
            else:
                hit = outline.contains_point(evt)

            if hit:
                return outline
        return False

    def _select_hit(self, outline):
        outline.set_edgecolor("yellow")
        outline.is_selected = True
        outline_info = database.getOutlineById(outline.outline_id)
        self.sub_ax.previous_id = outline_info.parent_id  # pylint: disable=no-member
        self.sub_ax.cell_id = outline_info.cell_id  # pylint: disable=no-member
        centre = outline_info.centre_y, outline_info.centre_x  # pylint: disable=no-member
        (self.sub_ax.offset_left, self.sub_ax.offset_top,
         centre_offset_left, centre_offset_top) = self.get_offsets(centre)
        roi = self.load_frame(channel_idx=0)[
            self.sub_ax.offset_left:self.sub_ax.offset_left + (self.sub_ax.region_halfwidth * 2),
            self.sub_ax.offset_top:self.sub_ax.offset_top + (self.sub_ax.region_halfheight * 2),
        ]
        self.sub_ax.outline_id = outline_info.outline_id  # pylint: disable=no-member
        current_nodes = np.load(outline_info.coords_path)  # pylint: disable=no-member
        self.fit_outline(
            roi,
            init_nodes=current_nodes - np.array([self.sub_ax.offset_left, self.sub_ax.offset_top]),
            centre_offset_left=centre_offset_left,
            centre_offset_top=centre_offset_top,
        )
        self.sub_ax.balloon.refining_cycles = 1

    def _motion_notify_event(self, evt):
        if self.parent().toolbar.mode:
            return

        if self.main_ax.drag.is_active() and evt.inaxes == self.main_ax:
            self.main_ax.drag.resize(evt)

        elif self.sub_ax.is_dragging and evt.inaxes == self.sub_ax:
            self.fig.canvas.restore_region(self.sub_ax.background)
            self.sub_ax.drag_node.center = evt.xdata, evt.ydata
            xdata, ydata = self.sub_ax.lines[0].get_data()
            xdata[self.sub_ax.drag_node.node_idx] = evt.xdata
            ydata[self.sub_ax.drag_node.node_idx] = evt.ydata
            self.sub_ax.lines[0].set_data(xdata, ydata)
            self.sub_ax.draw_artist(self.sub_ax.lines[0])
            self.sub_ax.draw_artist(self.sub_ax.drag_node)
            self.fig.canvas.blit(self.sub_ax.bbox)

    def home_event(self):
        frame = self.load_frame()
        self.main_ax.set_xlim([0, frame.shape[0]])
        self.main_ax.set_ylim([frame.shape[1], 0])
        self.sub_ax.set_xlim([0, self.sub_ax.region_halfwidth * 2])
        self.sub_ax.set_ylim([self.sub_ax.region_halfheight * 2, 0])
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

    def accept(self):
        if len(self.main_ax.selected_outlines) > 1:
            self.accept_multi()
        elif self.sub_ax.nodes:
            self.accept_event()

    def accept_event(self):
        if not self.sub_ax.nodes or self.sub_ax.outline_id is None:
            return

        self.set_status("Saving outline", status="working")
        # check balloon object has been refined at least once
        if self.sub_ax.balloon.refining_cycles == 0:
            add_confirm = self.question(
                ["The outline you are about to add has not been refined.",
                 "Are you sure you want to add it?"],
                "Add outline?",
            )
            if add_confirm != QtWidgets.QMessageBox.Yes:
                self.set_status("Save aborted")
                return

        # check overlap with other outlines?
        # objects already exist (red lines) so can compare perhaps?
        x_points, y_points = self.sub_ax.lines[0].get_data()
        existing_points = np.row_stack([x_points, y_points]).T
        for outline in self.sub_ax.outlines:
            contains_new = sum(outline.get_path().contains_points(
                existing_points,
            ))
            if contains_new > 0:
                add_confirm = self.question(
                    ["The outline you are about to add may overlap with existing outlines.",
                     "Are you sure you would like to add it?"],
                    "Add outline?",
                )
                if add_confirm != QtWidgets.QMessageBox.Yes:
                    self.set_status("Save aborted")
                    return

        if self.experiment.image_mode == "movie":
            self.accept_event_movie()

        elif self.experiment.image_mode == "static":
            self.save_outline()
            self.clear_sub_ax()

            # update existing outlines
            self.deselect_outlines()
            self.plot_existing_outlines()
            self.draw()
            self.set_status("Outline saved")

    def accept_multi(self):
        for outline_num, outline in enumerate(self.main_ax.selected_outlines):
            self.set_status(
                "Saving {0: 2d} of {1} outlines".format(
                    outline_num + 1,
                    len(self.main_ax.selected_outlines),
                ),
                status="working"
            )
            self.save_outline(explicit=outline)

        self.set_status(
            "Outlines saved ({0} outlines)".format(
                len(self.main_ax.selected_outlines)
            )
        )

        self.deselect_outlines()
        self.draw()

    def question(self, message_lines, title):
        message = "\n".join(message_lines)
        alert = QtWidgets.QMessageBox()
        response = alert.question(self.parent(), title, message)
        return response

    def accept_event_movie(self):
        area_diff = self._area_diff()
        if area_diff > 0.3:
            add_confirm = self.question(
                ["The outline you are about to add has a large difference "
                 "in cell area compared to its predecessor in the previous "
                 "frame ({0:.0f}%).".format(area_diff * 100),
                 "Are you sure you want to add it?"],
                "Add outline?",
            )
            if add_confirm != QtWidgets.QMessageBox.Yes:
                self.set_status("Save aborted")
                return

        self.save_outline()
        self.deselect_outlines()

        offset_centre = self.sub_ax.balloon.get_centre()

        # clear plot
        self.sub_ax.previous_id = str(self.sub_ax.outline_id)
        self.clear_sub_ax()

        # fit next
        centre = [offset_centre[0] + self.sub_ax.offset_left,
                  offset_centre[1] + self.sub_ax.offset_top]

        (self.sub_ax.offset_left, self.sub_ax.offset_top,
         centre_offset_left, centre_offset_top) = self.get_offsets(centre)

        if self.fig.frame_idx == self.num_frames - 1:
            bf_frame = self.load_frame()
            self.main_ax.image.set_data(bf_frame)
            self.sub_ax.previous_id = None
            self.sub_ax.outline_id = None
            self.plot_existing_outlines()
            self.clear_sub_outlines()
            self.draw()
            return

        self.fig.frame_idx += 1
        bf_frame = self.load_frame()
        self.main_ax.image.set_data(bf_frame)
        self.plot_existing_outlines()
        self.clear_sub_outlines()

        roi = self.load_frame(channel_idx=0)[
            self.sub_ax.offset_left:self.sub_ax.offset_left + (self.sub_ax.region_halfwidth * 2),
            self.sub_ax.offset_top:self.sub_ax.offset_top + (self.sub_ax.region_halfheight * 2)
        ]
        self.fit_outline(
            roi,
            centre_offset_left=centre_offset_left,
            centre_offset_top=centre_offset_top,
        )
        self.set_status(clear=True)

    def _area_diff(self):
        # check difference in total area is small
        if self.sub_ax.previous_id:
            previous_outline_entry = database.getOutlineById(self.sub_ax.previous_id)
            previous_outline = np.load(previous_outline_entry.coords_path)  # pylint: disable=no-member
            if previous_outline_entry.cell_id == self.sub_ax.cell_id:  # pylint: disable=no-member
                previous_area = self.get_area(previous_outline)
                current_area = self.get_area(np.array([
                    (n.x, n.y)
                    for n in self.sub_ax.balloon.nodes
                ]))
                area_diff = abs(previous_area - current_area) / previous_area
                return area_diff
        return 0

    def delete(self):
        if len(self.main_ax.selected_outlines) > 1:
            self.delete_multi()
        elif self.sub_ax.nodes:
            self.delete_event()

    def delete_event(self):
        if not self.sub_ax.nodes or self.sub_ax.outline_id is None:
            return

        delete_confirm = self.question(
            ["Are you really sure you want to delete this outline permanently?"],
            "Delete outline?",
        )
        if delete_confirm == QtWidgets.QMessageBox.Yes:
            self.set_status("Deleting outline", status="working")
            database.deleteOutlineById(self.sub_ax.outline_id)
            database.updateExperimentById(
                self.experiment.experiment_id,
                verified=False,
            )
            database.deleteCellById(self.sub_ax.cell_id)

            self.clear_sub_ax()
            self.plot_existing_outlines()
            self.draw()
            self.set_status("Outline deleted")

    def delete_multi(self):
        if len(self.main_ax.selected_outlines) < 2:
            return

        delete_confirm = self.question(
            ["Are you really sure you want to delete these outlines permanently?"],
            "Delete selected {0} outlines?".format(len(self.main_ax.selected_outlines)),
        )
        if delete_confirm == QtWidgets.QMessageBox.Yes:
            for outline_num, outline in enumerate(self.main_ax.selected_outlines):
                self.set_status(
                    "Deleting {0: 2d} of {1} outlines".format(
                        outline_num + 1,
                        len(self.main_ax.selected_outlines),
                    ),
                    status="working",
                )
                database.deleteOutlineById(outline.outline_id)
                database.updateExperimentById(
                    self.experiment.experiment_id,
                    verified=False,
                )
                database.deleteCellById(outline.cell_id)

            self.set_status("Outlines deleted ({0} outlines)".format(
                len(self.main_ax.selected_outlines)
            ))
            self.clear_sub_ax()
            self.plot_existing_outlines()
            self.draw()

    def refine(self, num_ref=10):
        if len(self.main_ax.selected_outlines) > 1:
            self.refine_multi(num_ref)
        elif self.sub_ax.nodes:
            self.refine_event(num_ref)

    def refine_event(self, num_ref=10):
        if not self.sub_ax.nodes or self.sub_ax.outline_id is None:
            return

        self.set_status("Refining outline", status="working")
        for i in range(num_ref):
            try:
                self.sub_ax.balloon.evolve(image_percentile=self.image_percentile)
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

        if self.main_ax.selected_outlines:
            self.main_ax.selected_outlines[0].is_modified = True

        self.set_status("Refinement complete ({0} cycles)".format(num_ref))
        self._plot_nodes()
        self.draw()

    def refine_multi(self, num_ref=10):
        if len(self.main_ax.selected_outlines) < 2:
            return

        self.sub_ax.outline_id = None
        for outline_num, outline in enumerate(self.main_ax.selected_outlines):
            self.set_status(
                "Refining {0: 2d} of {1} outlines".format(
                    outline_num + 1,
                    len(self.main_ax.selected_outlines),
                ),
                status="working"
            )

            outline_info = database.getOutlineById(outline.outline_id)
            offset_left, offset_top, _, _ = self.get_offsets((
                outline_info.centre_y,  # pylint: disable=no-member
                outline_info.centre_x,  # pylint: disable=no-member
            ))
            roi = self.load_frame(channel_idx=0)[
                offset_left:offset_left + (self.sub_ax.region_halfwidth * 2),
                offset_top:offset_top + (self.sub_ax.region_halfheight * 2),
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
                len(self.main_ax.selected_outlines)
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
        if self.canvas.count_selected() > 1:
            self.canvas.delete_multi()
        else:
            self.canvas.delete_event()

    def refine(self):
        if self.canvas.count_selected() > 1:
            self.canvas.refine_multi()
        else:
            self.canvas.refine_event()

    def refine_single(self):
        if self.canvas.count_selected() > 1:
            self.canvas.refine_multi(1)
        else:
            self.canvas.refine_event(1)

    def accept(self):
        if self.canvas.count_selected() > 1:
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
        self.experiment_data.image_loader = image_loader

        self.plot = None
        self.tolerance_widget = QtWidgets.QLineEdit()
        self.parent_window = parent_window
        self.window = QtWidgets.QDialog(self.parent_window)

        # must be set later
        self.screen_attributes = dict(
            screen_dpi=None,
            max_width_px=None,
            max_height_px=None,
        )

        self.plot = None
        self.tolerance_widget = QtWidgets.QLineEdit()
        self.parent_window = parent_window
        self.window = QtWidgets.QDialog(self.parent_window)

        # must be set later
        self.screen_attributes = dict(
            screen_dpi=None,
            max_width_px=None,
            max_height_px=None,
        )

    def set_screen_res(self, max_width_px, max_height_px, screen_dpi):
        self.screen_attributes.update({
            "max_width_px": max_width_px,
            "max_height_px": max_height_px,
            "screen_dpi": screen_dpi,
        })

    def _px_to_in(self, pixels):
        return pixels / self.screen_attributes["screen_dpi"]

    def start_outlining(self):
        if not self.screen_attributes["max_width_px"]:
            raise ValueError("Screen resolution has not been set")

        self.window.setModal(True)
        self.window.setGeometry(
            0, 60,
            self.screen_attributes["max_width_px"] * 0.5,
            self.screen_attributes["max_height_px"] - 120
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
            (self._px_to_in(self.screen_attributes["max_width_px"] * 0.5),
             self._px_to_in(self.screen_attributes["max_height_px"] * 0.5),
             self.screen_attributes["screen_dpi"]),
            experiment_data=self.experiment_data,
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
