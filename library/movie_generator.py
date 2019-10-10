#!/usr/bin/env python3

import itertools
import os

import matplotlib.patches
import matplotlib.figure
import numpy as np

from . import database


class MovieMaker:
    def __init__(self, **settings):
        self.experiment_id = settings["experiment_id"]
        self.cell_id = settings["cell_id"]
        self.image_loader = settings["image_loader"]
        self.screen_dpi = settings["screen_dpi"]

        outlines = self.get_outlines(self.cell_id, settings["descendants"])
        # this bit groups outlines by frame index
        self.outlines = [
            (frame_idx, list(outlines))
            for frame_idx, outlines in
            itertools.groupby(outlines, lambda x: x.frame_idx)
        ]

        self.create_movie(
            settings["draw_outlines"],
            settings["frames"],
            settings["scale"],
            settings["channels"],
        )

    def get_outlines(self, cell_id, descendants):
        outlines = database.getOutlinesByCellId(cell_id)

        if not descendants:
            return outlines

        if outlines[-1].child_id2:
            child1 = database.getOutlineById(outlines[-1].child_id1)
            child2 = database.getOutlineById(outlines[-1].child_id2)
            # pylint warning disabled, OutlineRow does have a cell_id attribute
            outlines.extend(self.get_outlines(child1.cell_id, descendants))  # pylint: disable=no-member
            outlines.extend(self.get_outlines(child2.cell_id, descendants))  # pylint: disable=no-member

        return sorted(
            outlines,
            key=lambda x: x.frame_idx
        )

    def get_edges(self):
        c_temp = self.image_loader.load_frame(self.outlines[0][0], 0)
        left_edge, right_edge = c_temp.shape[1], 0
        top_edge, bottom_edge = c_temp.shape[0], 0
        for frame_data in self.outlines:
            for outline in frame_data[1]:
                if outline.offset_left < left_edge:
                    left_edge = int(outline.offset_left)

                if outline.offset_left + 150 > right_edge:
                    right_edge = outline.offset_left + 150

                if outline.offset_top < top_edge:
                    top_edge = int(outline.offset_top)

                if outline.offset_top + 150 > bottom_edge:
                    bottom_edge = outline.offset_top + 150

                delta = (right_edge - left_edge) - (bottom_edge - top_edge)
                if delta > 0:
                    bottom_edge += delta / 2
                    top_edge -= delta / 2
                elif delta < 0:
                    right_edge -= delta / 2
                    left_edge += delta / 2

                if left_edge < 0:
                    left_edge = 0
                if right_edge > c_temp.shape[1]:
                    right_edge = c_temp.shape[1]
                if top_edge < 0:
                    top_edge = 0
                if bottom_edge > c_temp.shape[0]:
                    bottom_edge = c_temp.shape[0]

        return left_edge, right_edge, top_edge, bottom_edge

    def create_axes(self, frame_idx, edges, channels):
        channel1 = self.image_loader.load_frame(frame_idx, 0)
        channel2, channel3 = None, None

        width_in = (edges[3] - edges[2]) / self.screen_dpi
        if not channels or self.image_loader.num_channels == 1:
            fig = matplotlib.figure.Figure(figsize=(width_in, width_in))
            ax1 = fig.add_subplot(111)
            ax2 = None
            ax3 = None
        elif channels and self.image_loader.num_channels == 2:
            fig = matplotlib.figure.Figure(figsize=(width_in * 2, width_in))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            channel2 = self.image_loader.load_frame(frame_idx, 1)
            ax2.imshow(channel2, cmap="binary")
            ax3 = None
        elif channels and self.image_loader.num_channels == 3:
            fig = matplotlib.figure.Figure(figsize=(width_in * 3, width_in))
            ax1 = fig.add_subplot(131)
            ax2 = fig.add_subplot(132)
            channel2 = self.image_loader.load_frame(frame_idx, 1)
            ax2.imshow(channel2, cmap="binary")
            ax3 = fig.add_subplot(133)
            channel3 = self.image_loader.load_frame(frame_idx, 2)
            ax3.imshow(channel3, cmap="binary")

        ax1.imshow(channel1, cmap="gray")

        axes = [axis for axis in [ax1, ax2, ax3] if axis]
        for axis in axes:
            axis.set_aspect("equal")
            axis.axis("off")

        return fig, axes

    @staticmethod
    def add_outlines(outlines, axes):
        for outline in outlines:
            coords = (np.load(outline.coords_path) +
                      np.array([outline.offset_left, outline.offset_top]))
            for axis in axes:
                poly = matplotlib.patches.Polygon(
                    np.array([coords[:, 1], coords[:, 0]]).T,
                    edgecolor="r",
                    fill=False,
                    lw=1,
                )
                axis.add_patch(poly)

    @staticmethod
    def add_frame_text(frame_idx, axes):
        txt = "F{0}".format(int(frame_idx + 1))
        axes[0].text(
            0.01, 0.99,
            txt,
            transform=axes[0].transAxes,
            color="w",
            fontsize=15,
            ha="left",
            va="top",
        )

    def add_scale_bar(self, axes, edges):
        px_um = self.image_loader.get_pixel_conversion()
        if not px_um:
            px_um = 0.16

        pxs = 10 / px_um
        xleft = edges[0] + 10
        xright = edges[0] + 15
        ytop = edges[2] + 10
        ybottom = edges[2] + 10 + pxs
        scalebar = matplotlib.patches.Polygon(
            [(ytop, xleft),
             (ytop, xright),
             (ybottom, xright),
             (ybottom, xleft)],
            facecolor="w",
            edgecolor="w",
            lw=2,
        )
        axes[0].add_patch(scalebar)

    @staticmethod
    def set_limits(axes, edges):
        for axis in axes:
            axis.set_xlim(edges[2:])
            axis.set_ylim(edges[:2])

    def create_movie(self, draw_outlines, frames, scale, channels):
        out_path = os.path.join("data", "movies",
                                self.experiment_id, self.cell_id)
        if not os.path.exists(out_path):
            os.makedirs(out_path)

        for frame_idx, outlines in self.outlines:
            edges = self.get_edges()
            fig, axes = self.create_axes(frame_idx, edges, channels)
            if draw_outlines:
                self.add_outlines(outlines, axes)

            if frames:
                self.add_frame_text(frame_idx, axes)

            if scale:
                self.add_scale_bar(axes, edges)

            self.set_limits(axes, edges)

            fig.tight_layout()
            fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
            out_fn = os.path.join(out_path, "f{0:03d}.png".format(frame_idx + 1))
            fig.savefig(
                out_fn,
                dpi=self.screen_dpi,
            )
            print("Written {0}".format(out_fn))
            del fig

        print("Done writing movie")
