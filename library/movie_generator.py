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

        cells = self.get_cells(self.cell_id, settings["descendants"])
        # this bit groups outlines by frame index
        self.cells = [
            (frame_idx, list(cells))
            for frame_idx, cells in
            itertools.groupby(cells, lambda x: x.frame_idx)
        ]

        self.create_movie(
            settings["outlines"],
            settings["frames"],
            settings["scale"],
            settings["channels"],
        )

    def get_cells(self, cell_id, descendants):
        cells = database.getOutlinesByCellId(cell_id)

        if not descendants:
            return cells

        if cells[-1].child_id2:
            child1 = database.getOutlineById(cells[-1].child_id1)
            child2 = database.getOutlineById(cells[-1].child_id2)
            cells.extend(self.get_cells(child1.cell_id, descendants))
            cells.extend(self.get_cells(child2.cell_id, descendants))

        return sorted(
            cells,
            key=lambda x: x.frame_idx
        )

    def create_movie(self, outlines, frames, scale, channels):
        out_path = os.path.join("data", "movies", self.experiment_id, self.cell_id)
        if not os.path.exists(out_path):
            os.makedirs(out_path)

        c_temp = self.image_loader.load_frame(self.cells[0][0], 0)
        left_edge, right_edge = c_temp.shape[1], 0
        top_edge, bottom_edge = c_temp.shape[0], 0
        for frame_idx, cells in self.cells:
            for cell in cells:
                if cell.offset_left < left_edge:
                    left_edge = int(cell.offset_left)

                if cell.offset_left + 150 > right_edge:
                    right_edge = cell.offset_left + 150

                if cell.offset_top < top_edge:
                    top_edge = int(cell.offset_top)

                if cell.offset_top + 150 > bottom_edge:
                    bottom_edge = cell.offset_top + 150

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

        for frame_idx, cells in self.cells:
            channel1 = self.image_loader.load_frame(frame_idx, 0)
            channel2, channel3 = None, None

            width_in = (bottom_edge - top_edge) / self.screen_dpi
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
                ax3 = None
            elif channels and self.image_loader.num_channels == 3:
                fig = matplotlib.figure.Figure(figsize=(width_in * 3, width_in))
                ax1 = fig.add_subplot(131)
                ax2 = fig.add_subplot(132)
                ax3 = fig.add_subplot(133)
                channel2 = self.image_loader.load_frame(frame_idx, 1)
                channel3 = self.image_loader.load_frame(frame_idx, 2)

            for ax_idx, this_ax in enumerate([ax1, ax2, ax3]):
                if not this_ax:
                    break
                this_ax.set_aspect("equal")
                this_ax.axis("off")
                if ax_idx == 0:
                    this_ax.imshow(channel1, cmap="gray")
                elif ax_idx == 1:
                    this_ax.imshow(channel2, cmap="binary")
                elif ax_idx == 2:
                    this_ax.imshow(channel3, cmap="binary")

            if outlines:
                for cell in cells:
                    coords = (np.load(cell.coords_path) +
                              np.array([cell.offset_left, cell.offset_top]))
                    for this_ax in [ax1, ax2, ax3]:
                        if not this_ax:
                            break
                        poly = matplotlib.patches.Polygon(
                            np.array([coords[:, 1], coords[:, 0]]).T,
                            edgecolor="r",
                            fill=False,
                            lw=1,
                        )
                        this_ax.add_patch(poly)

            if frames:
                txt = "F{0}".format(int(frame_idx + 1))
                ax1.text(
                    0.01, 0.99,
                    txt,
                    transform=ax1.transAxes,
                    color="w",
                    fontsize=15,
                    ha="left",
                    va="top",
                )

            if scale:
                px_um = self.image_loader.get_pixel_conversion()
                if not px_um:
                    px_um = 0.16

                pxs = 10 / px_um
                xleft = left_edge + 10
                xright = left_edge + 15
                ytop = top_edge + 10
                ybottom = top_edge + 10 + pxs
                scalebar = matplotlib.patches.Polygon(
                    [(ytop, xleft),
                     (ytop, xright),
                     (ybottom, xright),
                     (ybottom, xleft)],
                    facecolor="w",
                    edgecolor="w",
                    lw=2,
                )
                ax1.add_patch(scalebar)

            for this_ax in [ax1, ax2, ax3]:
                if not this_ax:
                    break
                this_ax.set_xlim([bottom_edge, top_edge])
                this_ax.set_ylim([left_edge, right_edge])

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
