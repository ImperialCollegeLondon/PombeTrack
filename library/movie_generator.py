#!/usr/bin/env python3

import itertools
import matplotlib.patches
import matplotlib.figure
import numpy as np
import operator
import os
import tifffile

from . import database


class MovieMaker:
    def __init__(
        self,
        experiment_id,
        cell_id,
        descendants,
        outlines,
        frames,
        scale,
        channels,
        image_loader,
        screen_dpi,
    ):
        self.experiment_id = experiment_id
        self.cell_id = cell_id
        self.image_loader = image_loader
        self.screen_dpi = screen_dpi

        cells = sorted(
            self.get_cells(self.cell_id, descendants),
            key=operator.itemgetter(0),
        )
        self.cells = [
            (f, [x[1] for x in d])
            for f, d in itertools.groupby(cells, operator.itemgetter(0))
        ]

        self.create_movie(outlines, frames, scale, channels)

    def get_cells(self, cell_id, descendants):
        cells = [(x.frame_idx, x) for x in database.getOutlinesByCellId(cell_id)]

        if not descendants:
            return cells

        if cells[-1][1].child_id2:
            child1 = database.getOutlineById(cells[-1][1].child_id1)
            child2 = database.getOutlineById(cells[-1][1].child_id2)
            cells.extend(self.get_cells(child1.cell_id, descendants))
            cells.extend(self.get_cells(child2.cell_id, descendants))

        return cells

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

                d = (right_edge - left_edge) - (bottom_edge - top_edge)
                if d > 0:
                    bottom_edge += d / 2
                    top_edge -= d / 2
                elif d < 0:
                    right_edge -= d / 2
                    left_edge += d / 2
                d2 = (right_edge - left_edge) - (bottom_edge - top_edge)

                if left_edge < 0:
                    left_edge = 0
                if right_edge > c_temp.shape[1]:
                    right_edge = c_temp.shape[1]
                if top_edge < 0:
                    top_edge = 0
                if bottom_edge > c_temp.shape[0]:
                    bottom_edge = c_temp.shape[0]

        for frame_idx, cells in self.cells:
            frame_num = frame_idx + 1
            c1 = self.image_loader.load_frame(frame_idx, 0)
            c2, c3 = None, None

            width = bottom_edge - top_edge
            height = right_edge - left_edge
            width_in = width / self.screen_dpi
            if not channels or self.image_loader.num_channels == 1:
                fig = matplotlib.figure.Figure(figsize=(width_in, width_in))
                ax1 = fig.add_subplot(111)
                ax2 = None
                ax3 = None
            elif channels and self.image_loader.num_channels == 2:
                fig = matplotlib.figure.Figure(figsize=(width_in * 2, width_in))
                ax1 = fig.add_subplot(121)
                ax2 = fig.add_subplot(122)
                c2 = self.image_loader.load_frame(frame_idx, 1)
                ax3 = None
            elif channels and self.image_loader.num_channels == 3:
                fig = matplotlib.figure.Figure(figsize=(width_in * 3, width_in))
                ax1 = fig.add_subplot(131)
                ax2 = fig.add_subplot(132)
                ax3 = fig.add_subplot(133)
                c2 = self.image_loader.load_frame(frame_idx, 1)
                c3 = self.image_loader.load_frame(frame_idx, 2)

            for ax in [ax1, ax2, ax3]:
                if not ax:
                    break
                ax.set_aspect("equal")
                ax.axis("off")
                if ax == ax1:
                    ax.imshow(c1, cmap="gray")
                elif ax == ax2:
                    ax.imshow(c2, cmap="binary")
                elif ax == ax3:
                    ax.imshow(c3, cmap="binary")

            if outlines:
                for cell in cells:
                    c = np.load(cell.coords_path) + np.array([cell.offset_left, cell.offset_top])
                    coords = np.array([c[:, 1], c[:, 0]]).T
                    for ax in [ax1, ax2, ax3]:
                        if not ax:
                            break
                        poly = matplotlib.patches.Polygon(
                            coords,
                            edgecolor="r",
                            fill=False,
                            lw=1,
                        )
                        ax.add_patch(poly)

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

            for ax in [ax1, ax2, ax3]:
                if not ax:
                    break
                ax.set_xlim([bottom_edge, top_edge])
                ax.set_ylim([left_edge, right_edge])

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
