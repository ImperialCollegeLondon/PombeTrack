#!/usr/bin/env python3

import matplotlib.patches
import matplotlib.pyplot as plt
import tifffile

from . import database


class MovieMaker:
    def __init__(
        self, experiment_id, cell_id,
        descendants,
        outlines,
        frames,
        scale,
        channels,
        imagej,
    ):
        self.experiment_id = experiment_id
        self.cell_id = cell_id

        self.cells = self.get_cells(self.cell_id, descendants)
        if imagej:
            self.create_individual_files(outlines, frames, scale, channels)
        else:
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

    def create_individual_files(self, outlines, frames, scale, channels):
        pass

    def create_movie(self, outlines, frames, scale, channels):
        pass
