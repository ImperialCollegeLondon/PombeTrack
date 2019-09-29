#!/usr/bin/env python

import numpy as np
import os
import skimage.draw
import skimage.morphology
# import scipy.ndimage

class Node(object):
    def __init__(self, coord):
        self.x = coord[0]
        self.y = coord[1]
        self._connected = False
        self.buff_x = None
        self.buff_y = None

    def apply_changes(self):
        self.x = self.buff_x
        self.y = self.buff_y
        self.buff_x = None
        self.buff_y = None

    def _get_distance(self, n):
        return np.sqrt(
            (self.x - n.x) ** 2 +
            (self.y - n.y) ** 2
        )

    def connect(self, n1, n2):
        self.neighbour1 = n1
        self.neighbour2 = n2

    def _get_angle(self, p1, p2):
        ydelta = p2[0] - p1[0]
        xdelta = p2[1] - p1[1]
        if xdelta == 0:
            hypot = np.sqrt(xdelta ** 2 + ydelta ** 2)
            theta = np.arcsin(ydelta / hypot)
        elif ydelta == 0:
            hypot = np.sqrt(xdelta ** 2 + ydelta ** 2)
            theta = np.arccos(xdelta / hypot)
        else:
            theta = np.arctan(ydelta / xdelta)
        return theta

    def neighbour_force(self, new_pos):
        left_vector = np.array([
            self.neighbour1.x - new_pos[0],
            self.neighbour1.y - new_pos[1],
        ])
        right_vector = np.array([
            self.neighbour2.x - new_pos[0],
            self.neighbour2.y - new_pos[1],
        ])
        return 0.1 * (left_vector + right_vector)

    def get_skeleton_centre(self, skel_coords):
        # get nearest pixel coordinate from node to skel
        distances = np.sqrt(
            ((skel_coords - [self.x, self.y]) ** 2).sum(axis=1)
        )
        return skel_coords[np.argmin(distances)]

    def centroid_force(self, expansion_point):
        angle = self._get_angle(expansion_point, (self.x, self.y))
        force_x = np.cos(angle)
        force_y = np.sin(angle)
        scale_factor = 3
        if self.y < expansion_point[1]:
            return scale_factor * np.array([-force_y, -force_x])
        return scale_factor * np.array([force_y, force_x])

    def contract(self, centre):
        nforce = self.neighbour_force((self.x, self.y))
        self.set_position((
            self.x + nforce[0],
            self.y + nforce[1]
        ))

    def apply_force(self, skel_coords, image, percentile=5):
        if (self.x < 1 or
                self.x > image.shape[0] - 1 or
                self.y < 1 or
                self.y > image.shape[1] - 1):
            self.set_position((self.x, self.y))
            return
        expansion_point = self.get_skeleton_centre(skel_coords)
        cforce = self.centroid_force(expansion_point)
        node_updated = (self.x + cforce[0], self.y + cforce[1])
        nforce = self.neighbour_force(node_updated)
        vicinity = image[
            int(round(self.x)) - 1:int(round(self.x)) + 2,
            int(round(self.y)) - 1:int(round(self.y)) + 2,
        ].flatten()
        image_multiplier = np.log(vicinity.mean() / np.percentile(image, percentile))
        self.set_position((
            self.x + (image_multiplier * cforce[0]) + nforce[0],
            self.y + (image_multiplier * cforce[1]) + nforce[1],
        ))

    def set_position(self, updated):
        self.buff_x = updated[0]
        self.buff_y = updated[1]

    def __repr__(self):
        return "({0}, {1}) -- {2}, {3} -- ({4}, {5})".format(
            self.neighbour1.x, self.neighbour1.y,
            self.x, self.y,
            self.neighbour2.x, self.neighbour2.y,
        )


class Balloon(object):
    def __init__(self, initial_coords, base_image):
        self.mode = "initialise"
        self.nodes = self.create_nodes(initial_coords)
        self.base_image = base_image
        self.centre = self.get_centre()
        self.refining_cycles = 0

    def create_nodes(self, coords):
        nodes = []
        for coord in coords:
            nodes.append(Node(coord))

        for i in range(len(nodes)):
            n0 = nodes[i - 1]
            n1 = nodes[i]
            if i == len(nodes) - 1:
                n2 = nodes[0]
            else:
                n2 = nodes[i + 1]
            n1.connect(n0, n2)
        return nodes

    def get_centre(self):
        # just get the centroid
        # perhaps try something like:
        # https://github.com/mapbox/polylabel/blob/master/polylabel.js
        # in the future
        coords = np.array([(n.x, n.y) for n in self.nodes])
        Cx = coords[:, 0].mean()
        Cy = coords[:, 1].mean()
        return Cx, Cy

    def remove_node(self, n):
        self.nodes.pop(self.nodes.index(n))
        n1 = n.neighbour1
        n2 = n.neighbour2
        n1.neighbour2 = n2
        n2.neighbour1 = n1

    def prune_branches(self):
        if len(self.nodes) < 20:
            return self.nodes

        neighbour_min_angle = np.pi*2 / 3
        for n1 in self.nodes:
            n0 = n1.neighbour1
            n2 = n1.neighbour2
            n0vector = [n1.x - n0.x, n1.y - n0.y]
            n2vector = [n2.x - n1.x, n2.y - n1.y]
            angle = np.arctan2(
                np.linalg.norm(
                    np.cross(n0vector, n2vector)
                ),
                np.dot(n0vector, n2vector)
            )
            if angle > neighbour_min_angle:
                # prune it out!
                n0.neighbour2 = n2
                n2.neighbour1 = n0
                n1.pruned = True

        new_nodes = []
        for n in self.nodes:
            if not hasattr(n, "pruned"):
                new_nodes.append(n)

        return new_nodes

    def insert_nodes(self):
        neighbour_max_distance = 5
        new_nodes = []
        for n in self.nodes:
            left_distance = n._get_distance(n.neighbour1)
            right_distance = n._get_distance(n.neighbour2)
            if left_distance > neighbour_max_distance:
                # halfway
                half_point = (
                    n.x + (n.neighbour1.x - n.x) / 2,
                    n.y + (n.neighbour1.y - n.y) / 2
                )
                new_node = Node(half_point)
                n.neighbour1.connect(n.neighbour1.neighbour1, new_node)
                new_node.connect(n.neighbour1, n)
                n.connect(new_node, n.neighbour2)
                new_nodes.append(new_node)
            new_nodes.append(n)

            if right_distance > neighbour_max_distance:
                # halfway
                half_point = (
                    n.x + (n.neighbour2.x - n.x) / 2,
                    n.y + (n.neighbour2.y - n.y) / 2
                )
                new_node = Node(half_point)
                n.neighbour2.connect(new_node, n.neighbour2.neighbour2)
                new_node.connect(n, n.neighbour2)
                n.connect(n.neighbour1, new_node)
                new_nodes.append(new_node)

        return new_nodes

    def evolve(self, image_percentile=5):
        self.refining_cycles += 1
        node_positions = np.array([(x.x, x.y) for x in self.nodes])
        poly_rr, poly_cc = skimage.draw.polygon(
            node_positions[:, 0],
            node_positions[:, 1],
            shape=self.base_image.shape
        )
        poly_img = np.zeros_like(self.base_image)
        poly_img[poly_rr, poly_cc] = 1
        skel = skimage.morphology.skeletonize(poly_img)
        skel_coords = np.array(np.where(skel == True)).T

        for n in self.nodes:
            n.apply_force(skel_coords, self.base_image, image_percentile)

        original_positions = np.zeros((len(self.nodes), 2))
        new_positions = np.zeros((len(self.nodes), 2))
        for i, n in enumerate(self.nodes):
            original_positions[i, :] = (n.x, n.y)
            new_positions[i, :] = (n.buff_x, n.buff_y)
            n.apply_changes()

        original_area = self.get_area(original_positions)
        new_area = self.get_area(new_positions)
        delta_area = np.abs(new_area - original_area)
        if new_area < 50:
            for i, n in enumerate(self.nodes):
                n.x = original_positions[i, 0]
                n.y = original_positions[i, 1]
            raise ValueError("Area too small")

        self.nodes = self.prune_branches()
        self.nodes = self.insert_nodes()
        return delta_area

    def get_coordinates(self):
        return np.array([(n.x, n.y) for n in self.nodes])

    def get_area(self, coords=None):
        if coords is None:
            coords = np.array([(n.x, n.y) for n in self.nodes])

        area = np.dot(
            coords[:, 0],
            np.roll(coords[:, 1], 1)
        ) - np.dot(
            coords[:, 1],
            np.roll(coords[:, 0], 1)
        )
        return area * 0.5


def initial_nodes(centre, radius, num_nodes):
    initial_nodes = np.array([
        (centre[0] + radius * np.sin(x),
         centre[1] + radius * np.cos(x))
        for x in np.linspace(0, 2 * np.pi, num_nodes)[:-1]
    ])
    return initial_nodes


if __name__ == "__main__":
    base_image = np.load("test.npy")
    centre = (45, 43)
    radius = 5
    num_nodes = 20
    init = initial_nodes(centre, radius, num_nodes)
    B = Balloon(init, base_image)
    i = 0
    while True:
        if i % 50 == 0:
            B.evolve(True)
        else:
            B.evolve()
        i += 1
