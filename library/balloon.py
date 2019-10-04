#!/usr/bin/env python

""" Provide balloon expansion capability for discovering cell edges.

This module contains the majority of the maths involved in the balloon
algorithm.
Cell outlines are imagined as a series of points (nodes) attached
by strings under strong tension surrounding the cell.

Each node is pushed by an outwards force from the central spine of the cell
(expansion force).
This force is mediated by the local intensity of the image (4x4 pixel region),
where dark areas reduce the expansion force.
The darkness of the local vicinity is determined by comparing the mean
intensity within the region to the whole image, and applying the adjustable
sensitivity parameter.

The nodes are also subjected to a sideways force from each of their neighbours
(neighbour force).
This can be imagined as a form of spring tension force, with the result that
nodes that are further out of the plane of their two neighbours are pulled
backwards more forcefully that those that are close to the plane.
This has the effect of smoothing the outline.

A series of pruning and additive steps is also applied.
Nodes with an overly severe angle to their nearest neighbours are removed, this
somewhat prevents to formation of tangled outlines.
If the distance between two nodes is greater than a threshold value (5 pixels),
a node is inserted halfway between the two nodes on the line connecting them.

This process can be repeated iteratively multiple times to push an outline from
a central position until it matches the dark cell border very accurately.

I named this algorithm "balloon expansion", since I imagine it as a process in
which a balloon is inflated to fit the boundaries of a hard enclosing space,
conforming to any strange shape that could be envisioned.
The elastic nature of a balloon is reflected by the tension-like effect of
neighbours, the inflation of air by the expansion force from the central spine,
and the deforming to fit an arbitrary boundary by the modifying effect of the
image intensity.
"""

import numpy as np
import skimage.draw
import skimage.morphology

class Node:
    """Class representing a single node in an outline.

    Defines methods for calculating the forces acting on a node to determine
    which direction the node will move in the next step of the algorithm.

    The x and y coordinates recorded in this class are used to define the
    outline coordinates used throughout PombeTrack.

    Arguments:
        coord (tuple): x, y coordinate of the node position.

    Attributes:
        x (float):         x coordinate of the node position.
        y (float):         y coordinate of the node position.
        buff_x (float):    x coordinate of the node after applying forces.
        buff_y (float):    y coordinate of the node after applying forces.
        neighbour1 (Node): node to the left of this node.
        neighbour2 (Node): node to the right of this node.

    *left and right in this context are arbritrary
    """
    def __init__(self, coord):
        self.x = coord[0]
        self.y = coord[1]
        self.buff_x = None
        self.buff_y = None
        self.neighbour1 = None
        self.neighbour2 = None

    def apply_changes(self):
        """Update the x and y coordinates of the node."""
        self.x = self.buff_x
        self.y = self.buff_y
        self.buff_x = None
        self.buff_y = None

    def get_distance(self, node):
        """Calculate the distance between this node and another.

        Arguments:
            node (Node instance): node to compare to.

        Returns:
            distance (float): distance between the nodes (in pixels).
        """
        return np.sqrt(
            (self.x - node.x) ** 2 +
            (self.y - node.y) ** 2
        )

    def connect(self, node1, node2):
        """Sets the neighbours of the node.

        Arguments:
            node1 (Node instance): "left-hand" neighbour.
            node2 (node instance): "right-hand" neighbour.

        Returns:
            None
        """
        self.neighbour1 = node1
        self.neighbour2 = node2

    @staticmethod
    def _get_angle(point1, point2):
        """Calculate the angle between two points relative to the x axis.

        Arguments:
            point1 (tuple): x, y coordinates.
            point2 (tuple): x, y coordinates.

        Returns:
            angle (float): angle between the points.

        Uses simple trigonometry.
        """
        ydelta = point2[0] - point1[0]
        xdelta = point2[1] - point1[1]
        if xdelta == 0:
            hypot = np.sqrt(xdelta ** 2 + ydelta ** 2)
            theta = np.arcsin(ydelta / hypot)
        elif ydelta == 0:
            hypot = np.sqrt(xdelta ** 2 + ydelta ** 2)
            theta = np.arccos(xdelta / hypot)
        else:
            theta = np.arctan(ydelta / xdelta)
        return theta

    def neighbour_force(self, point):
        """Calculate the force applied by the neighbours to a hypothetical node
        position.

        Arguments:
            point (tuple): x, y coordinates of the hypothetical position.

        Returns:
            force (float): combination of the distance and position of each
                           neighbour.
        """
        left_vector = np.array([
            self.neighbour1.x - point[0],
            self.neighbour1.y - point[1],
        ])
        right_vector = np.array([
            self.neighbour2.x - point[0],
            self.neighbour2.y - point[1],
        ])
        return 0.1 * (left_vector + right_vector)

    def get_skeleton_centre(self, skel_coords):
        """Calculate the nearest point along the skeletonized image to the node.

        Arguments:
            skel_coords (Nx2 ndarray): x, y coordinates of the skeleton.

        Returns:
            skel_centre (numpy vector): x, y coordinates of the nearest
                                        skeleton point.

        Used to determine whether the expansion force originates from.
        """
        # get nearest pixel coordinate from node to skel
        distances = np.sqrt(
            ((skel_coords - [self.x, self.y]) ** 2).sum(axis=1)
        )
        return skel_coords[np.argmin(distances)]

    def expansion_force(self, expansion_point):
        """Calculate the expansion force from the cell skeleton.

        Arguments:
            expansion_point (vector): x, y coordinates of the expansion point.

        Returns:
            expansion_force (vector): expansion force acting on the node split
                                      into x and y directions.

        Determined by the angle from the point to the expansion point.
        """
        angle = self._get_angle(expansion_point, (self.x, self.y))
        force_x = np.cos(angle)
        force_y = np.sin(angle)
        scale_factor = 3
        if self.y < expansion_point[1]:
            return scale_factor * np.array([-force_y, -force_x])
        return scale_factor * np.array([force_y, force_x])

    def apply_force(self, skel_coords, image, percentile=5):
        """Determine the movement of a node after calculating all forces.

        Arguments:
            skel_coords (Nx2 ndarray): coordinates of the cell skeleton.
            image (NxM ndarray):       matrix of image intensity values.
            percentile (int or float): the level which is considered to be dark
                                       in the image; defaults to 5.

        Returns:
            None

        Determines the expansion point for the node, applies the expansion
        force, neighbour force, and calculates the image force.

        The image force is calculated from a 4x4 region around the node (in its
        current position) and the percentile setting for what is considered to
        be "dark" compared to the whole image.

        The calculated node position is updated in a buffer.
        This means that the order of expansion applied to each node in the
        outline doesn't matter.
        """
        if (self.x < 1 or
                self.x > image.shape[0] - 1 or
                self.y < 1 or
                self.y > image.shape[1] - 1):
            self.set_position((self.x, self.y))
            return
        expansion_point = self.get_skeleton_centre(skel_coords)
        cforce = self.expansion_force(expansion_point)
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
        """Set the buffered new position of the node after evolution is finished."""
        self.buff_x = updated[0]
        self.buff_y = updated[1]

    def __repr__(self):
        return "({0}, {1}) -- {2}, {3} -- ({4}, {5})".format(
            self.neighbour1.x, self.neighbour1.y,
            self.x, self.y,
            self.neighbour2.x, self.neighbour2.y,
        )


class Balloon:
    """ Class for handling all nodes in a balloon model of a cell outline.

    As mentioned in the module docstring, an outline is imagined as a series of
    points attached together that encircle the cell contents.
    Each point is handled by the `Node` class above.

    The Balloon class handles operations that affect all nodes, such as adding
    and removing nodes, node movements, the points from which the expansion
    force originates, and overall balloon evolution.

    Arguments:
        initial_coords (Nx2 array): initial x, y coordinates of balloon nodes.
        base_image (NxM array):     brightfield image containing the cell.

    Attributes:
        nodes (list):           all nodes defining the balloon in sequential order.
        base_image (NxM array): brightfield image.
        centre (tuple):         x, y coordinates of the centre of the balloon.
        refining_cycles (int):  number of iterations of the evolution.
    """
    def __init__(self, initial_coords, base_image):
        self.nodes = self.create_nodes(initial_coords)
        self.base_image = base_image
        self.centre = self.get_centre()
        self.refining_cycles = 0

    @staticmethod
    def create_nodes(coords):
        """Convert a series of coordinates into Node objects.

        Arguments:
            coords (Nx2 array): x y coordinates of each point.

        Returns:
            nodes (list): list of Node objects for each coordinate.

        This function also assumes that the coordinates are arranged in order,
        thus connects each Node to its preceding and succeeding Node (with
        wrapping).
        """
        nodes = []
        for coord in coords:
            nodes.append(Node(coord))

        for i, node1 in enumerate(nodes):
            node0 = nodes[i - 1]
            if i == len(nodes) - 1:
                node2 = nodes[0]
            else:
                node2 = nodes[i + 1]
            node1.connect(node0, node2)
        return nodes

    def get_centre(self):
        """Get the centre (centroid) of the area defined by the balloon.

        The centroid is defined simply as the mean x and mean y position of all
        nodes.

        Arguments:
            N/A

        Returns:
            centre_x (float): x coordinate of the centroid positions.
            centre_y (float): y coordinate of the centroid positions.
        """
        # just get the centroid
        # perhaps try something like:
        # https://github.com/mapbox/polylabel/blob/master/polylabel.js
        # in the future
        coords = np.array([(n.x, n.y) for n in self.nodes])
        centre_x = coords[:, 0].mean()
        centre_y = coords[:, 1].mean()
        return centre_x, centre_y

    def remove_node(self, node):
        """Remove a node from the balloon outline.

        Deletes a node and reconnects its neighbours to themselves.

        Arguments:
            node (Node): the node to be removed.

        Returns:
            None
        """
        self.nodes.pop(self.nodes.index(node))
        node1 = node.neighbour1
        node2 = node.neighbour2
        node1.neighbour2 = node2
        node2.neighbour1 = node1

    def prune_branches(self):
        """Remove errant nodes from the balloon outline.

        Checks whether the angle between the lines connecting the node to its
        two neighbours is less than a threshold value.
        This threshold is hard-coded to 2π/3 radians (120°).
        If the threshold is exceeded, the node is removed, and its neighbours
        are connected to one another.

        Arguments:
            N/A

        Returns:
            new_nodes (list):
        """
        if len(self.nodes) < 20:
            return self.nodes

        neighbour_min_angle = np.pi*2 / 3
        for node1 in self.nodes:
            node0 = node1.neighbour1
            node2 = node1.neighbour2
            n0vector = [node1.x - node0.x, node1.y - node0.y]
            n2vector = [node2.x - node1.x, node2.y - node1.y]
            angle = np.arctan2(
                np.linalg.norm(
                    np.cross(n0vector, n2vector)
                ),
                np.dot(n0vector, n2vector)
            )
            if angle > neighbour_min_angle:
                # prune it out!
                node0.neighbour2 = node2
                node2.neighbour1 = node0
                node1.pruned = True

        new_nodes = []
        for node in self.nodes:
            if not hasattr(node, "pruned"):
                new_nodes.append(node)

        return new_nodes

    def insert_nodes(self):
        """Add nodes to the outline.

        Checks the distance between each pair of nodes and adds a new node if
        the distance is larger than a threshold value, currently hard-coded to
        5 pixels.

        When adding a node, it is placed at the halfway point between the two
        nodes, and the existing nodes are connected appropriately as
        neighbours.

        Arguments:
            N/A

        Returns:
            new_nodes (list): record of all the new nodes.

        """
        neighbour_max_distance = 5
        new_nodes = []
        for node in self.nodes:
            left_distance = node.get_distance(node.neighbour1)
            right_distance = node.get_distance(node.neighbour2)
            if left_distance > neighbour_max_distance:
                # halfway
                half_point = (
                    node.x + (node.neighbour1.x - node.x) / 2,
                    node.y + (node.neighbour1.y - node.y) / 2
                )
                new_node = Node(half_point)
                node.neighbour1.connect(node.neighbour1.neighbour1, new_node)
                new_node.connect(node.neighbour1, node)
                node.connect(new_node, node.neighbour2)
                new_nodes.append(new_node)
            new_nodes.append(node)

            if right_distance > neighbour_max_distance:
                # halfway
                half_point = (
                    node.x + (node.neighbour2.x - node.x) / 2,
                    node.y + (node.neighbour2.y - node.y) / 2
                )
                new_node = Node(half_point)
                node.neighbour2.connect(new_node, node.neighbour2.neighbour2)
                new_node.connect(node, node.neighbour2)
                node.connect(node.neighbour1, new_node)
                new_nodes.append(new_node)

        return new_nodes

    def evolve(self, image_percentile=5):
        """Perform single iteration of balloon algorithm.

        Determines the central spine of the outline (skeleton) which become the
        expansion points.
        Iterates through each node and calls Node.apply_force.
        Applies each change after all the nodes have been adjusted.
        Rejects changes if the cell area becomes lower than 50 pixels^2.
        Errant nodes are pruned with Balloon.prune_branches.
        Finally, nodes are added with Balloon.insert_nodes.

        Arguments:
            image_percentile (int or float):
                The level which is considered to be dark in the image, see
                Node.apply_force; defaults to 5.

        Returns:
            delta_area (float):
                Change in area (pixels^2) for the original outline compared to
                the new outline.
        """
        self.refining_cycles += 1
        node_positions = np.array([(x.x, x.y) for x in self.nodes])
        poly_rr, poly_cc = skimage.draw.polygon(
            node_positions[:, 0],
            node_positions[:, 1],
            shape=self.base_image.shape
        )
        poly_img = np.zeros(self.base_image.shape)
        poly_img[poly_rr, poly_cc] = 1
        skel = skimage.morphology.skeletonize(poly_img)
        skel_coords = np.array(np.where(skel)).T

        for node in self.nodes:
            node.apply_force(skel_coords, self.base_image, image_percentile)

        original_positions = np.zeros((len(self.nodes), 2))
        new_positions = np.zeros((len(self.nodes), 2))
        for i, node in enumerate(self.nodes):
            original_positions[i, :] = (node.x, node.y)
            new_positions[i, :] = (node.buff_x, node.buff_y)
            node.apply_changes()

        new_area = self.get_area(new_positions)
        if new_area < 50:
            for i, node in enumerate(self.nodes):
                node.x = original_positions[i, 0]
                node.y = original_positions[i, 1]
            raise ValueError("Area too small")

        self.nodes = self.prune_branches()
        self.nodes = self.insert_nodes()
        return np.abs(new_area - self.get_area(original_positions))

    def get_coordinates(self):
        """Obtain the coordinates of all nodes in the outline.

        Arguments:
            N/A

        Returns:
            coords (np.array; Nx2): xy coordinates of each node.
        """
        return np.array([(n.x, n.y) for n in self.nodes])

    def get_area(self, coords=None):
        """Calculate the area within node coordinates.

        Uses a Shoelace algorithm.

        Arguments:
            coords (np.array; Nx2):
                xy coordinate positions of each node, defaults to the current
                outline coordinates if argument is omitted.

        Returns:
            area (float):
                Area within the coordinate bounds (pixels^2).
        """
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
    """Generate a set of coordinates around a single point.

    Creates a circular set of coordinate points around a given central
    coordinate.

    Arguments:
        centre (tuple):
            x, y coordinates of the central point.
        radius (float):
            radius in pixels of the desired circle.
        num_nodes (int):
            how many node positions to create in the circle.

    Returns:
        nodes (np.array; Nx2):
            x, y coordinates of each node generated.
    """
    nodes = np.array([
        (centre[0] + radius * np.sin(x),
         centre[1] + radius * np.cos(x))
        for x in np.linspace(0, 2 * np.pi, num_nodes)[:-1]
    ])
    return nodes
