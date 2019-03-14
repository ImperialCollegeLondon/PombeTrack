#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
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

        distances = []
        for x, y in skel_coords:
            distance = np.sqrt((x - self.x) ** 2 + (y - self.y) ** 2)
            distances.append(distance)
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
        if (self.x < 11 or
                self.x > image.shape[0] - 11 or
                self.y < 11 or
                self.y > image.shape[1] - 11):
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
        image_multiplier = vicinity.mean() / np.percentile(image, percentile) - 1
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

    def refine(self, display=False):
        # contract for 50 iterations
        for i in range(50):
            for n in self.nodes:
                n.contract(self.centre)

            for n in self.nodes:
                n.apply_changes()

        # evolve again
        self.run_evolution(display_step=1000)

    def remove_node(self, n):
        self.nodes.pop(self.nodes.index(n))
        n1 = n.neighbour1
        n2 = n.neighbour2
        n1.neighbour2 = n2
        n2.neighbour1 = n1

    def prune_branches(self):
        if len(self.nodes) < 20:
            return self.nodes

        neighbour_min_angle = np.pi / 3
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

    def evolve(self, display=False, image_percentile=5):
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
        # original_area = 0.5 * np.abs(
        #     np.dot(
        #         original_positions[:, 0],
        #         np.roll(original_positions[:, 1], 1)
        #     ) - np.dot(
        #         original_positions[:, 1],
        #         np.roll(original_positions[:, 0], 1)
        #     )
        # )
        new_area = self.get_area(new_positions)
        # new_area = 0.5 * np.abs(
        #     np.dot(
        #         new_positions[:, 0],
        #         np.roll(new_positions[:, 1], 1)
        #     ) - np.dot(
        #         new_positions[:, 1],
        #         np.roll(new_positions[:, 0], 1)
        #     )
        # )
        delta_area = np.abs(new_area - original_area)

        if display:
            f = plt.figure()
            ax = f.add_subplot(111)
            ax.imshow(self.base_image, cmap="gray")
            ax.plot(skel_coords[:, 1], skel_coords[:, 0], color="b", marker="o", ms=3, ls="none")
            ax.set_aspect("equal")
            ax.plot(original_positions[:, 1], original_positions[:, 0], marker="o", ms=5, color="k")
            ax.plot(new_positions[:, 1], new_positions[:, 0], marker="o", ms=5, color="r")
            ax.set_title("Change in area: {0:.5f}".format(delta_area))

            if type(display) is int:
                ax.set_title("Iteration {0}: Change in area: {1:.5f}".format(display, delta_area))
                if not os.path.exists("out"):
                    os.mkdir("out")
                f.savefig("out/{0}.png".format(display), dpi=600)
                plt.close()
            else:
                plt.show()
                plt.close()

        self.nodes = self.prune_branches()
        self.nodes = self.insert_nodes()
        return delta_area

    def run_evolution(
        self,
        display_step=None,
        min_iterations=100,
        max_iterations=2000,
        threshold=0.04,
        save_progress=False,
        display_final=False,
    ):
        i = 1
        threshold_hits = 0
        while True:
            if display_step and i % display_step == 0:
                if save_progress:
                    delta_area = self.evolve(i)
                else:
                    delta_area = self.evolve(display_final)
            else:
                delta_area = self.evolve()

            if delta_area < threshold and i > min_iterations:
                threshold_hits += 1
                if threshold_hits > 5:
                    print("Ending after {0} iterations".format(i))
                    self.evolve(display_final)
                    self.mode = "refinement"
                    return
            elif i == max_iterations:
                print("Ending after {0} iterations".format(i))
                self.evolve(display_final)
                self.mode = "refinement"
                return
            i += 1

    def _keypress(self, evt):
        if self.mode == "refinement":
            if evt.key == "enter":
                self.mode = "confirmed"
                plt.close()
            elif evt.key == "escape":
                self.mode = "adjust"
                self.adjust_coordinates()

        elif self.mode == "adjust":
            if evt.key == "escape":
                self.mode = "rejected"
            elif evt.key == "enter":
                self.mode = "confirmed"
            elif evt.key == ".":
                self.evolve()
                self._replot_nodes()
                plt.draw()
            elif evt.key == "r":
                for i in range(10):
                    self.evolve()
                self._replot_nodes()
                plt.draw()

            if evt.key == "escape" or evt.key == "enter":
                plt.close()

        else:
            print("Keypress:", evt.key, "mode ({0})".format(self.mode))

    def _click(self, evt):
        if self.mode != "adjust":
            return

        if evt.button == 1:
            for p in self.figure_patches:
                if p.contains_point((evt.x, evt.y)):
                    p.set_facecolor("r")
                    self.dragging = p
                    break
        elif evt.button == 3:
            for p in self.figure_patches:
                if p.contains_point((evt.x, evt.y)):
                    self.remove_node(p.this_node)
                    self._replot_nodes()
                    plt.draw()
                    break

    def _release(self, evt):
        if self.mode != "adjust":
            return

        if self.dragging:
            self.dragging.set_facecolor("y")
            # update node
            self.dragging.this_node.set_position((evt.ydata, evt.xdata))
            self.dragging.this_node.apply_changes()

            # update line
            ax = plt.gca()
            ax.lines.pop(0)
            ax.plot(
                [n.y for n in self.nodes],
                [n.x for n in self.nodes],
                color="y",
                lw=2,
            )

            self.dragging = False
            plt.draw()

    def _motion(self, evt):
        if self.mode != "adjust":
            return

        if not self.dragging:
            return

        self.dragging.center = evt.xdata, evt.ydata
        plt.draw()

    def get_coordinates(self, accept=False):
        if accept:
            return np.array([(n.x, n.y) for n in self.nodes])

        # adjust them first
        f = plt.figure(figsize=(15, 15))
        f.canvas.mpl_connect("key_press_event", self._keypress)
        f.canvas.mpl_connect("button_press_event", self._click)
        f.canvas.mpl_connect("button_release_event", self._release)
        f.canvas.mpl_connect("motion_notify_event", self._motion)
        ax = f.add_subplot(111)
        ax.imshow(self.base_image, cmap="gray")
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title("Confirm segmentation: ENTER/ESCAPE")
        self._replot_nodes(patches=False)
        plt.show()

        if self.mode == "confirmed":
            return np.array([(n.x, n.y) for n in self.nodes])
        else:
            return None

    def _replot_nodes(self, line=True, patches=True):
        # replot patches
        ax = plt.gca()
        if line:
            try:
                ax.lines.pop(0)
            except IndexError:
                pass
            ax.plot(
                [n.y for n in self.nodes],
                [n.x for n in self.nodes],
                color="y",
                lw=5,
            )

        if patches:
            # replot patches
            for i in range(len(self.figure_patches)):
                p = self.figure_patches.pop()
                p.remove()
               
            for n in self.nodes:
                patch = matplotlib.patches.Circle(
                    (n.y, n.x),
                    1,
                    fc="y",
                )
                patch.this_node = n
                self.figure_patches.append(patch)
                ax.add_artist(patch)

    def adjust_coordinates(self):
        f = plt.gcf()
        ax = plt.gca()
        # ax.lines[0].set_marker("o")
        # ax.lines[0].set_markersize(10)
        self.figure_patches = []
        self.dragging = False
        self._replot_nodes()
        ax.set_title("Edit segmentation: ENTER (accept changes)")
        plt.draw()

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
        return area



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
