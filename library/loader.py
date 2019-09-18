#!/usr/bin/env python3

import numpy as np
import os
import tifffile

class ImageLoaderMulti:
    """ This class stores images in memory with image series """
    def __init__(self, imagepaths):
        self.paths = [x["image_path"] for x in imagepaths]
        self.im_preload = {}
        self.load_metadata()
        self.load_frame(0, 0)
        self.green_factor = np.load(os.path.join("resources", "green_correction.npy"))
        self.red_factor = np.load(os.path.join("resources", "red_correction.npy"))

    def load_metadata(self):
        with tifffile.TiffFile(self.paths[0]) as im_frames:
            if not im_frames.is_imagej:
                raise IOError("Not an ImageJ file, ask Miles about this?")

            self.im_metadata = im_frames.imagej_metadata
            self.num_channels = 1
            self.num_slices = 1
            if "frames" in self.im_metadata:
                raise Exception("Image file {0} has more than 1 frame".format(
                    self.paths[0],
                ))

            self.num_frames = len(self.paths)
            if "channels" in self.im_metadata:
                self.num_channels  = int(self.im_metadata["channels"])

            if "slices" in self.im_metadata:
               self.num_slices = int(self.im_metadata["slices"])

    def get_pixel_conversion(self):
        with tifffile.TiffFile(self.paths[0]) as im_frames:
            try:
                xres = im_frames.pages[0].tags["XResolution"].value
            except KeyError:
                return None

        px_um = xres[1] / xres[0]
        return px_um

    def load_frame(self, frame_idx, slice_idx=0, channel_idx=0):
        page_idx = ((frame_idx * self.num_channels * self.num_slices) +
                    (slice_idx * self.num_channels) +
                    (channel_idx))

        if page_idx not in self.im_preload:
            if channel_idx == 1:
                self.load_slice(frame_idx, slice_idx, channel_idx, page_idx, self.green_factor)
            elif channel_idx == 2:
                self.load_slice(frame_idx, slice_idx, channel_idx, page_idx, self.red_factor)
            else:
                self.load_slice(frame_idx, slice_idx, channel_idx, page_idx)

        return self.im_preload[page_idx]

    def load_slice(self, frame_idx, slice_idx, channel_idx, page_idx, factor=None):
        with tifffile.TiffFile(self.paths[frame_idx]) as im_frames:
            subpage_idx = (slice_idx * self.num_channels) + channel_idx
            page = im_frames.pages[subpage_idx]
            im = page.asarray()
            if factor is not None:
                im = im * factor

            self.im_preload[page_idx] = im


class ImageLoaderSingle:
    """ This class help stores images in memory """
    def __init__(self, image_path):
        self.path = image_path
        self.im_preload = {}
        self.load_metadata()
        self.load_frame(0, 0)
        self.green_factor = np.load(os.path.join("resources", "green_correction.npy"))
        self.red_factor = np.load(os.path.join("resources", "red_correction.npy"))

    def load_metadata(self):
        with tifffile.TiffFile(self.path) as im_frames:
            if not im_frames.is_imagej:
                raise IOError("Not an ImageJ file, ask Miles about this?")

            self.im_metadata = im_frames.imagej_metadata
            self.num_frames = 1
            self.num_channels = 1
            self.num_slices = 1
            if "frames" in self.im_metadata:
                self.num_frames = int(self.im_metadata["frames"])

            if "channels" in self.im_metadata:
                self.num_channels = int(self.im_metadata["channels"])

            if "slices" in self.im_metadata:
                self.num_slices = int(self.im_metadata["slices"])

            if "slices" in self.im_metadata:
                self.num_frames = int(self.im_metadata["slices"])

    def get_pixel_conversion(self):
        with tifffile.TiffFile(self.path) as im_frames:
            try:
                xres = im_frames.pages[0].tags["XResolution"].value
            except KeyError:
                return None

        px_um = xres[1] / xres[0]
        return px_um

    def load_frame(self, frame_idx, slice_idx=0, channel_idx=0):
        page_idx = ((frame_idx * self.num_channels * self.num_slices) +
                    (slice_idx * self.num_channels) +
                    (channel_idx))

        if page_idx not in self.im_preload:
            if channel_idx == 1:
                self.load_slice(page_idx, self.green_factor)
            elif channel_idx == 2:
                self.load_slice(page_idx, self.red_factor)
            else:
                self.load_slice(page_idx)

        return self.im_preload[page_idx]

    def load_slice(self, page_idx, factor=None):
        with tifffile.TiffFile(self.path) as im_frames:
            page = im_frames.pages[page_idx]
            im = page.asarray()
            if factor is not None:
                im = im * factor

        self.im_preload[page_idx] = im
