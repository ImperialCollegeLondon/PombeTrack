#!/usr/bin/env python3


import tifffile

class ImageLoader:
    """ This class help stores images in memory """
    def __init__(self, image_path):
        self.path = image_path
        self.im_preload = {}
        self.load_metadata()
        self.load_frame(0, 0)

    def load_metadata(self):
        with tifffile.TiffFile(self.path) as im_frames:
            if not im_frames.is_imagej:
                raise IOError("Not an ImageJ file, ask Miles about this?")

            self.im_metadata = im_frames.imagej_metadata
            self.num_frames = 1
            self.num_channels = 1
            if "frames" in self.im_metadata:
                self.num_frames = int(self.im_metadata["frames"])

            if "channels" in self.im_metadata:
                self.num_channels  = int(self.im_metadata["channels"])

    def load_frame(self, frame_idx, channel_num):
        slice_idx = (frame_idx * self.num_channels) + channel_num
        if slice_idx not in self.im_preload:
            self.load_slice(slice_idx)

        return self.im_preload[slice_idx]

    def load_slice(self, slice_idx):
        with tifffile.TiffFile(self.path) as im_frames:
            page = im_frames.pages[slice_idx]
            im = page.asarray()

        self.im_preload[slice_idx] = im
