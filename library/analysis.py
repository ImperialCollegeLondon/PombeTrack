#!/usr/bin/env python3

import PyQt5.QtWidgets as QtWidgets
import PyQt5.QtGui as QtGui
import PyQt5.QtCore as QtCore
import PyQt5.Qt as Qt

from . import database

class Analyser:
    def __init__(self, experiment_data, image_loader):
        self.experiment_data = experiment_data
        self.image_loader = image_loader

    def set_screen_res(self, max_width_px, max_height_px, screen_dpi):
        self.max_width_px = max_width_px
        self.max_height_px = max_height_px
        self.screen_dpi = screen_dpi

    def start_analysing(self, parent_window):
        warning = QtWidgets.QMessageBox(parent_window)
        warning.setTextFormat(QtCore.Qt.RichText)
        warning.setWindowTitle("Not implemented")
        warning.setText(
            "Analysing hasn't been written quite yet"
        )
        warning.exec_()
