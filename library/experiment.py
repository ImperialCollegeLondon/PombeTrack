#!/usr/bin/env python

import datetime
import matplotlib
import numpy as np
import os
import pandas as pd
import PyQt5.QtWidgets as QtWidgets
import PyQt5.QtGui as QtGui
import PyQt5.QtCore as QtCore
import PyQt5.Qt as Qt
import time
import uuid

from . import database
from . import loader
from . import outline
from . import assignment
from . import analysis

from . import segmentation

class ExperimentView:
    def __init__(self, experiment_id):
        self._data = database.getExperimentById(experiment_id)

    def show_experiment(self, window=None):
        if not window:
            self.app = QtWidgets.QApplication([])
            self.app.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
            font = QtGui.QFont("Fira Sans", 11)
            self.app.setFont(font)
            self.window = QtWidgets.QWidget()
        else:
            self.window = window

        self.window.setGeometry(0, 60, 800, 100)
        self.window.setWindowTitle("Experiment #{0}".format(self._data.experiment_num))

        if self._data.file_mode == "single" and os.path.exists(self._data.image_path):
            self.image_loader = loader.ImageLoaderSingle(self._data.image_path)
        elif self._data.file_mode == "single":
            # get new image
            path = QtWidgets.QFileDialog.getOpenFileName(
                self.window,
                "Replace image for {0}".format(self._data.image_path),
                os.getcwd(),
            )[0]
            try:
                self.image_loader = loader.ImageLoaderSingle(path)
            except:
                alert = QtWidgets.QMessageBox()
                alert.warning(
                    self.window,
                    "Error",
                    "Uh oh... image file ({0}) cannot be loaded".format(
                        repr(path),
                    ),
                )
                return
            else:
                database.updateExperimentById(
                    self._data.experiment_id,
                    image_path=path,
                )
        elif self._data.file_mode == "multi":
            paths = database.getImagePaths(self._data.experiment_id)
            self.image_loader = loader.ImageLoaderMulti(paths)

        self.main_layout = QtWidgets.QVBoxLayout()
        self.createLayout()
        self.window.setLayout(self.main_layout)
        self.window.show()
        if not window:
            self.app.exec_()

    def createLayout(self):
        menubar = QtWidgets.QMenuBar()
        file_menu = menubar.addMenu("&File")
        delete_action = QtWidgets.QAction("&Delete experiment", menubar)
        delete_action.triggered[bool].connect(lambda: self.delete_experiment())
        file_menu.addAction(delete_action)

        quit_action = QtWidgets.QAction("&Close", menubar)
        quit_action.triggered[bool].connect(lambda: self.window.close())
        file_menu.addAction(quit_action)
        self.main_layout.setMenuBar(menubar)

        self._addDetails()
        self._addBatchAuto()
        self._addOutline()
        self._addLineageVerification()
        self._addAnalysis()

    def delete_experiment(self, close_window=True):
        # confirm first
        alert = QtWidgets.QMessageBox()
        delete_confirm = alert.question(
            self.window,
            "Delete experiment?",
            "Are you really sure you want to delete this experiment?"
        )
        if delete_confirm == QtWidgets.QMessageBox.Yes:
            database.deleteExperimentById(self._data.experiment_id)
            if close_window:
                self.window.close()

    def _addDetails(self):
        details_box = QtWidgets.QGroupBox("Details")
        layout = QtWidgets.QGridLayout()

        label_font = QtGui.QFont("Fira Sans", 11)
        label_font.setBold(True)

        row_num = 0
        col_num = 0
        detail_items = [
            ("Date", "date"),
            ("Strain", "strain"),
            ("Medium", "medium"),
            ("Image", "image_path"),
            ("Mode", "image_mode"),
            ("Files", "file_mode"),
            ("Frames", "num_frames"),
            ("Channels", "num_channels"),
            # ("Slices", "num_slices"),
        ]
        for elem_num, (label_str, kw) in enumerate(detail_items):
            label = QtWidgets.QLabel(label_str)
            label.setFont(label_font)
            val = QtWidgets.QLabel(str(self._data[kw]))
            layout.addWidget(label, row_num, col_num)
            layout.addWidget(val, row_num, col_num + 1)
            col_num += 2
            if elem_num % 2 == 1:
                row_num += 1
                col_num = 0

        details_box.setLayout(layout)
        self.main_layout.addWidget(details_box)


    # Indicate progress
    def set_status(self, text=None, clear=False, status=None):
        if not clear and text is not None:
            self.status_bar.showMessage(text)
        else:
            self.status_bar.clearMessage()

        if not status and self.current_status:
            QtGui.QGuiApplication.restoreOverrideCursor()
            # QtGui.QGuiApplication.setOverrideCursor(QtGui.QCursor(
            #     QtCore.Qt.ArrowCursor
            # ))
        elif status == "working" and self.current_status != status:
            QtGui.QGuiApplication.setOverrideCursor(QtGui.QCursor(
                QtCore.Qt.WaitCursor
            ))

        self.current_status = status
        self.status_bar.repaint()


    def _addBatchAuto(self):
        self.current_status = None
        auto_box = QtWidgets.QGroupBox("Automatic segmentation for all")
        layout = QtWidgets.QVBoxLayout()

        # Box for progression status
        status_label = QtWidgets.QLabel("Status:")
        self.status_bar = QtWidgets.QStatusBar()
        status_layout = QtWidgets.QHBoxLayout()
        status_layout.setAlignment(QtCore.Qt.AlignLeft)
        status_layout.addWidget(status_label)
        status_layout.addWidget(self.status_bar)

        # Box for button
        auto_btn = QtWidgets.QPushButton("Start Auto")
        auto_btn.clicked.connect(lambda: self.batch_auto())
        layout.addWidget(auto_btn)
        layout.addLayout(status_layout)

        auto_box.setLayout(layout)
        self.main_layout.addWidget(auto_box)

    def batch_auto(self):
        self.outline_store = os.path.join(
            "data", "outlines", self._data.experiment_id
        )
        if not os.path.exists(self.outline_store):
            os.makedirs(self.outline_store)
        startt=time.time()
        for ff in range(0,self.image_loader.num_frames):
            self.current_frame_idx=ff
            self.auto_one()
        self.set_status(text="Auto segmentation finished")
        print(time.time()-startt)

    def auto_one(self):
        self.set_status(
        text="Preprocessing image {0} of {1}".format(
            self.current_frame_idx+1, self.image_loader.num_frames), status="working"
            )
        # load_frame: frame, z-slice, channel
        im_mid = self.image_loader.load_frame(self.current_frame_idx, int(np.floor(self.image_loader.num_slices / 2)), 0)
        im_up = self.image_loader.load_frame(self.current_frame_idx, int(np.floor(self.image_loader.num_slices / 2) - 1), 0)
        im = np.maximum(im_mid, im_up)


        im_pp = segmentation.preprocessing(im)
        im_i = segmentation.find_cellinterior(im_pp)
        im_wat = segmentation.find_watershed(im_i)
        #  bd = segmentation.find_bd(im_wat)


        for index in range(1, im_wat.max()+1):
            self.set_status(
                    text="Processing image {0} of {1}, cell {2} of {3}".format(
                        self.current_frame_idx+1, self.image_loader.num_frames,
                        index, im_wat.max()), status="working"
                    )
            im_ii = im_wat == index
            if (np.any(np.asarray(im_ii.nonzero()) == 0) or
                    np.any(np.asarray(im_ii.nonzero()) == 2047)):
                continue

            im_ii_bd = segmentation.find_boundaries(im_ii, mode = 'inner')
            #  Sort in radial
            bd_ii_sorted = segmentation.sort_in_order(im_ii_bd)

            #  Define the balloon object
            balloon_obj, origin_y, origin_x, halfwidth = segmentation.find_balloon_obj(bd_ii_sorted.astype(int)[::5], im)


            #  Test if the cell exists
            overlap = False
            outline_data = database.getOutlinesByFrameIdx(self.current_frame_idx, self._data.experiment_id)
            for i, outline in enumerate(outline_data):
                if not os.path.exists(outline.coords_path):
                    continue
                if outline.centre_x not in range(origin_x, origin_x + 2 * halfwidth) or outline.centre_y not in range(origin_y, origin_y + 2 * halfwidth):
                    continue

                polygonpath = matplotlib.path.Path(np.append(balloon_obj.get_coordinates(accept = True),\
                        balloon_obj.get_coordinates(accept = True)[1, :].reshape(1, 2), axis = 0), closed = True)
                if polygonpath.contains_point([outline.centre_y-origin_y, outline.centre_x - origin_x]):
                    overlap = True

            if overlap:
                continue




            # Evolve the contour
            try:
                sensitivity = 0.4
                area_init = balloon_obj.get_area()
                for i in range(20):
                    balloon_obj.evolve(display = False, image_percentile = sensitivity)
                    if balloon_obj.get_area() > 1.5 * area_init or balloon_obj.get_area() < 0.5 * area_init:
                        raise ValueError()
                full_coords = balloon_obj.get_coordinates(accept = True) + [origin_y, origin_x]
                outline_id  =  str(uuid.uuid4())
                cell_id  =  str(uuid.uuid4())
                centre = [np.mean(full_coords[:, 0]).astype(int), np.mean(full_coords[:, 1]).astype(int)]
                centre_y, centre_x = centre

                # Save to database
                offset_left=0
                offset_top=0
                coords_path = os.path.join(
                    self.outline_store,
                    "{0}.npy".format(outline_id)
                )
                data = {
                    "outline_id": outline_id,
                    "cell_id": cell_id,
                    "experiment_num": self._data.experiment_num,
                    "experiment_id": self._data.experiment_id,
                    "image_path": self._data.image_path,
                    "frame_idx": self.current_frame_idx,
                    "coords_path": coords_path,
                    "offset_left": offset_left,
                    "offset_top": offset_top,
                    "parent_id": "",
                    "centre_x":int(centre_x),
                    "centre_y":int(centre_y),
                }
                #  if self._data.image_mode == "static":
                    #  data["parent_id"] = ""
    #
                if os.path.exists(coords_path):
                    os.remove(coords_path)
                    database.updateOutlineById(
                        outline_id,
                        centre_x=int(centre_x),
                        centre_y=int(centre_y),
                    )
                else:
                    database.insertOutline(**data)

                np.save(coords_path, full_coords)

                #  if self.previous_id:
                    #  database.addOutlineChild(self.previous_id, child1=self.outline_id)

                database.updateExperimentById(
                    self._data.experiment_id,
                    verified=False,
                )
                #  database.deleteCellById(self.cell_id)
            except ValueError:
                continue





    def _addOutline(self):
        outline_box = QtWidgets.QGroupBox("Outlines")
        layout = QtWidgets.QVBoxLayout()
        if self._data.outlined:
            outlines = database.getOutlinesByExperimentId(self._data.experiment_id)
            num_outlines = len(outlines)
            num_cells = len(pd.DataFrame(outlines).cell_id.unique())
            label_str = "{0} outlines have been defined arranged as {1} cells{2}".format(
                num_outlines,
                num_cells,
                self._data.verified and " " or " (unverified)"
            )
            label = QtWidgets.QLabel(label_str)
        else:
            label = QtWidgets.QLabel("No outlines have been created.")

        layout.addWidget(label)

        outline_btn = QtWidgets.QPushButton("Outline Cells")
        outline_btn.clicked.connect(lambda: self.outline_cells())
        layout.addWidget(outline_btn)

        outline_box.setLayout(layout)
        self.main_layout.addWidget(outline_box)

    def outline_cells(self):
        outliner = outline.Outliner(self._data, self.image_loader)
        desktop = QtWidgets.QDesktopWidget()
        outliner.set_screen_res(
            desktop.width(),
            desktop.height(),
            desktop.logicalDpiX(),
        )
        outliner.start_outlining(self.window)
        outliner.window.finished[int].connect(self.outline_finished)

    def outline_finished(self, *args):
        outlines = database.getOutlinesByExperimentId(self._data.experiment_id)
        if len(outlines) > 0:
            database.updateExperimentById(self._data.experiment_id, outlined=True)
        else:
            database.updateExperimentById(self._data.experiment_id, outlined=False)

        self._data = database.getExperimentById(self._data.experiment_id)
        self._refreshLayout()

    def _addLineageVerification(self):
        box = QtWidgets.QGroupBox("Lineage verification")
        layout = QtWidgets.QVBoxLayout()
        if self._data.verified:
            pass
        else:
            pass

        lineage_btn = QtWidgets.QPushButton("Verify Lineages")
        lineage_btn.clicked[bool].connect(lambda: self.assign_cell_lineages())
        layout.addWidget(lineage_btn)

        box.setLayout(layout)
        self.main_layout.addWidget(box)

    def assign_cell_lineages(self):
        desktop = QtWidgets.QDesktopWidget()
        assigner = assignment.Assigner(
            self._data, self.image_loader,
            (desktop.width(),
             desktop.height(),
             desktop.logicalDpiX()),
        )
        assigner.start_assigning(self.window)
        assigner.window.finished[int].connect(self.assignment_finished)

    def assignment_finished(self, *args):
        outlines = database.getOutlinesByExperimentId(self._data.experiment_id)
        if len(outlines) > 0:
            database.updateExperimentById(self._data.experiment_id, verified=True)
        else:
            database.updateExperimentById(self._data.experiment_id, verified=False)

        self._data = database.getExperimentById(self._data.experiment_id)
        self._refreshLayout()

    def _addAnalysis(self):
        if self._data.verified:
            box = QtWidgets.QGroupBox("Analysis")
            layout = QtWidgets.QVBoxLayout()
            self.analyser = analysis.Analyser(self, self._data, self.image_loader)
            desktop = QtWidgets.QDesktopWidget()
            self.analyser.set_screen_res(
                desktop.width(),
                desktop.height(),
                desktop.logicalDpiX(),
            )
            box.setLayout(layout)
            self.main_layout.addWidget(box)
            self.analyser.construct_box(layout)

    def _clearLayout(self, l):
        for i in reversed(range(l.count())):
            w = l.takeAt(i)
            if type(w) is QtWidgets.QHBoxLayout or type(w) is QtWidgets.QVBoxLayout:
                self._clear_layout(w)
            else:
                w.widget().deleteLater()

    def _refreshLayout(self, *args):
        self._clearLayout(self.main_layout)
        self.createLayout()


class Experiment:
    def __init__(self):
        self._check_database()
        self.reset_app()

    def _check_database(self):
        if not database.checkTable("experiments"):
            database.createExperimentsTable()

    def setDate(self, date):
        self.settings["date"] = (
            date.year(),
            date.month(),
            date.day()
        )

    def setMedium(self, medium):
        self.settings["medium"] = medium

    def setStrain(self, strain):
        self.settings["strain"] = strain

    def setImagePath(self, path):
        self.settings["image_path"] = path
        if not os.path.exists(path):
            self.settings["image_path"] = ""
            self.settings["num_frames"] = 1
            self.settings["num_channels"] = 1
            self.settings["num_slices"] = 1
            return

        if self.settings["file_mode"] == "multi":
            self._getPotentialFiles()

        self._assignDimensions()
        if self.settings["num_frames"] == 1:
            self.setImageMode("static", True, check=False)
        elif self.settings["num_frames"] > 1:
            self.setImageMode("movie", True, check=False)

    def _getNum(self, value, widget=None):
        try:
            num = int(value)
            assert num > 0
        except (ValueError, AssertionError):
            if widget:
                widget.setStyleSheet("QLineEdit { border: 1px solid red }")
            return
        else:
            if widget:
                widget.setStyleSheet("QLineEdit { }")
            return num

    def setNumF(self, value):
        widget = self.dim_widgets["frames"][1]
        num = self._getNum(value, widget)
        if not num:
            return

        widget.setText(str(num))
        self.settings["num_frames"] = num

    def setNumZ(self, value):
        widget = self.dim_widgets["zslices"][1]
        num = self._getNum(value, widget)
        if not num:
            return

        widget.setText(str(num))
        self.settings["num_slices"] = num

    def setNumC(self, value):
        widget = self.dim_widgets["channels"][1]
        num = self._getNum(value, widget)
        if not num:
            return

        widget.setText(str(num))
        self.settings["num_channels"] = num

    def _getPotentialFiles(self):
        # walk directory searching for TIFFs
        self.image_files = []
        image_dir = self.settings["image_path"]
        for root, dirnames, filenames in os.walk(image_dir):
            self.image_files.extend([
                os.path.join(root, f)
                for f in filenames
                if f.lower().endswith(".tif") or f.lower().endswith(".tiff")
            ])

    def _assignDimensions(self):
        if not self.settings["image_path"]:
            return

        if self.settings["file_mode"] == "single":
            image_path = self.settings["image_path"]
            if os.path.isdir(image_path):
                if len(self.image_files) > 0:
                    image_path = self.image_files[0]
                else:
                    image_path = ""
                self.image_path_form[1].setText(image_path)
                return

            if os.path.exists(image_path):
                meta = loader.ImageLoaderSingle(image_path).im_metadata
                if "frames" in meta:
                    self.setNumF(meta["frames"])
                else:
                    self.setNumF(1)

                if "channels" in meta:
                    self.setNumC(meta["channels"])
                else:
                    self.setNumC(1)

                if "slices" in meta:
                    self.setNumZ(meta["slices"])
                else:
                    self.setNumZ(1)
        else:
            image_dir = self.settings["image_path"]
            if not os.path.isdir(image_dir):
                image_dir = os.path.dirname(self.settings["image_path"])
                self.image_path_form[1].setText(image_dir)
                return

            if len(self.image_files) == 0:
                self.setNumC(0)
                self.setNumZ(0)
                self.setNumF(0)
            else:
                # use first image in set to define num_slices, num_channels
                meta = loader.ImageLoaderSingle(self.image_files[0]).im_metadata
                if "channels" in meta:
                    self.setNumC(meta["channels"])
                else:
                    self.setNumC(1)

                if "slices" in meta:
                    self.setNumZ(meta["slices"])
                else:
                    self.setNumZ(1)

                self.setNumF(len(self.image_files))

    def setImageMode(self, imagemode, state, check=True):
        static_test = (
            (imagemode == "static" and state == True) or
            (imagemode == "movie" and state == False)
        )
        movie_test = (
            (imagemode == "movie" and state == True) or
            (imagemode == "static" and state == False)
        )
        if static_test and self.settings["image_mode"] != "static":
            self.settings["image_mode"] = "static"
            # hide num frames input box
            self.dim_widgets["frames"][0].setVisible(False)
            self.dim_widgets["frames"][1].setVisible(False)
            self.image_mode_buttons[1].setChecked(True)
            if self.settings["image_path"] and check:
                self._assignDimensions()

        elif movie_test and self.settings["image_mode"] != "movie":
            self.settings["image_mode"] = "movie"
            self.dim_widgets["frames"][0].setVisible(True)
            self.dim_widgets["frames"][1].setVisible(True)
            self.image_mode_buttons[0].setChecked(True)
            if self.settings["image_path"] and check:
                self._assignDimensions()

    def setChannelGreen(self, state):
        self.settings["channels"]["green"] = state

    def setChannelRed(self, state):
        self.settings["channels"]["red"] = state

    def confirm_settings(self):
        alert_flag = False
        for k, v in self.settings.items():
            if v is None or (k.startswith("num_") and v == 0):
                alert_flag = True
                form_elements = getattr(self, "{0}_form".format(k))
                form_elements[1].setStyleSheet("QLineEdit { border: 1px solid red }")

        if not alert_flag:
            self.save_settings()

    def save_settings(self):
        dupl = database.checkExperimentDuplicate(
            self.settings["date"][0],
            self.settings["date"][1],
            self.settings["date"][2],
            self.settings["medium"],
            self.settings["strain"],
            self.settings["image_path"],
            )
        if dupl:
            alert = QtWidgets.QMessageBox()
            overwrite_confirm = alert.question(
                self.window,
                "Error",
                "An experiment with these settings has already been saved!\nOverwrite it?"
            )
            if overwrite_confirm == QtWidgets.QMessageBox.Yes:
                self.overwrite_settings(dupl)
            elif overwrite_confirm == QtWidgets.QMessageBox.No:
                return False
        else:
            self.write_settings()

    def overwrite_settings(self, old_id):
        # by definition, to overwrite, the date, medium, and strain must be the same
        database.updateExperimentById(
            old_id,
            image_path=self.settings["image_path"],
            image_mode=self.settings["image_mode"],
            num_channels=self.settings["num_channels"],
            num_slices=self.settings["num_slices"],
            num_frames=self.settings["num_frames"],
            file_mode=self.settings["file_mode"],
        )
        print("Experiment ID={0} overwritten".format(old_id))
        self.window.close()

    def write_settings(self):
        new_num, new_id = database.insertExperiment(
            *self.settings["date"],
            self.settings["medium"],
            self.settings["strain"],
            self.settings["image_path"],
            self.settings["image_mode"],
            self.settings["num_channels"],
            self.settings["num_slices"],
            self.settings["num_frames"],
            self.settings["file_mode"],
        )
        if self.settings["file_mode"] == "multi":
            if not database.checkTable("imagepath"):
                database.createImagePathTable()

            for i, image_path in enumerate(self.image_files):
                database.insertImagePath(
                    i + 1,
                    new_id,
                    image_path,
                )

        print("Experiment #{0} saved with ID={1}".format(new_num, new_id))
        self.window.close()

    def reset_app(self):
        if hasattr(self, "window"):
            self.date_form[1].setSelectedDate(
                QtCore.QDate(*self.settings["date"])
            )
            self.medium_form[1].setText("")
            self.medium_form[2].setCurrentIndex(0)
            self.strain_form[1].setText("")
            self.strain_form[2].setCurrentIndex(0)
            self.setImageMode("movie", True)
            # self.channels_form[1][0].setChecked(True)
            # self.channels_form[1][1].setChecked(True)

        self.settings = {
            "date": (datetime.datetime.today().year,
                     datetime.datetime.today().month,
                     datetime.datetime.today().day),
            "medium": None,
            "strain": None,
            "image_path": None,
            "image_mode": "movie",
            "num_channels": 1,
            "num_slices": 1,
            "num_frames": 1,
            "file_mode": "single",
        }
        self.image_files = []


    def create_new_experiment(self, window=None):
        """ This function will ask the user to define various parameters, and store them in a config file

        Needs to define:
            medium
            strain
            image file location
            image mode (static or movie)
            image dimensions (frames, slices, channels)

        Note, needs a ANALYSED boolean flag somewhere, or possibly a date stamp specifying when
        """
        if not window:
            self.app = QtWidgets.QApplication([])
            self.app.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
            font = QtGui.QFont("Fira Sans", 11)
            self.app.setFont(font)
            self.window = QtWidgets.QWidget()
        else:
            self.window = window

        self.window.setGeometry(0, 60, 800, 100)
        self.window.setWindowTitle("Add new experiment")

        main_layout = QtWidgets.QVBoxLayout()

        menubar = QtWidgets.QMenuBar(self.window)
        file_menu = menubar.addMenu("&File")
        quit_action = QtWidgets.QAction("&Close", menubar)
        quit_action.triggered[bool].connect(lambda: self.window.close())
        file_menu.addAction(quit_action)
        main_layout.setMenuBar(menubar)

        title_font = QtGui.QFont("Fira Sans", 20)
        title_font.setBold(True)
        title = QtWidgets.QLabel("Add new experiment")
        title.setFont(title_font)
        main_layout.addWidget(title)

        layout = QtWidgets.QGridLayout()
        # layout.setColumnStretch(0, 0)
        layout.setColumnStretch(1, 5)
        layout.setColumnStretch(2, 0)

        default_media = [
            "EMM",
            "EMMS",
            "YE",
            "YES",
        ]

        self.current_row = 0

        label = QtWidgets.QLabel("Date")
        today = datetime.datetime.today()
        date_widget = QtWidgets.QCalendarWidget()
        date_widget.clicked.connect(self.setDate)
        layout.addWidget(label, self.current_row, 0)
        layout.addWidget(date_widget, self.current_row, 1)
        self.date_form = label, date_widget
        self.current_row += 1

        self.medium_form = self._addFormEdit(
            layout,
            "Medium",
            change_callback=self.setMedium,
            options=default_media,
        )

        self.strain_form = self._addFormEdit(
            layout,
            "Strain",
            change_callback=self.setStrain,
        )

        self.image_path_form = self._addFormPath(
            layout,
            "Image data",
            change_callback=self.setImagePath,
            initial_path="/mnt/d",
        )

        label = QtWidgets.QLabel("Image mode")
        btngroup = QtWidgets.QHBoxLayout()
        movie = QtWidgets.QRadioButton("Movie")
        movie.setChecked(True)
        movie.clicked[bool].connect(lambda state: self.setImageMode("movie", state))
        static = QtWidgets.QRadioButton("Static")
        static.setChecked(False)
        static.clicked[bool].connect(lambda state: self.setImageMode("static", state))
        self.image_mode_buttons = [movie, static]
        btngroup.addWidget(movie)
        btngroup.addWidget(static)
        layout.addWidget(label, self.current_row, 0)
        layout.addLayout(btngroup, self.current_row, 1)
        self.current_row += 1

        label = QtWidgets.QLabel("Dimensionality")
        dim_layout = QtWidgets.QHBoxLayout()
        dim_frames = QtWidgets.QLabel("Frames")
        dim_frames_widget = QtWidgets.QLineEdit("1")
        dim_zslices = QtWidgets.QLabel("Z-slices")
        dim_zslices_widget = QtWidgets.QLineEdit("1")
        dim_channels = QtWidgets.QLabel("Channels")
        dim_channels_widget = QtWidgets.QLineEdit("1")
        self.dim_widgets = {
            "frames": (dim_frames, dim_frames_widget, self.setNumF),
            "zslices": (dim_zslices, dim_zslices_widget, self.setNumZ),
            "channels": (dim_channels, dim_channels_widget, self.setNumC),
        }
        for k, (wlabel, wwidg, wcb) in self.dim_widgets.items():
            wwidg.textChanged[str].connect(wcb)
            dim_layout.addWidget(wlabel)
            dim_layout.addWidget(wwidg)

        layout.addWidget(label, self.current_row, 0)
        layout.addLayout(dim_layout, self.current_row, 1)
        self.current_row += 1

        confirm_btn = QtWidgets.QPushButton("OK")
        confirm_btn.clicked.connect(self.confirm_settings)
        layout.addWidget(confirm_btn, self.current_row, 2)

        main_layout.addLayout(layout)
        self.window.setLayout(main_layout)
        self.window.show()
        if not window:
            self.app.exec_()

    def _addFormEdit(
        self, layout,
        label_text, label_args=[], label_kwargs={},
        widget_args=[], widget_kwargs={},
        options=None, options_callback=None,
        row_idx=None,
        change_callback=None,
    ):
        label = QtWidgets.QLabel(label_text, *label_args, **label_kwargs)
        widget = QtWidgets.QLineEdit(*widget_args, **widget_kwargs)

        if change_callback and type(change_callback) is list:
            for cb in change_callback:
                widget.textChanged[str].connect(cb)
        elif change_callback:
            widget.textChanged[str].connect(change_callback)

        def callback(index):
            if index == 0:
                widget.setText("")
            else:
                widget.setText(options[index - 1])
            widget.setStyleSheet("")

        if options:
            if not options_callback:
                options_callback = callback

            combo_box = QtWidgets.QComboBox()
            combo_box.addItem("--select--")
            for option in options:
                combo_box.addItem(option)

            combo_box.currentIndexChanged[int].connect(options_callback)
        else:
            combo_box = QtWidgets.QLabel("")

        layout.addWidget(label, self.current_row, 0)
        layout.addWidget(widget, self.current_row, 1)
        layout.addWidget(combo_box, self.current_row, 2)
        self.current_row += 1
        return label, widget, combo_box

    def _addFormPath(
        self, layout,
        label_text, label_args=[], label_kwargs={},
        widget_args=[], widget_kwargs={},
        initial_path=None,
        row_idx=None,
        change_callback=None,
    ):
        def _pathValidate(path):
            if os.path.exists(path):
                widget.setStyleSheet("QLineEdit { color: black }")
            else:
                widget.setStyleSheet("QLineEdit { color: rgb(157, 40, 20) }")

        label = QtWidgets.QLabel(label_text, *label_args, **label_kwargs)
        centralpart = QtWidgets.QGridLayout()
        centralpart.setColumnStretch(0, 5)
        centralpart.setColumnStretch(1, 1)

        widget = QtWidgets.QLineEdit(*widget_args, **widget_kwargs)
        widget.textChanged[str].connect(_pathValidate)
        checkbox = QtWidgets.QCheckBox("Image series?")
        def checkbox_cb(state):
            if state:
                self.settings["file_mode"] = "multi"
            else:
                self.settings["file_mode"] = "single"
            self._assignDimensions()

        checkbox.clicked[bool].connect(checkbox_cb)

        centralpart.addWidget(widget, 0, 0)
        centralpart.addWidget(checkbox, 0, 1)

        if change_callback and type(change_callback) is list:
            for cb in change_callback:
                widget.textChanged[str].connect(cb)
        elif change_callback:
            widget.textChanged[str].connect(change_callback)

        def callback():
            if self.settings["file_mode"] == "multi":
                path = QtWidgets.QFileDialog.getExistingDirectory(
                    self.window,
                    "Image root directory",
                    initial_path or os.getcwd(),
                )
                widget.setText(path)
            else:
                path = QtWidgets.QFileDialog.getOpenFileName(
                    self.window,
                    "Image file",
                    initial_path or os.getcwd(),
                )[0]
                widget.setText(path)

        btn = QtWidgets.QPushButton("Browse")
        btn.clicked.connect(callback)

        layout.addWidget(label, self.current_row, 0)
        layout.addLayout(centralpart, self.current_row, 1)
        layout.addWidget(btn, self.current_row, 2)
        self.current_row += 1
        return label, widget
