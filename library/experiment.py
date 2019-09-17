#!/usr/bin/env python

import datetime
import os
import pandas as pd
import PyQt5.QtWidgets as QtWidgets
import PyQt5.QtGui as QtGui
import PyQt5.QtCore as QtCore
import PyQt5.Qt as Qt

from . import database
from . import loader
from . import outline
from . import assignment
from . import analysis

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

        if os.path.exists(self._data.image_path):
            self.image_loader = loader.ImageLoader(self._data.image_path)
        else:
            # get new image
            path = QtWidgets.QFileDialog.getOpenFileName(
                self.window,
                "Replace image for {0}".format(self._data.image_path),
                os.getcwd(),
            )[0]
            try:
                self.image_loader = loader.ImageLoader(path)
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
        for elem_num, (label_str, kw) in enumerate([
            ("Date", "date"),
            ("Strain", "strain"),
            ("Medium", "medium"),
            ("Image", "image_path"),
        ]):
            label = QtWidgets.QLabel(label_str)
            label.setFont(label_font)
            val = QtWidgets.QLabel(self._data[kw])
            layout.addWidget(label, row_num, col_num)
            layout.addWidget(val, row_num, col_num + 1)
            col_num += 2
            if elem_num % 2 == 1:
                row_num += 1
                col_num = 0

        details_box.setLayout(layout)
        self.main_layout.addWidget(details_box)

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
        assigner = assignment.Assigner(self._data, self.image_loader)
        desktop = QtWidgets.QDesktopWidget()
        assigner.set_screen_res(
            desktop.width(),
            desktop.height(),
            desktop.logicalDpiX(),
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
        self._assignDimensions()
        if self.settings["num_frames"] == 1:
            self.setImageMode("static", True)
        elif self.settings["num_frames"] > 1:
            self.setImageMode("movie", True)

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

    def _assignDimensions(self):
        image_path = self.settings["image_path"]
        if os.path.exists(image_path):
            meta = loader.ImageLoader(image_path).im_metadata
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

    def setImageMode(self, imagemode, state):
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
            if self.settings["image_path"]:
                self._assignDimensions()

        elif movie_test and self.settings["image_mode"] != "movie":
            self.settings["image_mode"] = "movie"
            self.dim_widgets["frames"][0].setVisible(True)
            self.dim_widgets["frames"][1].setVisible(True)
            self.image_mode_buttons[0].setChecked(True)
            if self.settings["image_path"]:
                self._assignDimensions()

    def setChannelGreen(self, state):
        self.settings["channels"]["green"] = state

    def setChannelRed(self, state):
        self.settings["channels"]["red"] = state

    def confirm_settings(self):
        alert_flag = False
        for k, v in self.settings.items():
            if v is None:
                alert_flag = True
                form_elements = getattr(self, "{0}_form".format(k))
                form_elements[1].setStyleSheet("QLineEdit { border: 1px solid red }")

        if not alert_flag:
            self.save_settings()

    def save_settings(self):
        dupl = database.checkExperimentDuplicate(
            *self.settings["date"],
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
        )
        print("Experiment ID={0} overwritten".format(old_id))
        self.window.close()

    def write_settings(self):
        new_id = database.insertExperiment(
            *self.settings["date"],
            self.settings["medium"],
            self.settings["strain"],
            self.settings["image_path"],
            self.settings["image_mode"],
            self.settings["num_channels"],
            self.settings["num_slices"],
            self.settings["num_frames"],
        )
        print("Experiment saved with ID={0}".format(new_id))
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
        }


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
        widget = QtWidgets.QLineEdit(*widget_args, **widget_kwargs)
        widget.textChanged[str].connect(_pathValidate)

        if change_callback and type(change_callback) is list:
            for cb in change_callback:
                widget.textChanged[str].connect(cb)
        elif change_callback:
            widget.textChanged[str].connect(change_callback)

        def callback():
            path = QtWidgets.QFileDialog.getOpenFileName(
                self.window,
                "Image file",
                initial_path or os.getcwd(),
            )[0]
            widget.setText(path)

        btn = QtWidgets.QPushButton("Browse")
        btn.clicked.connect(callback)

        layout.addWidget(label, self.current_row, 0)
        layout.addWidget(widget, self.current_row, 1)
        layout.addWidget(btn, self.current_row, 2)
        self.current_row += 1
        return label, widget
