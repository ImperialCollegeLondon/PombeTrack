#!/usr/bin/env python

import datetime
import os
import PyQt5.QtWidgets as QtWidgets 
import PyQt5.QtGui as QtGui
import PyQt5.QtCore as QtCore
import PyQt5.Qt as Qt

from . import database

class ExperimentView:
    def __init__(self, experiment_id):
        self._data = database.getExperimentById(experiment_id)
        self.data_changed = False

    def show_experiment(self, window=None):
        if not window:
            self.app = QtWidgets.QApplication([])
            self.app.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
            font = QtGui.QFont("Fira Sans", 11)
            self.app.setFont(font)
            self.window = QtWidgets.QWidget()
        else:
            self.window = window

        self.window.setGeometry(0, 0, 800, 100)
        self.window.setWindowTitle("Experiment #{0}".format(self._data["experiment_id"]))

        main_layout = QtWidgets.QVBoxLayout()
        self._addDetails(main_layout)
        self._addOutline(main_layout)
        self._addLineageVerification(main_layout)
        self._addAnalysis(main_layout)
        self.window.setLayout(main_layout)
        self.window.show()
        if not window:
            self.app.exec_()

    def _addDetails(self, main_layout):
        details_box = QtWidgets.QGroupBox("Details")
        details_layout = QtWidgets.QGridLayout()

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
            details_layout.addWidget(label, row_num, col_num)
            details_layout.addWidget(val, row_num, col_num + 1)
            col_num += 2
            if elem_num % 2 == 1:
                row_num += 1
                col_num = 0

        details_box.setLayout(details_layout)
        main_layout.addWidget(details_box)

    def _addOutline(self, main_layout):
        pass

    def _addLineageVerification(self, main_layout):
        pass

    def _addAnalysis(self, main_layout):
        pass


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
            channel_green=self.settings["channels"]["green"],
            channel_red=self.settings["channels"]["red"],
        )
        print("Experiment ID={0} overwritten".format(old_id))
        self.window.close()

    def write_settings(self):
        new_id = database.insertExperiment(
            *self.settings["date"],
            self.settings["medium"],
            self.settings["strain"],
            self.settings["image_path"],
            self.settings["channels"]["green"],
            self.settings["channels"]["red"],
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
            self.channels_form[1][0].setChecked(True)
            self.channels_form[1][1].setChecked(True)

        self.settings = {
            "date": (datetime.datetime.today().year,
                     datetime.datetime.today().month,
                     datetime.datetime.today().day),
            "medium": None,
            "strain": None,
            "image_path": None,
            "channels": {"red": True, "green": True},
        }


    def create_new_experiment(self, window=None):
        """ This function will ask the user to define various parameters, and store them in a config file
        
        Needs to define:
            medium
            strain
            image file location
            channels

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

        self.window.setGeometry(0, 0, 800, 100)
        self.window.setWindowTitle("Add new experiment")

        main_layout = QtWidgets.QVBoxLayout()

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

        label = QtWidgets.QLabel("Fluorescent channels")
        btngroup = QtWidgets.QHBoxLayout()
        green = QtWidgets.QPushButton("Green")
        green.setCheckable(True)
        green.setChecked(True)
        green.clicked[bool].connect(self.setChannelGreen)
        red = QtWidgets.QPushButton("Red")
        red.setCheckable(True)
        red.setChecked(True)
        red.clicked[bool].connect(self.setChannelRed)
        btngroup.addWidget(green)
        btngroup.addWidget(red)
        layout.addWidget(label, self.current_row, 0)
        layout.addLayout(btngroup, self.current_row, 1)
        self.channels_form = label, (green, red)
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
