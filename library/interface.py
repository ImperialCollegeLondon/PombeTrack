#!/usr/bin/env python

import os
import pandas as pd
import PyQt5.QtWidgets as QtWidgets
import PyQt5.QtCore as QtCore
import PyQt5.QtGui as QtGui

from . import experiment
from . import database


class Interface:
    def __init__(self):
        self.app = QtWidgets.QApplication([])
        self.app.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
        main_font = QtGui.QFont("Fira Sans", 11)
        self.app.setFont(main_font)
        self.window = QtWidgets.QWidget()

        screen = self.app.primaryScreen()
        screen_size = screen.availableGeometry()
        self.window.setGeometry(0, 60, screen_size.width() / 2, screen_size.height() - 120)
        self.window.setWindowTitle("Static image analyser")
        self.window.setWindowIcon(QtGui.QIcon(os.path.abspath("resources/icon.png")))
        self.base_layout = QtWidgets.QVBoxLayout()
        self.base_layout.setAlignment(QtCore.Qt.AlignTop)

        self.add_menu()
        self.decorate_window()
        self.window.setLayout(self.base_layout)

    def add_menu(self):
        menubar = QtWidgets.QMenuBar(self.window)
        file_menu = menubar.addMenu("&File")

        refresh_action = QtWidgets.QAction("&Refresh", menubar)
        refresh_action.triggered[bool].connect(lambda: self._refresh_layout())
        file_menu.addAction(refresh_action)

        quit_action = QtWidgets.QAction("&Quit", menubar)
        quit_action.triggered[bool].connect(lambda: self.window.close())
        file_menu.addAction(quit_action)

        self.base_layout.setMenuBar(menubar)

    def initiate(self):
        self.window.show()
        self.app.exec_()

    def decorate_window(self):
        header_label = QtWidgets.QLabel("Experiments:")
        header_font = QtGui.QFont("Fira Sans", 13)
        header_font.setBold(True)
        header_label.setFont(header_font)
        self.base_layout.addWidget(header_label)

        existing_experiments = self.get_existing_experiments()
        if len(existing_experiments) > 0:
            self.experiment_table = QtWidgets.QTableWidget()
            self.experiment_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            self.experiment_table.setRowCount(len(existing_experiments))
            header_labels = [
                "Date", "Medium", "Strain", "Outlined", "Verified", "Analysed", "Cells detected",
            ]
            self.experiment_table.setColumnCount(len(header_labels))
            self.experiment_table.setHorizontalHeaderLabels(header_labels)
            table_rows = []
            for experiment in existing_experiments:
                date_item = QtWidgets.QTableWidgetItem()
                date_item.setText(experiment.date)
                date_item.setTextAlignment(QtCore.Qt.AlignCenter)
                date_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
                date_item._item_type = "experiment_date"
                date_item._experiment_id = experiment.experiment_id

                medium_item = QtWidgets.QTableWidgetItem()
                medium_item.setText(experiment.medium)
                medium_item.setTextAlignment(QtCore.Qt.AlignCenter)
                medium_item._item_type = "experiment_medium"
                medium_item._experiment_id = experiment.experiment_id

                strain_item = QtWidgets.QTableWidgetItem()
                strain_item.setText(experiment.strain)
                strain_item.setTextAlignment(QtCore.Qt.AlignCenter)
                strain_item._item_type = "experiment_strain"
                strain_item._experiment_id = experiment.experiment_id

                outlined_item = QtWidgets.QTableWidgetItem()
                outlined_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)

                if experiment.outlined:
                    outlined_item.setData(
                        QtCore.Qt.BackgroundRole,
                        QtGui.QBrush(QtGui.QColor("green"))
                    )
                    num_outlines = len(database.getOutlinesByExperimentId(experiment.experiment_id))
                    outlined_item.setText(str(num_outlines))

                else:
                    outlined_item.setData(
                        QtCore.Qt.BackgroundRole,
                        QtGui.QBrush(QtGui.QColor("red"))
                    )
                    outlined_item.setText("NO")

                outlined_item.setTextAlignment(QtCore.Qt.AlignCenter)
                outlined_item._experiment_id = experiment.experiment_id
                outlined_item._item_type = "experiment_outlined"

                verified_item = QtWidgets.QTableWidgetItem()
                verified_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
                if experiment.verified:
                    verified_item.setData(
                        QtCore.Qt.BackgroundRole,
                        QtGui.QBrush(QtGui.QColor("green"))
                    )
                    verified_item.setText("YES")

                else:
                    verified_item.setData(
                        QtCore.Qt.BackgroundRole,
                        QtGui.QBrush(QtGui.QColor("red"))
                    )
                    verified_item.setText("NO")

                verified_item.setTextAlignment(QtCore.Qt.AlignCenter)
                verified_item._experiment_id = experiment.experiment_id
                verified_item._item_type = "experiment_verified"

                cells_counted_item = QtWidgets.QTableWidgetItem()
                cells_counted_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
                analysed_item = QtWidgets.QTableWidgetItem()
                analysed_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
                if experiment.analysed:
                    analysed_item.setData(
                        QtCore.Qt.BackgroundRole,
                        QtGui.QBrush(QtGui.QColor("green"))
                    )
                    analysed_item.setText("YES")
                else:
                    analysed_item.setData(
                        QtCore.Qt.BackgroundRole,
                        QtGui.QBrush(QtGui.QColor("red"))
                    )
                    analysed_item.setText("NO")
                analysed_item.setTextAlignment(QtCore.Qt.AlignCenter)
                analysed_item._experiment_id = experiment.experiment_id
                analysed_item._item_type = "experiment_analysed"

                if experiment.verified:
                    cells_counted_item.setText(self.get_cell_count(experiment.experiment_id))
                else:
                    cells_counted_item.setText("0")

                cells_counted_item.setTextAlignment(QtCore.Qt.AlignCenter)
                cells_counted_item._experiment_id = experiment.experiment_id
                cells_counted_item._item_type = "experiment_cells_counted"

                table_rows.append({
                    "date": (experiment.date, date_item),
                    "medium": (experiment.medium, medium_item),
                    "strain": (experiment.strain, strain_item),
                    "outlined": (None, outlined_item),
                    "verified": (None, verified_item),
                    "analysed": (None, analysed_item),
                    "counted": (None, cells_counted_item),
                })

            table_rows.sort(key=lambda x: (x["date"][0], x["strain"][0]))

            for row_num, table_row in enumerate(table_rows):
                for col_num, table_col in enumerate(["date", "medium", "strain", "outlined", "verified", "analysed", "counted"]):
                    self.experiment_table.setItem(row_num, col_num, table_row[table_col][1])

            self.experiment_table.itemClicked.connect(self.table_click_event)
            self.experiment_table.itemDoubleClicked.connect(self.view_experiment)
            self.experiment_table.itemChanged.connect(self.table_change_event)
            self.experiment_table.setCurrentCell(0, 0, QtCore.QItemSelectionModel.Clear)
            self.base_layout.addWidget(self.experiment_table)
        else:
            self.experiment_table = QtWidgets.QLabel("No experiments")
            self.base_layout.addWidget(self.experiment_table)

        self.btn_row = QtWidgets.QHBoxLayout()
        btn = QtWidgets.QPushButton("Add new experiment")
        btn.clicked[bool].connect(self.add_experiment)
        self.btn_row.addWidget(btn)

        self.edit_btn = QtWidgets.QPushButton("View experiment")
        self.edit_btn.clicked[bool].connect(self.view_experiment)
        self.btn_row.addWidget(self.edit_btn)

        self.btn_row.setAlignment(QtCore.Qt.AlignLeft)
        self.base_layout.addLayout(self.btn_row)

    def get_cell_count(self, experiment_id):
        cells = database.getCellsByExperimentId(experiment_id)
        cells_wanted = database.getCellsByExperimentId(experiment_id, True, True)
        return "{0} ({1})".format(len(cells_wanted), len(cells))

    def table_change_event(self, item):
        if not hasattr(item, "_item_type"):
            return

        permitted_types = ["experiment_medium", "experiment_strain"]
        if item._item_type not in permitted_types:
            return

        data = {}
        item_name = item._item_type.split("experiment_")[1]
        data[item_name] = item.data(QtCore.Qt.DisplayRole)
        database.updateExperimentById(
            item._experiment_id,
            **data
        )
        self._refresh_layout()

    def table_click_event(self, item):
        self.edit_btn._experiment_id = item._experiment_id

    def get_existing_experiments(self):
        experiments = database.getExperiments()
        return experiments

    def add_experiment(self):
        e = experiment.Experiment()
        dialog = QtWidgets.QDialog(self.window)
        dialog.setModal(True)
        dialog.finished[int].connect(self._refresh_layout)
        e.create_new_experiment(window=dialog)

    def view_experiment(self, item=None):
        forbidden_types = ["experiment_medium", "experiment_strain"]
        if item and item._item_type in forbidden_types:
            return

        if (not hasattr(self.experiment_table, "selectedItems") or
                len(self.experiment_table.selectedItems()) == 0):
            return

        e = experiment.ExperimentView(self.edit_btn._experiment_id)
        dialog = QtWidgets.QDialog(self.window)
        dialog.setModal(True)
        dialog.finished[int].connect(self._refresh_layout)
        e.show_experiment(window=dialog)

    def _clear_layout(self, l):
        for i in reversed(range(l.count())):
            w = l.takeAt(i)
            if type(w) is QtWidgets.QHBoxLayout or type(w) is QtWidgets.QVBoxLayout:
                self._clear_layout(w)
            else:
                w.widget().deleteLater()

    def _refresh_layout(self, *args):
        # refresh table
        self._clear_layout(self.base_layout)
        self.decorate_window()


