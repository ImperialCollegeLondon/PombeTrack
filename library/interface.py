#!/usr/bin/env python

""" Provide the base interface for the entire application.

Creates the Qt5 application and all its immediate daughter windows.
Retrieves and displays a summarised view of all experiments, and provides
access to adding, viewing, and removing experiments.
"""

import os

import PyQt5.QtWidgets as QtWidgets
import PyQt5.QtCore as QtCore
import PyQt5.QtGui as QtGui

from . import experiment
from . import database


class Interface:
    """ Base interface for PombeTrack.

    Creates the first window of the interface, with the experiment table
    displaying existing experiments with a subset of their status/settings;
    buttons to add experiments, remove experiments, and view/edit experiments.

    Viewing/editing experiment functionality is also provided through the
    experiment table, by single/doubling clicking on experiment elements.
    """
    VERSION = (0, 5)
    def __init__(self):
        self.check_database()
        QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
        self.app = QtWidgets.QApplication([])
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
        self.experiment_table = None

        self.add_menu()
        self.decorate_window()
        self.window.setLayout(self.base_layout)

    def check_database(self):
        """Determine the database version and update it if necessary.

        Arguments: N/A
        Returns: None
        """
        version_check = database.checkTable("version")
        if not version_check:
            database.createVersionTable()
            if database.checkTable("experiments"):
                database.insertVersion(0, 0)
            else:
                database.insertVersion(*self.VERSION)

        known_version = database.getVersion()
        if known_version != self.VERSION:
            database.run_database_updates(known_version, self.VERSION)

    def add_menu(self):
        """Create the File menu for PombeTrack in the base window.

        Arguments: N/A
        Returns: None
        """
        menubar = QtWidgets.QMenuBar(self.window)
        file_menu = menubar.addMenu("&File")

        refresh_action = QtWidgets.QAction("&Refresh", menubar)
        refresh_action.triggered.connect(self._refresh_layout)
        file_menu.addAction(refresh_action)

        quit_action = QtWidgets.QAction("&Quit", menubar)
        quit_action.triggered.connect(self.window.close)
        file_menu.addAction(quit_action)

        self.base_layout.setMenuBar(menubar)

    def initiate(self):
        """Display the window and start the application.

        Arguments: N/A
        Returns: None
        """
        self.window.show()
        self.app.exec_()

    @staticmethod
    def get_date_item(exp):
        """Create a widget displaying the date of an experiment.

        Arguments:
            exp (database.ExperimentRow):
                Data for the experiment in question.

        Returns:
            date_item (QtWidgets.QTableWidgetItem):
                Widget displaying the date of the experiment.
        """
        date_item = QtWidgets.QTableWidgetItem()
        date_item.setText(exp.date)
        date_item.setTextAlignment(QtCore.Qt.AlignCenter)
        date_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        date_item.item_type = "experiment_date"
        date_item.experiment_id = exp.experiment_id
        return date_item

    @staticmethod
    def get_medium_item(exp):
        """Create a widget displaying the medium of an experiment.

        Arguments:
            exp (database.ExperimentRow):
                Data for the experiment in question.

        Returns:
            medium_item (QtWidgets.QTableWidgetItem):
                Widget displaying the medium of the experiment.
        """
        medium_item = QtWidgets.QTableWidgetItem()
        medium_item.setText(exp.medium)
        medium_item.setTextAlignment(QtCore.Qt.AlignCenter)
        medium_item.item_type = "experiment_medium"
        medium_item.experiment_id = exp.experiment_id
        return medium_item

    @staticmethod
    def get_strain_item(exp):
        """Create a widget displaying the strain of an experiment.

        Arguments:
            exp (database.ExperimentRow):
                Data for the experiment in question.

        Returns:
            strain_item (QtWidgets.QTableWidgetItem):
                Widget displaying the strain of the experiment.
        """
        strain_item = QtWidgets.QTableWidgetItem()
        strain_item.setText(exp.strain)
        strain_item.setTextAlignment(QtCore.Qt.AlignCenter)
        strain_item.item_type = "experiment_strain"
        strain_item.experiment_id = exp.experiment_id
        return strain_item

    @staticmethod
    def get_outlined_item(exp):
        """Create a widget displaying the outlined status of an experiment.

        Arguments:
            exp (database.ExperimentRow):
                Data for the experiment in question.

        Returns:
            outlined_item (QtWidgets.QTableWidgetItem):
                Widget displaying the outlined status of the experiment.

        The widget is a simple cell with the text YES, or NO, depending on
        whether cells have been outlined for the experiment.
        Its background colour is green or red respectively.
        """
        outlined_item = QtWidgets.QTableWidgetItem()
        outlined_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        if exp.outlined:
            outlined_item.setData(
                QtCore.Qt.BackgroundRole,
                QtGui.QBrush(QtGui.QColor("green"))
            )
            num_outlines = len(database.getOutlinesByExperimentId(exp.experiment_id))
            outlined_item.setText(str(num_outlines))

        else:
            outlined_item.setData(
                QtCore.Qt.BackgroundRole,
                QtGui.QBrush(QtGui.QColor("red"))
            )
            outlined_item.setText("NO")

        outlined_item.setTextAlignment(QtCore.Qt.AlignCenter)
        outlined_item.experiment_id = exp.experiment_id
        outlined_item.item_type = "experiment_outlined"
        return outlined_item

    @staticmethod
    def get_verified_item(exp):
        """Create a widget displaying the verification status of an experiment.

        Arguments:
            exp (database.ExperimentRow):
                Data for the experiment in question.

        Returns:
            verified_item (QtWidgets.QTableWidgetItem):
                Widget displaying the verification status of the experiment.

        The widget is a simple cell with the text YES, or NO, depending on
        whether cells have been verified.
        Its background colour is green or red respectively.
        The verification status of an experiment is reset if any outlines are
        added, removed, or modified.
        """
        verified_item = QtWidgets.QTableWidgetItem()
        verified_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        if exp.verified:
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
        verified_item.experiment_id = exp.experiment_id
        verified_item.item_type = "experiment_verified"
        return verified_item

    @staticmethod
    def get_analysed_item(exp):
        """Create a widget displaying the analysis status of an experiment.

        Arguments:
            exp (database.ExperimentRow):
                Data for the experiment in question.

        Returns:
            analysed_item (QtWidgets.QTableWidgetItem):
                Widget displaying the analysis status of the experiment.

        The widget is a simple cell with the text YES, or NO, depending on
        whether cells have been analysed.
        Its background colour is green or red respectively.
        """
        analysed_item = QtWidgets.QTableWidgetItem()
        analysed_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        if exp.analysed:
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
        analysed_item.experiment_id = exp.experiment_id
        analysed_item.item_type = "experiment_analysed"
        return analysed_item

    @staticmethod
    def get_cells_counted_item(exp):
        """Create a widget displaying the number of cells in an experiment.

        Arguments:
            exp (database.ExperimentRow):
                Data for the experiment in question.

        Returns:
            cells_counted_item (QtWidgets.QTableWidgetItem):
                Widget displaying the number of cells in the experiment.

        The number of cells is defined as verified cells in which both birth
        and division are observed.
        In brackets in this widget, the total number of cells is also shown,
        including those in which birth or division are not observed.
        """
        cells = database.getCellsByExperimentId(exp.experiment_id)
        cells_wanted = database.getCellsByExperimentId(
            exp.experiment_id,
            birth_observed=True,
            division_observed=True,
        )
        cell_count = "{0} ({1})".format(len(cells_wanted), len(cells))

        cells_counted_item = QtWidgets.QTableWidgetItem()
        cells_counted_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
        if exp.verified:
            cells_counted_item.setText(cell_count)
        else:
            cells_counted_item.setText("0")

        cells_counted_item.setTextAlignment(QtCore.Qt.AlignCenter)
        cells_counted_item.experiment_id = exp.experiment_id
        cells_counted_item.item_type = "experiment_cells_counted"
        return cells_counted_item

    def create_experiment_table(self):
        """Create the experiments table and populate it.

        Arguments: N/A

        Returns:
            experiment_table (QtWidgets.QTableWidget):
                The experiments table.

        The table has the following columns:
            Date, Medium, Strain, Outlined, Verified, Analysed, and Cells
            detected

        The first 3 of these columns (Date, Medium, Strain) are defined when
        adding a new experiment (or when edited).

        The next 3 columns (Outlined, Verified, Analysed) are defined by user
        actions.

        The final column depends on the verification process to define cells.

        When a table row is clicked, the view/edit and delete buttons are
        re-associated to that experiment (table_click_event function).

        When a table row is double clicked on a non-editable cell, the view
        interface for that experiment is opened (view_experiment function).
        Editable cells are Medium and Strain.
        If these are double clicked, the user can edit the contents of the cell
        and submit it (table_change_event function).
        """
        existing_experiments = sorted(
            database.getExperiments(),
            key=lambda x: (x.date, x.strain),
        )

        experiment_table = QtWidgets.QTableWidget()
        experiment_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        experiment_table.setRowCount(len(existing_experiments))
        header_labels = [
            "Date", "Medium", "Strain", "Outlined", "Verified", "Analysed", "Cells detected",
        ]
        item_methods = [
            self.get_date_item,
            self.get_medium_item,
            self.get_strain_item,
            self.get_outlined_item,
            self.get_verified_item,
            self.get_analysed_item,
            self.get_cells_counted_item,
        ]
        experiment_table.setColumnCount(len(header_labels))
        experiment_table.setHorizontalHeaderLabels(header_labels)
        for row_idx, exp in enumerate(existing_experiments):
            for col_idx, method in enumerate(item_methods):
                experiment_table.setItem(row_idx, col_idx, method(exp))

        experiment_table.itemClicked.connect(self.table_click_event)
        experiment_table.itemDoubleClicked.connect(self.view_experiment)
        experiment_table.itemChanged.connect(self.table_change_event)
        experiment_table.setCurrentCell(0, 0, QtCore.QItemSelectionModel.Clear)

        return experiment_table

    def decorate_window(self):
        """Create the base layout for PombeTrack.

        Arguments: N/A
        Returns: None

        Calls the create_experiment_table method, then adds buttons for adding,
        viewing, and deleting experiments.
        """
        header_label = QtWidgets.QLabel("Experiments:")
        header_font = QtGui.QFont("Fira Sans", 13)
        header_font.setBold(True)
        header_label.setFont(header_font)
        self.base_layout.addWidget(header_label)

        self.experiment_table = self.create_experiment_table()
        self.base_layout.addWidget(self.experiment_table)

        self.btn_row = QtWidgets.QHBoxLayout()
        btn = QtWidgets.QPushButton("Add new experiment")
        btn.clicked.connect(self.add_experiment)
        self.btn_row.addWidget(btn)

        self.edit_btn = QtWidgets.QPushButton("View experiment")
        self.edit_btn.clicked.connect(self.view_experiment)
        self.btn_row.addWidget(self.edit_btn)

        self.delete_btn = QtWidgets.QPushButton("Delete experiment")
        self.delete_btn.clicked.connect(self.delete_experiment)
        self.btn_row.addWidget(self.delete_btn)

        self.btn_row.setAlignment(QtCore.Qt.AlignLeft)
        self.base_layout.addLayout(self.btn_row)

    def table_change_event(self, item):
        """Respond to an item edit in the experiment table.

        Arguments:
            item (QtWidgets.QTableWidgetItem):
                The item that has been edited.

        Returns: None

        Should be triggered by Qt5 itself by signal binding.

        Updates the database accordingly and refresh the interface.
        """
        if not hasattr(item, "item_type"):
            return

        permitted_types = ["experiment_medium", "experiment_strain"]
        if item.item_type not in permitted_types:
            return

        data = {}
        item_name = item.item_type.split("experiment_")[1]
        data[item_name] = item.data(QtCore.Qt.DisplayRole)
        database.updateExperimentById(
            item.experiment_id,
            **data
        )
        self._refresh_layout()

    def table_click_event(self, item):
        """Respond to a single click of an experiment table item.

        Arguments:
            item (QtWidgets.QTableWidgetItem):
                The item that has been clicked.

        Returns: None

        Should be triggered by a Qt5 signal.

        Simply assigns the experiment_id of the clicked row to the view/edit
        and delete buttons.
        """
        self.edit_btn.experiment_id = item.experiment_id
        self.delete_btn.experiment_id = item.experiment_id

    def add_experiment(self):
        """Create dialog for creating new experiments.

        Arguments: N/A
        Returns: None

        See the experiment.Experiment.create_new_experiment method for more
        details.

        Will update the interface when the creation dialog is closed.
        """
        exp = experiment.Experiment()
        dialog = QtWidgets.QDialog(self.window)
        dialog.setModal(True)
        dialog.finished.connect(self._refresh_layout)
        exp.create_new_experiment(window=dialog)

    def view_experiment(self, item=None):
        """Create dialog for viewing and interacting with an experiment.

        Arguments:
            item (QtWidgets.QTableWidgetItem):
                optional argument

        Returns: None

        Can be called by a Qt5 event (in which case an item is needed), or
        explicitly by pressing the view button.

        See the experiment.ExperimentView.show_experiment method for more
        details.

        Will update the interface when the dialog is closed.
        """
        forbidden_types = ["experiment_medium", "experiment_strain"]
        if item and item.item_type in forbidden_types:
            return

        if (not hasattr(self.experiment_table, "selectedItems") or
                not self.experiment_table.selectedItems()):
            return

        exp_view = experiment.ExperimentView(self.edit_btn.experiment_id)
        dialog = QtWidgets.QDialog(self.window)
        dialog.setModal(True)
        dialog.finished.connect(self._refresh_layout)
        exp_view.show_experiment(window=dialog)

    def delete_experiment(self):
        """Delete an experiment.

        Arguments: N/A
        Returns: None

        Uses the experiment.ExperimentView.delete_experiment method.

        Refreshes the interface when complete.
        """
        if (not hasattr(self.experiment_table, "selectedItems") or
                not self.experiment_table.selectedItems()):
            return

        exp_view = experiment.ExperimentView(self.delete_btn.experiment_id)
        exp_view.window = self.window
        exp_view.delete_experiment(close_window=False)
        self._refresh_layout()

    def _clear_layout(self, layout):
        """Remove all elements from a layout.

        Arguments:
            layout (QtWidgets.QHBoxLayout or
                    QtWidgets.QVBoxLayout):
                The layout to be cleared (may also accept other layouts, simply
                needs the methods `count` and `takeAt`, though recursion may be
                affected if children are other types).

        Returns: None

        Accesses each widget in a layout and deletes it.
        Will recurse if the widget is itself a QtWidgets.QHBoxLayout or
        QtWidgets.QVBoxLayout.
        """

        for i in reversed(range(layout.count())):
            element = layout.takeAt(i)
            if isinstance(element, (QtWidgets.QHBoxLayout,
                                    QtWidgets.QVBoxLayout)):
                self._clear_layout(element)
            else:
                element.widget().deleteLater()

    def _refresh_layout(self):
        """Recreate the base layout and all its contents.

        Arguments: N/A
        Returns: None
        """
        # refresh table
        self._clear_layout(self.base_layout)
        self.decorate_window()
