import traceback
from PySide2.QtWidgets import QLabel, QComboBox
from PySide2.QtWidgets import (
    QWidget, QGridLayout, QScrollArea
)
from PySide2.QtCore import Qt, Slot, QAbstractItemModel, QModelIndex, QTimer

import sys
import os
import time
sys.path.append('C:/Program Files/STKO')  # Ensure current directory is in path for imports
sys.path.append(os.path.dirname(__file__) + '/../../')  # Ensure current directory is in path for imports
from NLTH.odb_utils.mpco_evaluator import NLTHEvaluator
from NLTH.odb_utils.mpco_document_utils import get_all_selection_sets
from NLTH.gui.mpl_line_plot_widget import MplLinePlotContainer
import math
from typing import Dict, List, Tuple

class _IsolatorSelectionSetComboModel(QAbstractItemModel):
    """Model for selection sets combo box."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._selection_sets = {}  # {selection_set_id: selection_set_name}
        self._db = None
    
    def load_selection_sets(self, db=None):
        """Load selection sets from odb_utils."""
        self.beginResetModel()
        self._db = db
        if self._db:
            self._selection_sets = get_all_selection_sets(db)
        else:
            self._selection_sets = {}
        self.endResetModel()
    
    def index(self, row, column, parent=QModelIndex()):
        if not self.hasIndex(row, column, parent):
            return QModelIndex()
        
        if not parent.isValid() and column == 0:  # Single column
            if row < len(self._selection_sets):
                sset_id = list(self._selection_sets.keys())[row]
                return self.createIndex(row, column, sset_id)
        
        return QModelIndex()
    
    def parent(self, index):
        return QModelIndex()  # Flat structure, no parent-child relationships
    
    def rowCount(self, parent=QModelIndex()):
        if not parent.isValid():
            return len(self._selection_sets)
        return 0
    
    def columnCount(self, parent=QModelIndex()):
        return 1  # Single column for combo box
    
    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        
        sset_id = index.internalId()
        sset_name = self._selection_sets.get(sset_id, "")
        
        if role == Qt.DisplayRole:
            return sset_name
        elif role == Qt.UserRole:  # Store the ID for retrieval
            return sset_id
        
        return None
    
    def get_selection_set_id(self, index):
        """Get selection set ID for given index."""
        if index >= 0 and index < len(self._selection_sets):
            return list(self._selection_sets.keys())[index]
        return None
    
    def get_selection_set_name(self, index):
        """Get selection set name for given index."""
        if index >= 0 and index < len(self._selection_sets):
            sset_id = list(self._selection_sets.keys())[index]
            return self._selection_sets[sset_id]
        return None

class _IsolatorPlotMode:
    TIME_DISP = 0
    TIME_FORCE = 1
    DISP_FORCE = 2

def _plot_labels_for_mode(mode : _IsolatorPlotMode, component: int):
    """Get axis labels for given plot mode and component."""
    component_labels = ('X', 'Y', 'Z')
    comp_label = component_labels[component]
    if mode == _IsolatorPlotMode.TIME_DISP:
        xlabel = "Time"
        ylabel = f"Displacement {comp_label}"
    elif mode == _IsolatorPlotMode.TIME_FORCE:
        xlabel = "Time"
        ylabel = f"Reaction {comp_label}"
    else:  # DISP_FORCE
        xlabel = f"Displacement {comp_label}"
        ylabel = f"Reaction {comp_label}"
    return xlabel, ylabel

class _IsolatorPlotWidget(QWidget):
    """
    A widget that manages multiple MplLinePlotContainer widgets in an adaptive grid layout.
    Creates one plot container for each database ID in the NLTHEvaluator.
    """
    
    def __init__(self, 
                 parent = None, 
                 plot_mode : _IsolatorPlotMode = _IsolatorPlotMode.DISP_FORCE, 
                 component = 0):
        """
        Initialize the NLTH plot widget.
        
        Args:
            parent: Parent widget
            plot_mode: Plot mode (time vs displacement, time vs force, displacement vs force)
            component: 0=X, 1=Y, 2=Z component
        """
        super().__init__(parent)
        
        self.plot_mode = plot_mode
        self.component = component
        
        # Main layout
        main_layout = QGridLayout()
        self.setLayout(main_layout)
        
        # Create scroll area for plots
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        main_layout.addWidget(self.scroll_area, 0, 0)
        
        # Container widget for plots
        self.plots_container = QWidget()
        self.plots_layout = QGridLayout()
        self.plots_container.setLayout(self.plots_layout)
        self.scroll_area.setWidget(self.plots_container)
        
        # Storage for plot widgets
        self.plot_widgets : Dict[int, MplLinePlotContainer] = {}  # {db_id: MplLinePlotContainer}
        self.node_line_map : Dict[int, Dict[int, int]] = {}  # {db_id: {node_id: line_index}}
        self.evaluator : NLTHEvaluator = None
    
    def _calculate_grid_layout(self, num_plots):
        """
        Calculate optimal grid layout for the given number of plots.
        
        Args:
            num_plots: Number of plots to arrange
            
        Returns:
            tuple: (rows, cols)
        """
        if num_plots == 0:
            return (0, 0)
        elif num_plots == 1:
            return (1, 1)
        elif num_plots == 2:
            return (1, 2)
        elif num_plots <= 4:
            return (2, 2)
        elif num_plots <= 6:
            return (2, 3)
        elif num_plots <= 9:
            return (3, 3)
        else:
            # For more than 9 plots, try to keep roughly square
            cols = math.ceil(math.sqrt(num_plots))
            rows = math.ceil(num_plots / cols)
            return (rows, cols)
    
    def _clear_plots(self):
        """
        Clear all existing plot widgets.
        """
        # Remove all widgets from layout
        while self.plots_layout.count():
            child = self.plots_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        
        # Clear storage
        self.plot_widgets.clear()
        self.node_line_map.clear()

    def _create_plot_widgets(self):
        """
        Create plot widgets for each database ID.
        
        Args:
            evaluator: NLTHEvaluator object
        """
        
        # quick return
        if not self.evaluator:
            return

        # Get database IDs
        db_ids = self.evaluator.db_ids
        if not db_ids:
            return
        
        # Calculate grid layout
        rows, cols = self._calculate_grid_layout(len(db_ids))
        
        # Component label for axes
        xlabel, ylabel = _plot_labels_for_mode(self.plot_mode, self.component)
        
        # Create plot widgets
        for i, db_id in enumerate(db_ids):
            row = i // cols
            col = i % cols
            
            # Create plot widget
            plot_widget = MplLinePlotContainer(
                parent=self.plots_container,
                width=6, height=4, dpi=80,
                title=f"Database {db_id}",
                labelx=xlabel,
                labely=ylabel,
                show_toolbar=True
            )
            
            # Add to layout
            self.plots_layout.addWidget(plot_widget, row, col)
            
            # Store reference
            self.plot_widgets[db_id] = plot_widget
    
    def _populate_plot(self, db_id):
        """
        Populate a single plot based on current evaluator.
        
        Args:
            db_id: Database ID to populate
        """

        # get plot_widget for this db_id
        plot_widget = self.plot_widgets.get(db_id, None)
        if plot_widget is None:
            raise ValueError(f"No plot widget found for database ID {db_id}")
        
        # Get data source
        TIME = self.evaluator.db_times.get(db_id, None)
        DISP_dict = self.evaluator.isolator_displacements.get(db_id, None)
        REAC_dict = self.evaluator.isolator_reactions.get(db_id, None)
        if TIME is None or DISP_dict is None or REAC_dict is None:
            raise ValueError(f"Missing data for database ID {db_id}")
        
        # Initialize node line mapping for this database if not exists
        db_line_map = {}
        self.node_line_map[db_id] = db_line_map
        
        # one plot for each isolator (top_node_id, bot_node_id)
        for lindex_index, (top_node_id, bot_node_id) in enumerate(zip(DISP_dict.keys(), REAC_dict.keys())):

            # get data arrays
            if self.plot_mode == _IsolatorPlotMode.TIME_DISP:
                XSOURCE = TIME[:]
                YSOURCE = DISP_dict[top_node_id][:, self.component]
            elif self.plot_mode == _IsolatorPlotMode.TIME_FORCE:
                XSOURCE = TIME[:]
                YSOURCE = REAC_dict[bot_node_id][:, self.component]
            else:  # DISP_FORCE
                XSOURCE = DISP_dict[top_node_id][:, self.component]
                YSOURCE = REAC_dict[bot_node_id][:, self.component]

            # map top_node_id to line index
            db_line_map[top_node_id] = lindex_index

            # add new line to plot
            plot_widget.add_line(XSOURCE, YSOURCE, ls = '-')
        
        # Update plot using the optimized async method
        plot_widget.update_plot()   

    def _populate_plots(self):
        """
        Populate plots based on current evaluator.
        """
        if not self.evaluator:
            return
        for db_id in self.evaluator.db_ids:
            self._populate_plot(db_id)

    def _update_plots(self):
        
        try:

            # clear existing plots and recreate
            self._clear_plots()
            self._create_plot_widgets()

            # populate all plots
            self._populate_plots()

        except Exception as ex:
            print(f"Exception in _update_plots:\n{ex}")
            print(traceback.format_exc())

    def set_evaluator(self, evaluator : NLTHEvaluator):
        """Set the NLTHEvaluator for the plot widget."""
        if evaluator is not self.evaluator:
            self.evaluator = evaluator
            self._update_plots()
    
    def set_plot_settings(self, 
                        plot_mode : _IsolatorPlotMode = None, 
                        component: int = None):
        """Set plot settings."""
        do_update = False
        if plot_mode is not None:
            self.plot_mode = plot_mode
            do_update = True
        if component is not None:
            self.component = component
            do_update = True
        # Recreate plots with new settings
        if do_update:
            self._update_plots()


class IsolatorWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QGridLayout()
        self.setLayout(layout)

        label = QLabel("Isolators selection set")
        self.isolators_combo = QComboBox()
        
        # Create and set model
        self.model = _IsolatorSelectionSetComboModel(self)
        self.isolators_combo.setModel(self.model)

        layout.addWidget(label, 0, 0)
        layout.addWidget(self.isolators_combo, 0, 1)
        layout.setColumnStretch(2, 1)

        self.plot_widget = _IsolatorPlotWidget(parent=self)
        layout.addWidget(self.plot_widget, 1, 0, 1, 3)

        self.refresh()
    
    def refresh(self, db=None):
        """Reload selection sets from odb_utils."""
        self.model.load_selection_sets(db)
    
    def get_selected_set_id(self):
        """Get the currently selected selection set ID."""
        current_index = self.isolators_combo.currentIndex()
        return self.model.get_selection_set_id(current_index)
    
    def get_selected_set_name(self):
        """Get the currently selected selection set name."""
        current_index = self.isolators_combo.currentIndex()
        return self.model.get_selection_set_name(current_index)
    
    def set_evaluator(self, evaluator : NLTHEvaluator):
        """Set the NLTHEvaluator for the plot widget."""
        self.plot_widget.set_evaluator(evaluator)



