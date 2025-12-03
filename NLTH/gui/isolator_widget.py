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
    
    def __init__(self, parent=None, plot_mode : _IsolatorPlotMode = _IsolatorPlotMode.DISP_FORCE, component=0):
        """
        Initialize the NLTH plot widget.
        
        Args:
            parent: Parent widget
            plot_mode: 'displacement' or 'reaction'
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
        self.initialized = False
        
        # Performance throttling
        self._last_update_time = 0
        self._min_update_interval = 0.1  # Minimum 100ms between updates
        self._update_pending = False
        
        # Timer for throttled updates
        self._update_timer = QTimer()
        self._update_timer.setSingleShot(True)
        self._update_timer.timeout.connect(self._perform_pending_update)
        self._pending_update_args = None
    
    def set_plot_mode(self, mode : _IsolatorPlotMode, component=0):
        """
        Set what to plot and update all existing plots.
        
        Args:
            mode: 'displacement' or 'reaction'
            component: 0=X, 1=Y, 2=Z component
        """
        self.plot_mode = mode
        self.component = component
        xlabel, ylabel = _plot_labels_for_mode(mode, component)        
        for plot_widget in self.plot_widgets.values():
            plot_widget.set_labels(xlabel=xlabel, ylabel=ylabel)
    
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
    
    def _create_plot_widgets(self, evaluator : NLTHEvaluator):
        """
        Create plot widgets for each database ID.
        
        Args:
            evaluator: NLTHEvaluator object
        """
        # Clear existing widgets
        self._clear_plots()
        
        # Get database IDs
        db_ids = evaluator.db_ids
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
                show_toolbar=False  # Don't show toolbar to save space
            )
            
            # Add to layout
            self.plots_layout.addWidget(plot_widget, row, col)
            
            # Store reference
            self.plot_widgets[db_id] = plot_widget
    
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
    
    def update_plot(self, evaluator : NLTHEvaluator, db_id, index, count):
        """
        Update plots with data from NLTHEvaluator (with throttling for performance).
        
        Args:
            evaluator: NLTHEvaluator object containing the data
            db_id: Database ID to update
            index: Starting index for the update
            count: Number of points to update
        """
        current_time = time.time()
        
        # Store the update arguments
        self._pending_update_args = (evaluator, db_id, index, count)
        
        # Throttle updates to prevent GIL contention and UI lag
        if current_time - self._last_update_time < self._min_update_interval:
            if not self._update_pending:
                self._update_pending = True
                remaining_time = self._min_update_interval - (current_time - self._last_update_time)
                self._update_timer.start(int(remaining_time * 1000))
            return
        
        self._perform_update(evaluator, db_id, index, count)
    
    def _perform_pending_update(self):
        """Perform the pending update from timer."""
        if self._pending_update_args:
            evaluator, db_id, index, count = self._pending_update_args
            self._perform_update(evaluator, db_id, index, count)
        self._update_pending = False
    
    def _perform_update(self, evaluator : NLTHEvaluator, db_id, index, count):
        """
        Actually perform the plot update.
        
        Args:
            evaluator: NLTHEvaluator object containing the data
            db_id: Database ID to update
            index: Starting index for the update
            count: Number of points to update
        """
        try:
            # Initialize on first call
            if not self.initialized or self.evaluator != evaluator:
                self._create_plot_widgets(evaluator)
                self.evaluator = evaluator
                self.initialized = True

            # get plot_widget for this db_id
            if db_id not in self.plot_widgets:
                return
            plot_widget = self.plot_widgets[db_id]
            
            # Get data source
            TIME = evaluator.db_times.get(db_id, None)
            DISP_dict = evaluator.isolator_displacements.get(db_id, None)
            REAC_dict = evaluator.isolator_reactions.get(db_id, None)
            if TIME is None or DISP_dict is None or REAC_dict is None:
                return

            # Determine time range to update (from input)
            end_index = index + count
            
            # Initialize node line mapping for this database if not exists
            if db_id not in self.node_line_map:
                self.node_line_map[db_id] = {}
            db_line_map = self.node_line_map[db_id]
            
            # For first call (index=0), clear and initialize lines
            if index == 0:
                plot_widget.clear_plot()
                db_line_map.clear()
                # Create lines for all nodes
                for node_id in DISP_dict.keys(): # always map with top node IDs
                    if node_id not in db_line_map:
                        # Add empty line for this node
                        plot_widget.add_line([], [], ls='-')
                        line_index = plot_widget.get_line_count() - 1
                        db_line_map[node_id] = line_index
            
            # Update existing lines for each node
            for top_node_id, bot_node_id in zip(DISP_dict.keys(), REAC_dict.keys()):

                # get data arrays
                x_1d = True
                if self.plot_mode == _IsolatorPlotMode.TIME_DISP:
                    XSOURCE = TIME
                    YSOURCE = DISP_dict[top_node_id]
                elif self.plot_mode == _IsolatorPlotMode.TIME_FORCE:
                    XSOURCE = TIME
                    YSOURCE = REAC_dict[bot_node_id]
                else:  # DISP_FORCE
                    x_1d = False
                    XSOURCE = DISP_dict[top_node_id]
                    YSOURCE = REAC_dict[bot_node_id]

                # get line index
                line_index = db_line_map[top_node_id]

                # get data from plot so we can update it
                x_data, y_data = plot_widget.get_line_data(line_index)

                # now we have to copy SOURCE[index:end_index, component] to the plot line
                # compute the start index in case index is beyond current data length
                start_index = min(index, len(x_data))
                for row_index in range(start_index, end_index):
                    ix = XSOURCE[row_index] if x_1d else XSOURCE[row_index, self.component]
                    iy = YSOURCE[row_index, self.component]
                    if row_index < len(x_data):
                        x_data[row_index] = ix
                        y_data[row_index] = iy
                    else:
                        x_data.append(ix)
                        y_data.append(iy)
            
            # Update plot using the optimized async method
            plot_widget.update_plot()
            
            # Update last update time
            self._last_update_time = time.time()
            
        except Exception as ex:
            print(f"Exception in _IsolatorPlotWidget._perform_update: {ex}")
    
    def clear_all_plots(self):
        """
        Clear all plots and reset the widget.
        """
        for plot_widget in self.plot_widgets.values():
            plot_widget.clear_plot()
        
        self.node_line_map.clear()
        self.initialized = False
        self.evaluator = None
    
    def get_plot_widget(self, db_id):
        """
        Get the plot widget for a specific database ID.

        Args:
            db_id: Database ID
            
        Returns:
            MplLinePlotContainer or None
        """
        return self.plot_widgets.get(db_id, None)

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
    
    @Slot(int, int, int, NLTHEvaluator)
    def on_update_range_requested(self, db_id: int, index: int, count: int, evaluator: NLTHEvaluator):
        """Slot to handle update range requests."""
        self.plot_widget.update_plot(evaluator, db_id=db_id, index=index, count=count)
        print(f"      IsolatorWidget received update range request: DB ID={db_id}, Index={index}, Count={count}")



