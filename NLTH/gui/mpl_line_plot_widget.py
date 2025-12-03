import matplotlib
# Make sure that we are using QT5
matplotlib.use('Qt5Agg')
from matplotlib.lines import Line2D
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from PySide2.QtCore import Signal, Slot, QEventLoop, QTimer
from PySide2.QtWidgets import QWidget, QVBoxLayout, QSizePolicy
import sys
import time
from collections import deque
from typing import List, Tuple

class MplLinePlotWidget(FigureCanvas):
    """
    A simple matplotlib-based widget for line plotting using PySide2.
    This widget provides basic line plotting functionality with matplotlib integration.
    """
    # signals
    requestUpdateFigure = Signal()
    updateFigureDone = Signal()
    
    def __init__(self, parent=None, width=5, height=4, dpi=100,
                 title=None, labelx=None, labely=None):
        """
        Initialize the matplotlib canvas widget.
        
        Args:
            parent: Parent widget
            width: Figure width in inches
            height: Figure height in inches
            dpi: Dots per inch for the figure
            title: Plot title
            labelx: X-axis label
            labely: Y-axis label
        """
        self.figure = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, self.figure)
        self.control = parent
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        # Set up the axes
        self.axes = self.figure.add_subplot(111)
        self.axes.grid(True, ls=':', color=(0.6, 0.6, 0.6))
        if title:
            self.axes.set_title(title)
        if labelx:
            self.axes.set_xlabel(labelx)
        if labely:
            self.axes.set_ylabel(labely)
        
        self.setParent(parent)
        FigureCanvas.updateGeometry(self)
        
        # Data storage
        self._line_data_x: List[List[float]] = []
        self._line_data_y: List[List[float]] = []
        self._line_styles: List[str] = []
        self._line_objects: List[Line2D] = []
        
        # Performance optimizations
        self._update_pending = False
        self._last_update_time = 0
        self._min_update_interval = 0.05  # Minimum 50ms between updates
        self._update_queue = deque(maxlen=1)  # Only keep latest update
        
        # Timer for throttled updates
        self._update_timer = QTimer()
        self._update_timer.setSingleShot(True)
        self._update_timer.timeout.connect(self._perform_update)
        
        # Set up connections
        self.requestUpdateFigure.connect(self.onUpdateFigureRequested)
    
    def add(self, x, y, ls='-'):
        """
        Add a new line to the plot.
        
        Args:
            x: X data points (list or array)
            y: Y data points (list or array)
            ls: Line style (default: '-')
        """
        self._line_data_x.append(x)
        self._line_data_y.append(y)
        self._line_styles.append(ls)
        line = self.axes.plot(x, y, ls)[0]
        self._line_objects.append(line)
        # Use draw_idle() for better performance
        self.draw_idle()
    
    def clear(self):
        """Clear all lines from the plot."""
        # Cancel any pending updates
        if self._update_timer.isActive():
            self._update_timer.stop()
        self._update_pending = False
        
        self.axes.clear()
        self.axes.grid(True, ls=':', color=(0.6, 0.6, 0.6))
        self._line_data_x.clear()
        self._line_data_y.clear()
        self._line_styles.clear()
        self._line_objects.clear()
        # Use draw_idle() for better performance
        self.draw_idle()
    
    def set_labels(self, title=None, xlabel=None, ylabel=None):
        """
        Set plot labels.
        
        Args:
            title: Plot title
            xlabel: X-axis label
            ylabel: Y-axis label
        """
        if title:
            self.axes.set_title(title)
        if xlabel:
            self.axes.set_xlabel(xlabel)
        if ylabel:
            self.axes.set_ylabel(ylabel)
        # Use draw_idle() for better performance
        self.draw_idle()
    
    @Slot()
    def onUpdateFigureRequested(self):
        """Handle figure update requests with throttling."""
        current_time = time.time()
        
        # Throttle updates to prevent GIL contention
        if current_time - self._last_update_time < self._min_update_interval:
            if not self._update_pending:
                self._update_pending = True
                remaining_time = self._min_update_interval - (current_time - self._last_update_time)
                self._update_timer.start(int(remaining_time * 1000))
            return
        
        self._perform_update()
    
    def _perform_update(self):
        """Actually perform the figure update."""
        try:
            self._update_pending = False
            self._last_update_time = time.time()
            
            # Batch all line updates
            for i in range(len(self._line_objects)):
                self._line_objects[i].set_data(self._line_data_x[i], self._line_data_y[i])
            
            # Single relim/autoscale call
            self.axes.relim()
            self.axes.autoscale_view()
            
            # Use draw_idle() instead of draw() for better performance
            self.draw_idle()
        except Exception as ex:
            print(f"Exception in _perform_update:\n{ex}")
        finally:
            self.updateFigureDone.emit()
    
    def updateFigure(self):
        """Update the figure synchronously."""
        loop = QEventLoop()
        self.updateFigureDone.connect(loop.quit)
        self.requestUpdateFigure.emit()
        loop.exec_()
    
    def updateFigureAsync(self):
        """Update the figure asynchronously (non-blocking)."""
        self.requestUpdateFigure.emit()


class MplLinePlotContainer(QWidget):
    """
    A container widget that includes the matplotlib canvas and optional navigation toolbar.
    """
    
    def __init__(self, parent=None, width=5, height=4, dpi=100,
                 title=None, labelx=None, labely=None, show_toolbar=True):
        """
        Initialize the container widget.
        
        Args:
            parent: Parent widget
            width: Figure width in inches
            height: Figure height in inches
            dpi: Dots per inch for the figure
            title: Plot title
            labelx: X-axis label
            labely: Y-axis label
            show_toolbar: Whether to show the matplotlib navigation toolbar
        """
        super().__init__(parent)
        
        # Create layout
        layout = QVBoxLayout()
        
        # Create the matplotlib widget
        self.canvas = MplLinePlotWidget(self, width, height, dpi, title, labelx, labely)
        layout.addWidget(self.canvas)
        
        # Optionally add navigation toolbar
        if show_toolbar:
            self.toolbar = NavigationToolbar(self.canvas, self)
            layout.addWidget(self.toolbar)
        else:
            self.toolbar = None
        
        self.setLayout(layout)
    
    def add_line(self, x, y, ls='-'):
        """Add a line to the plot."""
        self.canvas.add(x, y, ls)

    def get_line_count(self) -> int:
        """Get the number of lines in the plot."""
        return len(self.canvas._line_objects)
    
    def get_line_data(self, index: int) -> Tuple[List[float], List[float]]:
        """Get the data of a specific line by index.
        
        Args:
            index: Index of the line
            
        Returns:
            Tuple of (x data, y data)
        """
        if 0 <= index < len(self.canvas._line_objects):
            return (self.canvas._line_data_x[index], self.canvas._line_data_y[index])
        else:
            raise IndexError("Line index out of range.")
    
    def clear_plot(self):
        """Clear the plot."""
        self.canvas.clear()
    
    def set_labels(self, title=None, xlabel=None, ylabel=None):
        """Set plot labels."""
        self.canvas.set_labels(title, xlabel, ylabel)
    
    def update_plot(self):
        """Update the plot using async method for better performance."""
        self.canvas.updateFigureAsync()


