from PyMpc import App, IO
from PySide2.QtCore import (
    QObject, QThread, Signal, QTimer
)
from PySide2.QtWidgets import (
    QDialog, QTabWidget, 
    QVBoxLayout, QHBoxLayout, 
    QPushButton, QProgressBar,
)
from NLTH.gui.isolator_widget import IsolatorWidget
from NLTH.gui.odb_widget import DatabaseWidget
from NLTH.odb_utils.mpco_evaluator import NLTHEvaluator, NLTHEvaluatorCallbackBase
import traceback

class NLTHEvaluatorCallback(NLTHEvaluatorCallbackBase, QObject):
    """Callback class that inherits from both NLTHEvaluatorCallbackBase and QObject to emit signals."""

    # Define signals

    # Emitted when progress changes (progress : float 0.0 to 1.0)
    progressUpdated = Signal(float)

    def __init__(self, parent=None):
        NLTHEvaluatorCallbackBase.__init__(self)
        QObject.__init__(self, parent)
        self._last_emitted_progress = 0.0
        self._progress_threshold = 0.01  # Only emit every 1% change
    
    def on_progress(self, progress: float):
        """Called to report progress."""
        # Only emit signal if progress changed significantly to avoid GUI flooding
        if abs(progress - self._last_emitted_progress) >= self._progress_threshold or progress >= 1.0:
            self._last_emitted_progress = progress
            self.progressUpdated.emit(progress)

class EvaluationWorker(QObject):
    """Worker object for running the evaluation in the background."""
    
    # Worker signals
    finished = Signal()
    error = Signal(str)
    
    def __init__(self, evaluator : NLTHEvaluator, callback: NLTHEvaluatorCallback):
        super().__init__()
        self.callback = callback
        self.evaluator = evaluator
    
    def run_evaluation(self):
        """Run the evaluation in the background thread."""
        try:
            self.evaluator.evaluate()
        except Exception as e:
            error_msg = f"Error during evaluation: {e}\n{traceback.format_exc()}"
            self.error.emit(error_msg)
        finally:
            self.finished.emit()

class NLTHMainWindow(QDialog):
    """Main window for NLTH post-processing tools."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("NLTH Post-Processing Tools")
        self.setModal(True)
        
        # Create tab widget
        self.tab_widget = QTabWidget()
        
        # Create odb tab
        self.odb_widget = DatabaseWidget()
        self.tab_widget.addTab(self.odb_widget, "Databases")
        
        # Create isolator tab
        self.isolator_widget = IsolatorWidget()
        self.tab_widget.addTab(self.isolator_widget, "Isolation System")
        
        # Connect signals
        self.odb_widget.checkedItemsChanged.connect(self.on_checked_items_changed)
        
        # Set up layout
        layout = QVBoxLayout()
        layout.addWidget(self.tab_widget)
        
        # Create horizontal layout for buttons
        button_layout = QHBoxLayout()
        
        # Create Run button
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.on_run_clicked)
        
        # Create Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)  # Initially hidden
        
        # Add button and progress bar to layout, aligned left
        button_layout.addWidget(self.run_button)
        button_layout.addWidget(self.progress_bar)
        button_layout.addStretch()  # Pushes widgets to the left
        
        # Add button layout to main layout
        layout.addLayout(button_layout)

        self.setLayout(layout)
        
        # Set dialog size
        self.resize(600, 400)
        
        # Initialize worker thread and callback
        self.worker_thread = None
        self.worker = None
        self.callback = None
        self.evaluator = None
    
    def on_checked_items_changed(self, checked_ids):
        """Slot called when checked database items change."""
        try:
            if checked_ids:
                # Use the first checked database ID
                first_db_id = checked_ids[0]
                first_db = App.postDocument().databases.get(first_db_id, None)
                self.isolator_widget.refresh(first_db)
            else:
                # No databases checked, refresh with None
                self.isolator_widget.refresh(None)
        except Exception as e:
            print(f"Error refreshing isolator widget: {e}")
    
    def on_run_clicked(self):
        """Slot called when Run button is clicked."""
        try:
            # Prevent multiple runs
            if self.worker_thread is not None and self.worker_thread.isRunning():
                print("Evaluation is already running...")
                return
            
            # get selected database IDs
            db_ids = self.odb_widget.get_checked_database_ids()
            if not db_ids:
                print("No databases selected for evaluation.")
                return
                
            # get isolator selection set ID
            isolator_set_id = self.isolator_widget.get_selected_set_id()

            # Create callback and connect signals
            self.callback = NLTHEvaluatorCallback(self)
            self.callback.progressUpdated.connect(self.on_progress_updated)

            # Create evaluator
            self.evaluator = NLTHEvaluator(db_ids, isolator_set_id, self.callback)

            # Create worker and thread
            self.worker_thread = QThread()
            self.worker = EvaluationWorker(self.evaluator, self.callback)
            
            # Move worker to thread
            self.worker.moveToThread(self.worker_thread)
            
            # Connect thread and worker signals with proper cleanup order
            self.worker_thread.started.connect(self.worker.run_evaluation)
            self.worker.finished.connect(self.worker_thread.quit)
            
            # Connect cleanup signals - these will be called after thread quits
            self.worker_thread.finished.connect(self.on_evaluation_finished)
            
            # Connect error handling - but only from worker, not callback
            self.worker.error.connect(self.on_evaluation_error)
            
            # Start the thread
            self.worker_thread.start()
            
            # Update UI
            self.run_button.setEnabled(False)
            self.progress_bar.setVisible(True)
            self.progress_bar.setValue(0)
            print("Evaluation started...")
            
        except Exception as e:
            print(f"Error starting evaluation: {e}")
            print(traceback.format_exc())
    
    def on_progress_updated(self, progress: float):
        """Slot called when progress is updated."""
        percentage = int(progress * 100)
        self.progress_bar.setValue(percentage)

    def on_evaluation_finished(self):
        """Slot called when evaluation is finished."""
        print("Evaluation completed successfully.")
        self.run_button.setEnabled(True)
        self.progress_bar.setVisible(False)
        
        # Clean up worker and thread
        if self.worker:
            self.worker.deleteLater()
            self.worker = None
        if self.worker_thread:
            self.worker_thread.deleteLater()
            self.worker_thread = None
        self.callback = None

        # TODO: set evalutator to sub-widgets if needed
    
    def on_evaluation_error(self, error_message: str):
        """Slot called when an error occurs during evaluation."""
        IO.write_cerr(error_message)


