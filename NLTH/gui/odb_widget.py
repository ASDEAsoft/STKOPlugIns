from PySide2.QtWidgets import QWidget, QTreeView, QVBoxLayout, QHeaderView
from PySide2.QtCore import Qt, QAbstractItemModel, QModelIndex, Signal
from PySide2.QtGui import QStandardItemModel, QStandardItem
import os
from NLTH.odb_utils.mpco_document_utils import get_all_databases

class DatabaseTreeModel(QAbstractItemModel):
    """Model for database tree with checkable items and two columns."""
    
    # Signal emitted when the list of checked items changes
    checkedItemsChanged = Signal(list)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._databases = {}  # {database_id: database_name}
        self._checked_ids = set()
        self._headers = ["Filename", "Directory"]
    
    def load_databases(self):
        """Load databases from odb_utils."""
        self.beginResetModel()
        self._databases = get_all_databases()
        self._checked_ids.clear()
        self.endResetModel()
    
    def index(self, row, column, parent=QModelIndex()):
        if not self.hasIndex(row, column, parent):
            return QModelIndex()
        
        if not parent.isValid():  # Root items
            if row < len(self._databases) and 0 <= column < len(self._headers):
                db_id = list(self._databases.keys())[row]
                return self.createIndex(row, column, db_id)
        
        return QModelIndex()
    
    def parent(self, index):
        return QModelIndex()  # Flat structure, no parent-child relationships
    
    def rowCount(self, parent=QModelIndex()):
        if not parent.isValid():
            return len(self._databases)
        return 0
    
    def columnCount(self, parent=QModelIndex()):
        return len(self._headers)
    
    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if 0 <= section < len(self._headers):
                return self._headers[section]
        return None
    
    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        
        db_id = index.internalId()
        filepath = self._databases.get(db_id, "")
        
        if role == Qt.DisplayRole:
            if index.column() == 0:  # Filename column
                return os.path.basename(filepath)
            elif index.column() == 1:  # Directory column
                return os.path.dirname(filepath)
        elif role == Qt.CheckStateRole and index.column() == 0:  # Only first column checkable
            return Qt.Checked if db_id in self._checked_ids else Qt.Unchecked
        
        return None
    
    def setData(self, index, value, role=Qt.EditRole):
        if not index.isValid():
            return False
        
        if role == Qt.CheckStateRole and index.column() == 0:
            db_id = index.internalId()
            if value == Qt.Checked:
                self._checked_ids.add(db_id)
            else:
                self._checked_ids.discard(db_id)
            self.dataChanged.emit(index, index, [Qt.CheckStateRole])
            # Emit signal with current list of checked IDs
            self.checkedItemsChanged.emit(list(self._checked_ids))
            return True
        
        return False
    
    def flags(self, index):
        if not index.isValid():
            return Qt.NoItemFlags
        
        flags = Qt.ItemIsEnabled | Qt.ItemIsSelectable
        if index.column() == 0:  # Only first column is checkable
            flags |= Qt.ItemIsUserCheckable
        return flags
    
    def get_checked_ids(self):
        """Return list of checked database IDs."""
        return list(self._checked_ids)


class DatabaseWidget(QWidget):
    """Widget displaying a tree view of checkable databases with filename and directory columns."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Create layout
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Create tree view
        self.tree_view = QTreeView()
        self.tree_view.setRootIsDecorated(False)  # No expand/collapse indicators
        self.tree_view.setAlternatingRowColors(True)
        
        # Create and set model
        self.model = DatabaseTreeModel(self)
        self.tree_view.setModel(self.model)
        
        # Configure column widths
        header = self.tree_view.header()
        header.setStretchLastSection(True)
        header.resizeSection(0, 200)  # Filename column width

        # Add to layout
        layout.addWidget(self.tree_view)
        
        # Load databases
        self.model.load_databases()
    
    @property
    def checkedItemsChanged(self):
        """Signal emitted when checked items change. Returns the model's signal."""
        return self.model.checkedItemsChanged
    
    def get_checked_database_ids(self):
        """Return list of checked database IDs."""
        return self.model.get_checked_ids()
    
    def refresh(self):
        """Reload databases from odb_utils."""
        self.model.load_databases()