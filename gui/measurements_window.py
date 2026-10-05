from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QTabWidget, QWidget, QTableWidget,QFileDialog,QMessageBox,
    QTableWidgetItem, QHeaderView, QAbstractItemView, QPushButton, QHBoxLayout,QCheckBox,QWidget,QHBoxLayout
)
from PyQt5.QtCore import Qt
import csv

class MeasurementsWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Measurements")
        self.setMinimumSize(420, 420)
        self.setWindowFlags(Qt.Window | Qt.WindowStaysOnTopHint)

        self.setStyleSheet("""
            QDialog { background-color: #f8f9fa; }
            QTableWidget {
                background-color: #ffffff; color: #1f2937;
                border: 1px solid #e5e7eb; border-radius: 8px;
                font-size: 12px; gridline-color: transparent; outline: none;
            }
            QTableWidget::item {
                padding: 8px 12px; color: #1f2937;
                border-bottom: 1px solid #f3f4f6;
            }
            QTableWidget::item:selected { background-color: #f0f4ff; color: #1f2937; }
            QHeaderView::section {
                background-color: #f9fafb; color: #6b7280; font-weight: 600;
                font-size: 12px; padding: 7px 12px; border: none;
                border-bottom: 1px solid #e5e7eb;
            }
            QTabWidget::pane {
                border: 1px solid #e5e7eb; border-radius: 8px;
                background-color: #ffffff; top: -1px;
            }
            QTabBar::tab {
                color: #1f2937; padding: 12px 17px;
                background-color: #f3f4f6;
                border: 1px solid #e5e7eb; border-bottom: none;
                border-top-left-radius: 8px; border-top-right-radius: 8px;
                margin-right: 2px;
            }
            QTabBar::tab:selected { background-color: #ffffff; }
            QPushButton {
                background-color: #ffffff; color: #1f2937; font-weight: 600;
                border: 1px solid #e5e7eb; border-radius: 8px; padding: 6px 12px;
            }
            QPushButton:hover { background-color: #f3f4f6; }
        """)

        layout = QVBoxLayout(self)
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        self.dist_table = self._make_table(["Atoms", "Distance (Å)"])
        self.angle_table = self._make_table(["Atoms", "Angle (°)"])
        self.dihedral_table = self._make_table(["Atoms", "Dihedral (°)"])

        self.tabs.addTab(self._wrap(self.dist_table), "Distances")
        self.tabs.addTab(self._wrap(self.angle_table), "Angles")
        self.tabs.addTab(self._wrap(self.dihedral_table), "Dihedral")

        
        btn_row = QHBoxLayout()
        export_btn = QPushButton("Export CSV")
        export_btn.clicked.connect(self._export_csv)
        clear_btn = QPushButton("Clear current tab")
        clear_btn.clicked.connect(self._clear_current)
        btn_row.addWidget(export_btn)
        btn_row.addStretch()
        btn_row.addWidget(clear_btn)
        layout.addLayout(btn_row)

    def _make_table(self, headers):
        t = QTableWidget()
        t.setColumnCount(len(headers)+1)
        t.setHorizontalHeaderLabels(['Show'] + headers)
        t.setEditTriggers(QAbstractItemView.NoEditTriggers)
        t.verticalHeader().setVisible(False)
        t.setShowGrid(False)
        t.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        t.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        t.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        return t

    def _wrap(self, table):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.addWidget(table)
        return w

    def add_distance(self, atoms_label, value, mid):
        self._add_row(self.dist_table, atoms_label, f"{value:.3f}", mid)
        self.tabs.setCurrentIndex(0)

    def add_angle(self, atoms_label, value, mid):
        self._add_row(self.angle_table, atoms_label, f"{value:.2f}", mid)
        self.tabs.setCurrentIndex(1)

    def add_dihedral(self, atoms_label, value, mid):
        self._add_row(self.dihedral_table, atoms_label, f"{value:.2f}", mid)
        self.tabs.setCurrentIndex(2)

    def _add_row(self, table, label, value,mid):
        row = table.rowCount()
        table.insertRow(row)
        table.setItem(row, 0, QTableWidgetItem(label))
        table.setItem(row, 1, QTableWidgetItem(value))
        chk = QCheckBox()
        chk.stateChanged.connect(
            lambda state, m=mid: self.parent().set_measure_visible(m, state == Qt.Checked))
        holder = QWidget()
        lay = QHBoxLayout(holder)
        lay.addWidget(chk)
        lay.setAlignment(Qt.AlignCenter)
        lay.setContentsMargins(0, 0, 0, 0)
        table.setCellWidget(row, 0, holder)

        table.setItem(row, 1, QTableWidgetItem(label))
        table.setItem(row, 2, QTableWidgetItem(value))
        item_atoms = QTableWidgetItem(label)
        item_atoms.setData(Qt.UserRole, mid)
        table.setItem(row, 1, item_atoms)
        table.setItem(row, 2, QTableWidgetItem(value))

    def _clear_current(self):
        idx = self.tabs.currentIndex()
        table = [self.dist_table, self.angle_table, self.dihedral_table][idx]
        for row in range(table.rowCount()):
            item =table.item(row, 1)
            if item is not None:
                mid = item.data(Qt.UserRole)
                if mid is not None:
                    self.parent().remove_measure(mid)
        table.setRowCount(0)

    def clear_all(self):
        for t in (self.dist_table, self.angle_table, self.dihedral_table):
            t.setRowCount(0)

    def _export_csv(self):
        tables = [
            ("Distances (Å)", self.dist_table),
            ("Angles (°)", self.angle_table),
            ("Dihedrals (°)", self.dihedral_table),
        ]
        if all(t.rowCount() == 0 for _, t in tables):
            QMessageBox.information(self, "Export", "No measurements to export yet.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export measurements", "measurements.csv", "CSV Files (*.csv)")
        if not file_path:
            return

        try:
            with open(file_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f, delimiter=';')
                for title, table in tables:
                    if table.rowCount() == 0:
                        continue
                    writer.writerow([title])
                    writer.writerow(["Atoms", "Value"])
                    for row in range(table.rowCount()):
                        atoms = table.item(row, 1).text()
                        value = table.item(row, 2).text()
                        writer.writerow([atoms, value])
                    writer.writerow([])   
        except Exception as e:
            QMessageBox.warning(self, "Export", f"Could not export:\n{e}")