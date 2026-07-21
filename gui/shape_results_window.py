from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QTableWidget,QPushButton,
    QTableWidgetItem, QHeaderView, QAbstractItemView,QColorDialog
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont


class ShapeResultsWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Shape Measures")
        self.setMinimumSize(400, 350)
        self.setWindowFlags(Qt.Window | Qt.WindowStaysOnTopHint)
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(10)
        self.setLayout(layout)

        self.header_label = QLabel("")
        self.header_label.setStyleSheet(
            "font-size: 13px; font-weight: bold; color: #1d4ed8;"
        )
        layout.addWidget(self.header_label)

        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Shape Reference", "CShM"])
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.setShowGrid(False)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setStyleSheet("""
            QTableWidget {
                background-color: #ffffff;
                alternate-background-color: #2e3e69;
                border: 1px solid #d1d5db;
                border-radius: 6px;
                font-size: 12px;
                gridline-color: transparent;
            }
            QTableWidget::item {
                padding: 6px 10px;
                border: none;
                color: #1f2937;
            }
            QTableWidget::item:!selected {
                background-color: #fcfdff ;
                color: #000000;
            }
            QHeaderView::section {
                background-color: #2e3e69;
                color: white;
                font-weight: bold;
                font-size: 12px;
                padding: 6px 10px;
                border: none;
            }
        """)
        hdr = self.table.horizontalHeader()
        hdr.setSectionResizeMode(0, QHeaderView.Stretch)
        hdr.setSectionResizeMode(1, QHeaderView.Stretch)
        layout.addWidget(self.table)

        # En _build_ui, debajo de la tabla:
        self.show_poly_button = QPushButton("Show Polyhedra")
        self.show_poly_button.setVisible(False)
        self.show_poly_button.setStyleSheet("""
            QPushButton {
            background-color: #74034e;
            color: white;
            font-weight: bold;
            border-radius: 4px;
            border: 2px solid black;
            }
            QPushButton:hover { background-color: #74034e }
            """)
        layout.addWidget(self.show_poly_button)

        #boton color poliedro
        self.poly_color_button = QPushButton("Set Polyhedron Color")
        self.poly_color_button.setVisible(False)
        self.poly_color_button.setStyleSheet("""
            QPushButton {
            background-color: #74034e;
            color: white;
            font-weight: bold;
            border-radius: 4px;
            border: 2px solid black;
            }
            QPushButton:hover { background-color: #74034e }
            """)
        layout.addWidget(self.poly_color_button)


        self.distortion_path_button = QPushButton("Distortion Path")
        self.distortion_path_button.setVisible(False)
        self.distortion_path_button.setStyleSheet("""
            QPushButton {
                background-color: #b45309;
                color: white;
                font-weight: bold;
                border-radius: 4px;
                border: 2px solid black;
            }
            QPushButton:hover { background-color: #92400e }
            """)
        layout.addWidget(self.distortion_path_button)

    def update_results(self, shape_results, atom_symbol, atom_idx, cn):
        self.header_label.setText(f"Atom: {atom_symbol}{atom_idx}  |  CN: {cn}")
        self.table.setRowCount(len(shape_results))
        for row, (name, value) in enumerate(shape_results.items()):
            item_name = QTableWidgetItem(name)
            item_name.setTextAlignment(Qt.AlignCenter)
            item_name.setFont(QFont("", -1, QFont.Bold))

            item_value = QTableWidgetItem(f"{value:.4f}")
            item_value.setTextAlignment(Qt.AlignCenter | Qt.AlignCenter)

            self.table.setItem(row, 0, item_name)
            self.table.setItem(row, 1, item_value)

        self.table.resizeRowsToContents()
        self.show()
        self.raise_()
        self.show_poly_button.setVisible(True)
        self.distortion_path_button.setVisible(True)
        
    def get_selected_shape(self):
        selected = self.table.selectedItems()
        if not selected:
            return None
        row = self.table.currentRow()
        return self.table.item(row, 0).text()