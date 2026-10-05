# Copyright (C) 2026 MOLEKTRA
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QTableWidget,QPushButton,
    QTableWidgetItem, QHeaderView, QAbstractItemView,QColorDialog
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QColor 
from gui.shape_data import SHAPE_FULL_NAMES

class ShapeResultsWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Shape Measures")
        self.setMinimumSize(500, 500)
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
                border: 1px solid #e5e7eb;
                border-radius: 8px;
                font-size: 12px;
                gridline-color: transparent;
                outline: none;
            }
            QTableWidget::item {
                padding: 8px 12px;
                border: none;
                color: #1f2937;
                border-bottom: 1px solid #f3f4f6;
            }
            QTableWidget::item:selected {
                background-color: #f0f4ff;
                color: #1f2937;
            }
            QHeaderView::section {
                background-color: #f9fafb;
                color: #6b7280;
                font-weight: 600;
                font-size: 12px;
                padding: 7px 12px;
                border: none;
                border-bottom: 1px solid #e5e7eb;
            }
        """)
        hdr = self.table.horizontalHeader()
        hdr.setSectionResizeMode(0, QHeaderView.Stretch)
        hdr.setSectionResizeMode(1, QHeaderView.Stretch)
        layout.addWidget(self.table)

        
        self.show_poly_button = QPushButton("Show Polyhedra")
        self.show_poly_button.setVisible(False)
        self.show_poly_button.setStyleSheet("""
            QPushButton {
                background-color: #D0F5D1;
                color: black;
                font-size: 13px;
                font-weight: 600;
                border: 1px solid #e5e7eb;
                border-radius: 8px;
                padding: 7px 12px;
            }
            QPushButton:hover { background-color: #f3f4f6; }
        """)
        layout.addWidget(self.show_poly_button)

        
        self.poly_color_button = QPushButton("Set Polyhedron Color")
        self.poly_color_button.setVisible(False)
        self.poly_color_button.setStyleSheet("""
            QPushButton {
                background-color: #ffffff;
                color: #1f2937;
                font-size: 13px;
                font-weight: 600;
                border: 1px solid #e5e7eb;
                border-radius: 8px;
                padding: 7px 12px;
            }
            QPushButton:hover { background-color: #f3f4f6; }
        """)
        layout.addWidget(self.poly_color_button)


        self.distortion_path_button = QPushButton("Distortion Path")
        self.distortion_path_button.setVisible(False)
        self.distortion_path_button.setStyleSheet("""
            QPushButton {
                background-color: #FBE1BC;
                color: black;
                font-size: 13px;
                font-weight: 600;
                border: 1px solid #e5e7eb;
                border-radius: 8px;
                padding: 7px 12px;
            }
            QPushButton:hover { background-color: #f3f4f6; }
        """)
        layout.addWidget(self.distortion_path_button)

        self.custom_shape_button = QPushButton("Load custom polyhedron")
        self.custom_shape_button.setVisible(False)
        self.custom_shape_button.setStyleSheet("""
            QPushButton {
                background-color: #B8FFE1;
                color: black;
                font-size: 13px;
                font-weight: 600;
                border: 1px solid #e5e7eb;
                border-radius: 8px;
                padding: 7px 12px;
            }
            QPushButton:hover { background-color: #f3f4f6; }
        """)
        self.custom_shape_button.clicked.connect(self.on_custom_clicked)
        layout.addWidget(self.custom_shape_button)
        
    def add_result_row(self, name, value):
        row = self.table.rowCount()
        self.table.insertRow(row)

        item_name = QTableWidgetItem(name)
        item_name.setTextAlignment(Qt.AlignCenter)
        item_name.setFont(QFont("", -1, QFont.Bold))
        item_name.setBackground(QColor("#f0f4ff"))

        item_value = QTableWidgetItem(f"{value:.4f}")
        item_value.setTextAlignment(Qt.AlignCenter)
        item_value.setBackground(QColor("#f0f4ff"))

        self.table.setItem(row, 0, item_name)
        self.table.setItem(row, 1, item_value)
        self.table.resizeRowsToContents()


    def update_results(self, shape_results, atom_symbol, atom_idx, cn):
        self.header_label.setText(f"Atom: {atom_symbol}{atom_idx}  |  CN: {cn}")
        sorted_results = sorted(shape_results.items(), key=lambda x: x[1])
        self.table.setRowCount(len(sorted_results))
        for row, (name, value) in enumerate(sorted_results):
            item_name = QTableWidgetItem(name)
            item_name.setTextAlignment(Qt.AlignCenter)
            item_name.setFont(QFont("", -1, QFont.Bold))
            item_name.setToolTip(SHAPE_FULL_NAMES.get(name, name))

            item_value = QTableWidgetItem(f"{value:.4f}")
            item_value.setTextAlignment(Qt.AlignCenter | Qt.AlignCenter)

            self.table.setItem(row, 0, item_name)
            self.table.setItem(row, 1, item_value)

        self.table.resizeRowsToContents()
        self.show()
        self.raise_()
        self.show_poly_button.setVisible(True)
        self.distortion_path_button.setVisible(True)
        self.custom_shape_button.setVisible(True)

    def on_custom_clicked(self):
        main = self.parent()
        if main is not None:
            main.load_custom_polyhedron()
        
    def get_selected_shape(self):
        selected = self.table.selectedItems()
        if not selected:
            return None
        row = self.table.currentRow()
        return self.table.item(row, 0).text()