# Copyright (C) 2026 MOLEKTRA
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

from PyQt5.QtWidgets import QDialog, QGridLayout, QVBoxLayout, QLabel, QPushButton
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor
from matplotlib import colors as mcolors
from gui.atom_colors import atom_color_dict

Lanthanide_row = 8

def _hex(color):
    try:
        return mcolors.to_hex(color)
    except ValueError:
        return "#000000"  

def _text_color(hex_color):
    c = QColor(hex_color)
    luminance = (0.299 * c.red() + 0.587 * c.green() + 0.114 * c.blue()) 
    return 'black' if luminance > 140 else 'white'

class PeriodicTableDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Add atom")
        self.selected_symbol = None

        main = QVBoxLayout(self)
        main.setContentsMargins(18,18,18,18)
        main.setSpacing(10)
        title = QLabel("Choose an element to add")
        title.setStyleSheet("font-size: 16px; font-weight: bold; color: #333;")
        main.addWidget(title)

        grid = QGridLayout()
        grid.setSpacing(3)
        main.addLayout(grid)
        lanthanide_col = 3
        for sym, props in atom_color_dict.items():
            group, period = props.get('group'), props.get('period')
            if group == 'Lanthanides':
                row, col= Lanthanide_row, lanthanide_col
                lanthanide_col += 1
            elif isinstance(group, int) and isinstance(period, int):
                row, col = period, group
            else:
                continue
            grid.addWidget(self.element_button(sym, props), row, col)

        marker = QLabel("La-Lu")
        marker.setAlignment(Qt.AlignCenter)
        marker.setStyleSheet("font-size: 14px; font-weight: bold; color: #333;")
        grid.addWidget(marker,6,3)
        grid.setRowMinimumHeight(6,9)

    def element_button(self, symbol, props):
        bg = _hex(props.get('color', 'grey'))
        fg = _text_color(bg)
        btn = QPushButton(symbol)
        btn.setFixedSize(50,50)
        btn.setCursor(Qt.PointingHandCursor)
        btn.setToolTip(f"{symbol}\n{props.get('name', '')}\n{props.get('atomic_number', '')}")
        btn.setStyleSheet(f"background-color: {bg}; color: {fg}; font-weight: bold; border: 1px solid black;")
        btn.clicked.connect(lambda _, s=symbol: self.get_position(s))
        return btn

    def get_position(self, symbol):
        self.selected_symbol = symbol
        self.accept()