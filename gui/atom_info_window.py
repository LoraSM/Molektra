# Copyright (C) 2026 MOLEKTRA
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QHeaderView,
    QLabel, QFrame, QTableWidget, QTableWidgetItem, QWidget
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QColor
from gui.atom_colors import atom_color_dict
from matplotlib import colors as mcolors

ELEMENT_NAMES = {
    "H": "Hydrogen", "He": "Helium", "Li": "Lithium", "Be": "Beryllium",
    "B": "Boron", "C": "Carbon", "N": "Nitrogen", "O": "Oxygen",
    "F": "Fluorine", "Ne": "Neon", "Na": "Sodium", "Mg": "Magnesium",
    "Al": "Aluminium", "Si": "Silicon", "P": "Phosphorus", "S": "Sulfur",
    "Cl": "Chlorine", "Ar": "Argon", "K": "Potassium", "Ca": "Calcium",
    "Ti": "Titanium","V": "Vanadium","Cr": "Chromium","Mn": "Manganese",
    "Fe": "Iron", "Co":"Cobalt", "Ni":"Nickel", "Cu": "Copper", "Zn": "Zinc",
    "Ga":"Gallium", "Ge": "Germanium", "As": "Arsenic", "Se": "Selenium",
    "Br": "Bromine", "Kr": "Krypton", "Rb": "Rubidium", "Sr": "Strontium",
    "Y": "Yttrium", "Zr": "Zirconium", "Nb": "Niobium", "Mo": "Molybdenum",
    "Tc": "Technetium", "Ru": "Ruthenium", "Rh": "Rhodium", "Pd": "Palladium",
    "Ag": "Silver", "Cd": "Cadmium", "In": "Indium", "Sn": "Tin", "Sb": "Antimony",
    "Te": "Tellurium", "I": "Iodine", "Xe": "Xenon", "Cs": "Cesium", "Ba": "Barium",
    "La": "Lanthanum", "Ce": "Cerium", "Pr": "Praseodymium", "Nd": "Neodymium",
    "Pm": "Promethium", "Sm": "Samarium", "Eu": "Europium", "Gd": "Gadolinium", "Tb": "Terbium",
    "Dy": "Dysprosium", "Ho": "Holmium", "Er": "Erbium", "Tm": "Thulium", "Yb": "Ytterbium","Lu": "Lutetium",
    "Hf": "Hafnium", "Ta": "Tantalum", "W": "Tungsten", "Re": "Rhenium",
    "Os": "Osmium", "Ir": "Iridium", "Pt": "Platinum", "Au": "Gold", "Hg": "Mercury", "Tl": "Thallium",
    "Pb": "Lead", "Bi": "Bismuth", "Po": "Polonium", "At": "Astatine", "Rn": "Radon",
    "Fr": "Francium", "Ra": "Radium", "Ac": "Actinium", "Th": "Thorium", "Pa": "Protactinium","U": "Uranium",
    "Np": "Neptunium", "Pu": "Plutonium", "Am": "Americium", "Cm": "Curium","Bk": "Berkelium",
    "Cf": "Californium", "Es": "Einsteinium", "Fm": "Fermium", "Md": "Mendelevium","No": "Nobelium", "Lr": "Lawrencium",
    "Rf": "Rutherfordium", "Db": "Dubnium", "Sg": "Seaborgium", "Bh": "Bohrium", "Hs": "Hassium", "Mt": "Meitnerium", "Ds": "Darmstadtium",
    "Rg": "Roentgenium", "Cn": "Copernicium", "Nh": "Nihonium", "Fl": "Flerovium", "Mc": "Moscovium", "Lv": "Livermorium", "Ts": "Tennessine", "Og": "Oganesson"
}

def get_element_category(symbol, props):
    group = props.get('group')
    period = props.get('period', 0)
    an = props.get('atomic_number', 0)

    if symbol in ('H',):
        return "Non-metal"
    if group == 18:
        return "Noble gas"
    if group == 1:
        return "Alkali metal"
    if group == 2:
        return "Alkaline earth"
    if group in (3,4,5,6,7,8,9,10,11,12):
        return "Transition metal"
    if symbol in ('C','N','O','P','S','Se'):
        return "Non-metal"
    if symbol in ('F','Cl','Br','I','At'):
        return "Halogen"
    if symbol in ('B','Si','Ge','As','Sb','Te','Po'):
        return "Metalloid"
    if 57 <= an <= 71:
        return "Lanthanide"
    if 89 <= an <= 103:
        return "Actinide"
    if group in (13,14,15,16):
        return "Post-transition"
    return "Metal"

CATEGORY_COLORS = {
    "Non-metal":       ("#fce8e8", "#c0392b", "#f5b7b1"),
    "Noble gas":       ("#e8f4fc", "#1a6fa0", "#b5d4f4"),
    "Alkali metal":    ("#f3e8fc", "#6c3483", "#d7a3f5"),
    "Alkaline earth":  ("#e8fce8", "#1a7a1a", "#a3f5a3"),
    "Transition metal":("#f0f0f0", "#555555", "#cccccc"),
    "Halogen":         ("#e8fcf0", "#0e6655", "#a3f5c8"),
    "Metalloid":       ("#fef9e7", "#7d6608", "#f9e79f"),
    "Lanthanide":      ("#fde8f5", "#8e2461", "#f5a3d4"),
    "Actinide":        ("#fde8e8", "#922b21", "#f5a3a3"),
    "Post-transition": ("#f5f5f5", "#444444", "#cccccc"),
    "Metal":           ("#eaf4fb", "#1a5276", "#85c1e9"),
}


class AtomInfoWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Atom Information")
        self.setMinimumSize(460, 520)
        self.setWindowFlags(Qt.Window | Qt.WindowStaysOnTopHint)
        self.current_element_symbol = None
        self.element_color = "#e03b3b"
        self._build_ui()

    def _build_ui(self):
        main = QVBoxLayout(self)
        main.setContentsMargins(20, 20, 20, 20)
        main.setSpacing(14)

        header_row = QHBoxLayout()
        header_row.setSpacing(14)
        header_row.setAlignment(Qt.AlignVCenter)

        self.symbol_label = QLabel("?")
        self.symbol_label.setFixedSize(72, 72)
        self.symbol_label.setAlignment(Qt.AlignCenter)
        self._update_symbol_style(self.element_color)

        name_col = QVBoxLayout()
        name_col.setSpacing(4)

        self.name_label = QLabel("Unknown")
        self.name_label.setStyleSheet("font-size: 20px; font-weight: 600; color: #111827;")

        self.index_label = QLabel("Index: —")
        self.index_label.setStyleSheet("font-size: 11px; color: #9ca3af; font-style: italic;")

        self.category_tag = QLabel("—")
        self.category_tag.setStyleSheet(
            "font-size: 11px; padding: 2px 10px; border-radius: 10px;"
            "background: #fce8e8; color: #c0392b; border: 1px solid #f5b7b1;"
        )
        self.category_tag.setSizePolicy(self.category_tag.sizePolicy().Preferred,
                                        self.category_tag.sizePolicy().Fixed)

        tag_row = QHBoxLayout()
        tag_row.addWidget(self.category_tag)
        tag_row.addStretch()

        name_col.addWidget(self.name_label)
        name_col.addWidget(self.index_label)
        name_col.addLayout(tag_row)

        header_row.addWidget(self.symbol_label)
        header_row.addLayout(name_col)
        header_row.addStretch()
        main.addLayout(header_row)


        stats_row = QHBoxLayout()
        stats_row.setSpacing(8)

        self.atomic_number_card = self._stat_card("ATOMIC NUMBER", "—")
        self.cov_radius_card     = self._stat_card("COVALENT RADIUS", "—")
        self.cn_card             = self._stat_card("COORDINATION", "—")

        stats_row.addWidget(self.atomic_number_card["widget"])
        stats_row.addWidget(self.cov_radius_card["widget"])
        stats_row.addWidget(self.cn_card["widget"])
        main.addLayout(stats_row)

        coords_section_label = QLabel("COORDINATES")
        coords_section_label.setStyleSheet(
            "font-size: 10px; font-weight: 600; color: #9ca3af; letter-spacing: 1px;"
        )
        main.addWidget(coords_section_label)

        coords_row = QHBoxLayout()
        coords_row.setSpacing(8)
        self.x_label = self._coord_cell("X", "#e03b3b", coords_row)
        self.y_label = self._coord_cell("Y", "#3b82f6", coords_row)
        self.z_label = self._coord_cell("Z", "#10b981", coords_row)
        main.addLayout(coords_row)

        info_section_label = QLabel("ADDITIONAL INFORMATION")
        info_section_label.setStyleSheet(
            "font-size: 10px; font-weight: 600; color: #9ca3af; letter-spacing: 1px;"
        )
        main.addWidget(info_section_label)

        self.info_table = QTableWidget()
        self.info_table.setColumnCount(2)
        self.info_table.setHorizontalHeaderLabels(["Property", "Value"])
        self.info_table.verticalHeader().setVisible(False)
        self.info_table.setShowGrid(False)
        self.info_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.info_table.setAlternatingRowColors(False)
        self.info_table.setMaximumHeight(160)
        self.info_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.info_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.info_table.setStyleSheet("""
            QTableWidget {
                background-color: #ffffff;
                border: 1px solid #e5e7eb;
                border-radius: 8px;
                font-size: 13px;
                outline: none;
            }
            QTableWidget::item {
                padding: 8px 14px;
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
                padding: 7px 14px;
                border: none;
                border-bottom: 1px solid #e5e7eb;
            }
        """)
        main.addWidget(self.info_table)
        main.addStretch()


    def _update_symbol_style(self, color):
        darker = self._darken_color(color, 15)
        qc = QColor(color)
        luminance = 0.299 * qc.red() + 0.587 * qc.green() + 0.114 * qc.blue()
        text_col = "black" if luminance > 140 else "white"
        self.symbol_label.setStyleSheet(f"""
            background-color: {color};
            color: {text_col};
            border-radius: 36px;
            font-size: 26px;
            font-weight: bold;
        """)

    def _stat_card(self, label_text, value_text):
        """Returns a dict with the card widget and references to update labels."""
        card = QFrame()
        card.setStyleSheet("""
            QFrame {
                background-color: #f9fafb;
                border-radius: 10px;
                border: 1px solid #e5e7eb;
            }
        """)
        layout = QVBoxLayout(card)
        layout.setContentsMargins(10, 12, 10, 12)
        layout.setSpacing(4)
        layout.setAlignment(Qt.AlignCenter)

        lbl = QLabel(label_text)
        lbl.setAlignment(Qt.AlignCenter)
        lbl.setStyleSheet("font-size: 9px; font-weight: 600; color: #9ca3af; letter-spacing: 0.8px; background: transparent; border: none;")
        lbl.setWordWrap(True)

        val = QLabel(value_text)
        val.setAlignment(Qt.AlignCenter)
        val.setStyleSheet("font-size: 20px; font-weight: 600; color: #111827; background: transparent; border: none;")

        layout.addWidget(lbl)
        layout.addWidget(val)
        return {"widget": card, "value_label": val}

    def _coord_cell(self, axis, color, parent_layout):
        cell = QFrame()
        c = QColor(color)
        cell.setStyleSheet(f"""
            QFrame {{
        background: white;
        border: 1px solid rgba({c.red()}, {c.green()}, {c.blue()}, 85);
        border-radius: 8px;
    }}
        """)
        cl = QVBoxLayout(cell)
        cl.setContentsMargins(10, 8, 10, 8)
        cl.setSpacing(2)

        axis_lbl = QLabel(axis)
        axis_lbl.setAlignment(Qt.AlignLeft)
        axis_lbl.setStyleSheet(f"color: {color}; font-size: 10px; font-weight: 700; background: transparent; border: none;")

        val_lbl = QLabel("—")
        val_lbl.setAlignment(Qt.AlignLeft)
        val_lbl.setFont(QFont("Courier New", 11, QFont.Bold))
        val_lbl.setStyleSheet("color: #111827; background: transparent; border: none;")

        cl.addWidget(axis_lbl)
        cl.addWidget(val_lbl)
        parent_layout.addWidget(cell)
        return val_lbl

    def _darken_color(self, hex_color, percent=15):
        try:
            h = hex_color.lstrip('#')
            rgb = tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
            d = tuple(max(0, int(c * (1 - percent / 100))) for c in rgb)
            return '#{:02x}{:02x}{:02x}'.format(*d)
        except Exception:
            return hex_color


    def update_atom(self, symbol, index, x, y, z, cn=None, types=None):
        self.current_element_symbol = symbol
        name = ELEMENT_NAMES.get(symbol, symbol)

        props = atom_color_dict.get(symbol, {'color': '#e03b3b'})
        color_raw = props.get('color', '#e03b3b')
        try:
            if isinstance(color_raw, str) and color_raw.startswith('#'):
                self.element_color = color_raw
            else:
                rgba = mcolors.to_rgba(color_raw)
                self.element_color = '#{:02x}{:02x}{:02x}'.format(
                    int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255)
                )
        except Exception:
            self.element_color = '#e03b3b'

        # Symbol circle
        self._update_symbol_style(self.element_color)
        self.symbol_label.setText(symbol)

        # Name + index
        self.name_label.setText(name)
        self.index_label.setText(f"Index: {index}")

        # Category tag
        category = get_element_category(symbol, props)
        bg, fg, border = CATEGORY_COLORS.get(category, ("#f0f0f0", "#555555", "#cccccc"))
        self.category_tag.setText(category)
        self.category_tag.setStyleSheet(
            f"font-size: 11px; padding: 2px 10px; border-radius: 10px;"
            f"background: {bg}; color: {fg}; border: 1px solid {border};"
        )

        # Stat cards
        an = props.get('atomic_number', '—')
        self.atomic_number_card["value_label"].setText(str(an))

        rad = props.get('cov_radius', None)
        rad_text = f"{rad:.2f} Å" if rad is not None else "—"
        self.cov_radius_card["value_label"].setText(rad_text)
        self.cov_radius_card["value_label"].setStyleSheet(
            "font-size: 17px; font-weight: 600; color: #111827; background: transparent; border: none;"
        )

        cn_text = str(cn) if cn is not None else "—"
        self.cn_card["value_label"].setText(cn_text)

        # Coordinates
        self.x_label.setText(f"{x:.4f}")
        self.y_label.setText(f"{y:.4f}")
        self.z_label.setText(f"{z:.4f}")

        # Additional info table
        self._update_info_table(symbol, props)

        self.show()
        self.raise_()
        self.move(20, 20)

    def _update_info_table(self, symbol, props):
        self.info_table.setRowCount(0)
        rows = []

        en = props.get('electronegativity', None)
        if en is not None:
            rows.append(("Electronegativity", f"{en:.2f}"))

        config = props.get('electron_config', None)
        if config:
            rows.append(("Electron config", config))

        oxidation = props.get('oxidation_state', None)
        if oxidation:
            rows.append(("Common ox. state", str(oxidation)))

        density = props.get('density', None)
        if density:
            rows.append(("Density (g/cm³)", f"{density:.5f}"))

        melting = props.get('melting_point', None)
        if melting:
            rows.append(("Melting point (K)", str(melting)))

        vdw_radius = props.get('vdw', None)
        if vdw_radius:
            rows.append(("Van der Waals radius (Å)", str(vdw_radius)))

        self.info_table.setRowCount(len(rows))
        for r, (key, val) in enumerate(rows):
            key_item = QTableWidgetItem(key)
            key_item.setFlags(key_item.flags() & ~Qt.ItemIsEditable)
            key_item.setTextAlignment(Qt.AlignVCenter | Qt.AlignLeft)

            val_item = QTableWidgetItem(val)
            val_item.setFlags(val_item.flags() & ~Qt.ItemIsEditable)
            val_item.setTextAlignment(Qt.AlignVCenter | Qt.AlignRight)
            val_item.setForeground(Qt.gray)

            self.info_table.setItem(r, 0, key_item)
            self.info_table.setItem(r, 1, val_item)
            self.info_table.setRowHeight(r, 38)

    def update_shapes(self, shape_results: dict):
        pass

    def clear_shapes(self):
        pass

