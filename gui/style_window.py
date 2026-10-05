from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QColorDialog, QSlider, QComboBox
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor


class StyleWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.main = parent
        self.setWindowTitle("Atom")
        self.resize(320, 340)
        self.setStyleSheet("""
            QDialog { background-color: #f8f9fa; }
            QGroupBox {
                font-weight: 600; font-size: 13px; color: #374151;
                border: 1px solid #e5e7eb; border-radius: 10px;
                margin-top: 10px; padding: 12px 10px 10px 10px;
                background-color: #ffffff;
            }
            QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 4px; }
            QLabel { color: #1f2937; font-size: 12px; }
            QPushButton {
                background-color: #ffffff; color: black; font-size: 13px;
                font-weight: 600; border: 1px solid #e5e7eb; border-radius: 8px;
                padding: 7px 12px;
            }
            QPushButton:hover { background-color: #f3f4f6; }
            QComboBox {
                background-color: #ffffff;
                color: #1f2937;
                border: 1px solid #e5e7eb;
                border-radius: 6px;
                padding: 5px 8px;
                font-size: 13px;
            }
            QComboBox QAbstractItemView {
                background-color: #ffffff;
                color: #1f2937;
                selection-background-color: #f0f4ff;
                selection-color: #1f2937;
            }
        """)

        layout = QVBoxLayout(self)

        box = QGroupBox("By element")
        v = QVBoxLayout(box)

        v.addWidget(QLabel("Element:"))
        self.elem_combo = QComboBox()
        self.elem_combo.currentTextChanged.connect(self.on_element_changed)
        v.addWidget(self.elem_combo)

        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Color:"))
        self.color_btn = QPushButton("Choose color")
        self.color_btn.clicked.connect(self.choose_color)
        color_row.addWidget(self.color_btn)
        v.addLayout(color_row)

        v.addWidget(QLabel("Size:"))
        self.size_slider = QSlider(Qt.Horizontal)
        self.size_slider.setMinimum(5)
        self.size_slider.setMaximum(60)
        self.size_slider.valueChanged.connect(self.on_size_changed)
        v.addWidget(self.size_slider)

        layout.addWidget(box)


        box_bg = QGroupBox("Background")
        h = QHBoxLayout(box_bg)
        self.bg_color_btn = QPushButton("Choose background color")
        self.bg_color_btn.clicked.connect(self.choose_bg_color)
        h.addWidget(self.bg_color_btn)
        layout.addWidget(box_bg)

        layout.addStretch()

    def populate(self):
        self.elem_combo.blockSignals(True)
        self.elem_combo.clear()
        elements = sorted(set(self.main.atom_types))
        self.elem_combo.addItems(elements)
        self.elem_combo.blockSignals(False)
        if elements:
            self.on_element_changed(elements[0])

    def on_element_changed(self, elem):
        if not elem:
            return
        size = self.main.get_element_size(elem)
        self.size_slider.blockSignals(True)
        self.size_slider.setValue(int(size * 100))
        self.size_slider.blockSignals(False)

    def choose_color(self):
        elem = self.elem_combo.currentText()
        if not elem:
            return
        qcolor = QColorDialog.getColor(parent=self, title=f"Color for {elem}")
        if qcolor.isValid():
            rgb = (qcolor.red()/255, qcolor.green()/255, qcolor.blue()/255)
            self.main.set_element_color(elem, rgb)

    def choose_bg_color(self):
        qcolor = QColorDialog.getColor(parent=self, title="Background color")
        if qcolor.isValid():
            self.main.set_background_color(qcolor.name())   

    def on_size_changed(self, value):
        elem = self.elem_combo.currentText()
        if elem:
            self.main.set_element_size(elem, value / 100.0)





class StyleWindow2(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.main = parent
        self.setWindowTitle("Bond")
        self.resize(320, 340)
        self.setStyleSheet("""
            QDialog { background-color: #f8f9fa; }
            QGroupBox {
                font-weight: 600; font-size: 13px; color: #374151;
                border: 1px solid #e5e7eb; border-radius: 10px;
                margin-top: 10px; padding: 12px 10px 10px 10px;
                background-color: #ffffff;
            }
            QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 4px; }
            QLabel { color: #1f2937; font-size: 12px; }
            QPushButton {
                background-color: #ffffff; color: black; font-size: 13px;
                font-weight: 600; border: 1px solid #e5e7eb; border-radius: 8px;
                padding: 7px 12px;
            }
            QPushButton:hover { background-color: #f3f4f6; }
            QComboBox {
                background-color: #ffffff;
                color: #1f2937;
                border: 1px solid #e5e7eb;
                border-radius: 6px;
                padding: 5px 8px;
                font-size: 13px;
            }
            QComboBox QAbstractItemView {
                background-color: #ffffff;
                color: #1f2937;
                selection-background-color: #f0f4ff;
                selection-color: #1f2937;
            }
        """)

        layout = QVBoxLayout(self)

        box = QGroupBox("All bonds")
        v = QVBoxLayout(box)

        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Color:"))
        self.color_btn_bond = QPushButton("Choose color")
        self.color_btn_bond.clicked.connect(self.choose_color_bond)
        color_row.addWidget(self.color_btn_bond)
        v.addLayout(color_row)

        v.addWidget(QLabel("Radius:"))
        self.radius_slider = QSlider(Qt.Horizontal)
        self.radius_slider.setMinimum(5)
        self.radius_slider.setMaximum(60)
        self.radius_slider.valueChanged.connect(self.on_radius_changed)
        v.addWidget(self.radius_slider)

        layout.addWidget(box)
        layout.addStretch()


    def populate(self):
        current = getattr(self.main, 'bond_radius', 0.1)
        self.radius_slider.blockSignals(True)
        self.radius_slider.setValue(int(current * 100))
        self.radius_slider.blockSignals(False)

    def choose_color_bond(self):
        qcolor = QColorDialog.getColor(parent=self, title=f"Bond color")
        if qcolor.isValid():
            rgb = (qcolor.red()/255, qcolor.green()/255, qcolor.blue()/255)
            self.main.set_bond_color(rgb) 

    def on_radius_changed(self, value):
            self.main.set_bond_radius(value / 100.0)