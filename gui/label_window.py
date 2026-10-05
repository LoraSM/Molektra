from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QListWidget, QListWidgetItem, QTabWidget, QWidget
)
from PyQt5.QtCore import Qt


class LabelWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.main = parent
        self.setWindowTitle("Atom Labels")
        self.resize(350, 460)
        self.setStyleSheet("""
            QDialog { background-color: #f8f9fa; }
            QLabel { color: #1f2937; font-size: 12px; }
            QPushButton {
                background-color: #ffffff; color: black; font-size: 13px;
                font-weight: 600; border: 1px solid #e5e7eb; border-radius: 8px;
                padding: 7px 12px;
            }
            QPushButton:hover { background-color: #f3f4f6; }
            QLineEdit {
                background-color: #ffffff; color: #1f2937;
                border: 1px solid #e5e7eb; border-radius: 6px; padding: 6px;
            }
            QListWidget {
                background-color: #ffffff; color: #1f2937;
                border: 1px solid #e5e7eb; border-radius: 6px;
            }
            QListWidget::item:selected {
                background-color: #f0f4ff; color: #1f2937;
            }
            QTabWidget::pane {
                border: 1px solid #e5e7eb; border-radius: 6px;
                background-color: #ffffff;
            }
            QTabBar::tab {
                color: #1f2937; padding: 8px 17px;
                background-color: #f3f4f6;
                border: 1px solid #e5e7eb; border-bottom: none;
                border-top-left-radius: 4px; border-top-right-radius: 4px;
            }
            QTabBar::tab:selected {
                background-color: #ffffff; 
            }
        """)

        layout = QVBoxLayout(self)
        tabs = QTabWidget()
        layout.addWidget(tabs)


        tab_manual = QWidget()
        m_layout = QVBoxLayout(tab_manual)
        m_layout.addWidget(QLabel("Indices (e.g. 0, 5, 12-18):"))
        self.index_input = QLineEdit()
        self.index_input.setPlaceholderText("0, 3, 7-10")
        m_layout.addWidget(self.index_input)
        apply_manual = QPushButton("Apply")
        apply_manual.clicked.connect(self.apply_manual)
        m_layout.addWidget(apply_manual)
        m_layout.addStretch()
        tabs.addTab(tab_manual, "By index")


        tab_list = QWidget()
        l_layout = QVBoxLayout(tab_list)
        self.atom_list = QListWidget()
        l_layout.addWidget(self.atom_list)
        apply_list = QPushButton("Apply selection")
        apply_list.clicked.connect(self.apply_list)
        l_layout.addWidget(apply_list)
        tabs.addTab(tab_list, "By atom")


        tab_elem = QWidget()
        e_layout = QVBoxLayout(tab_elem)
        self.elem_list = QListWidget()
        e_layout.addWidget(self.elem_list)
        apply_elem = QPushButton("Apply elements")
        apply_elem.clicked.connect(self.apply_elements)
        e_layout.addWidget(apply_elem)
        tabs.addTab(tab_elem, "By element")

        btn_row = QHBoxLayout()
        show_all = QPushButton("Show all labels")
        show_all.clicked.connect(self.show_all)
        clear = QPushButton("Clear labels")
        clear.clicked.connect(self.clear_labels)
        btn_row.addWidget(show_all)
        btn_row.addWidget(clear)
        layout.addLayout(btn_row)

    def populate(self):
        self.atom_list.clear()
        self.elem_list.clear()
        types = self.main.atom_types

        for i, sym in enumerate(types):
            item = QListWidgetItem(f"{sym}{i}")
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Unchecked)
            item.setData(Qt.UserRole, i)
            self.atom_list.addItem(item)

        for elem in sorted(set(types)):
            item = QListWidgetItem(elem)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Unchecked)
            self.elem_list.addItem(item)

    def apply_manual(self):
        text = self.index_input.text().strip()
        indices = []
        n = len(self.main.atom_types)
        for part in text.split(','):
            part = part.strip()
            if '-' in part:
                try:
                    a, b = part.split('-')
                    indices.extend(range(int(a), int(b) + 1))
                except ValueError:
                    pass
            elif part.isdigit():
                indices.append(int(part))
        indices = [i for i in indices if 0 <= i < n]
        self.main.set_labeled_atoms(indices)

    def apply_list(self):
        indices = []
        for row in range(self.atom_list.count()):
            item = self.atom_list.item(row)
            if item.checkState() == Qt.Checked:
                indices.append(item.data(Qt.UserRole))
        self.main.set_labeled_atoms(indices)

    def apply_elements(self):
        chosen = set()
        for row in range(self.elem_list.count()):
            item = self.elem_list.item(row)
            if item.checkState() == Qt.Checked:
                chosen.add(item.text())
        indices = [i for i, s in enumerate(self.main.atom_types) if s in chosen]
        self.main.set_labeled_atoms(indices)

    def show_all(self):
        n = len(self.main.atom_types)
        if n > 300:
            from PyQt5.QtWidgets import QMessageBox
            reply = QMessageBox.question(
                self, "Many atoms",
                f"This molecule has {n} atoms. Showing all labels may be slow. Continue?",
                QMessageBox.Yes | QMessageBox.No
            )
            if reply != QMessageBox.Yes:
                return
        self.main.set_labeled_atoms(list(range(n)))

    def clear_labels(self):
        self.main.set_labeled_atoms([])