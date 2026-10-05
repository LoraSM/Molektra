# Copyright (C) 2026 MOLEKTRA
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import numpy as np, os, pyqtgraph as pg , csv
import io, matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from cosymlib import Cosymlib
from PyQt5.QtWidgets import (
     QMainWindow, QWidget, QPushButton,QTableWidget,QMessageBox,QHeaderView,QMenu,QLabel,
     QVBoxLayout, QHBoxLayout, QFileDialog,QTableWidgetItem,QColorDialog,QDialog,QSizePolicy
)
from PyQt5.QtCore import QCoreApplication, QUrl,Qt
QCoreApplication.processEvents()

from PyQt5.QtCore import Qt,QSize,QPoint
from PyQt5.QtGui import QPalette, QColor, QVector4D,QDesktopServices,QPixmap,QIcon,QPainter, QPen
from cosymlib import Geometry
from pyqtgraph.opengl import GLViewWidget
from pyqtgraph.opengl import GLMeshItem, MeshData
from pyqtgraph.opengl import GLLinePlotItem
from matplotlib import colors 
from utils.graphics import create_cylinder
from gui.atom_colors import atom_color_dict
from gui.atom_info_window import AtomInfoWindow
from gui.shape_data import SHAPE_REFERENCES
from gui.coordination_geometry import SHAPE_SYMMETRIES, POINT_GROUPS_TO_TEST
from gui.shape_results_window import ShapeResultsWindow
from posym import SymmetryMolecule
from molecule.loader import load_structure

###For selection 
class SelectionOverlay(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAttribute(Qt.WA_TransparentForMouseEvents) 
        self.mode = None
        self.path = []

    def update_path(self, mode, path):
        self.mode = mode
        self.path = path
        self.update()  

    def clear(self):
        self.mode = None
        self.path = []
        self.update()

    def paintEvent(self, event):
        if not self.path or self.mode is None:
            return
        painter = QPainter(self)
        pen = QPen(QColor(31, 65, 74), 2, Qt.DashLine)
        painter.setPen(pen)

        pts = self.path
        if self.mode == 'rect':
            x0, y0 = pts[0]
            x1, y1 = pts[-1]
            painter.drawRect(int(min(x0, x1)), int(min(y0, y1)),
                             int(abs(x1 - x0)), int(abs(y1 - y0)))
        elif self.mode == 'free':
            for a, b in zip(pts, pts[1:]):
                painter.drawLine(int(a[0]), int(a[1]), int(b[0]), int(b[1]))


#SYMMETRY/SHAPE
def calculate_symmetry_posym(coords, symbols, groups=None):
    if groups is None:
        groups = POINT_GROUPS_TO_TEST
    if isinstance(groups, str):
        groups = [groups]

    results = {}
    operations = []
    center = np.zeros(3)
    coordinates = [list(c) for c in coords]

    for group in groups:
        try:
            sym = SymmetryMolecule(group=group, coordinates=coordinates, symbols=symbols)
            results[group] = round(float(sym.measure_pos), 4)
            operations = sym.get_oriented_operations()
            center = np.array(sym.center)  
        except Exception:
            pass

    return dict(sorted(results.items(), key=lambda x: x[1])), operations, center

def get_ideal_coords(symbols, coords, shape_code, central_atom=1):
    geo = Geometry(symbols=symbols, positions=np.array(coords))
    mol = Cosymlib([geo])

    buffer = io.StringIO()
    mol.print_shape_structure([shape_code], central_atom=central_atom, output=buffer)
    output = buffer.getvalue()

    ideal_coords = []
    lines = output.strip().split('\n')
    found = False
    for line in lines:
        if shape_code in line:
            found = True
            continue
        if found:
            parts = line.split()
            if len(parts) == 4:
                ideal_coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return np.array(ideal_coords) if ideal_coords else None

def calculate_distortion_path(symbols, coords, shape_code_A, shape_code_B,
                               central_atom=1, n_points=50):
    geo = Geometry(symbols=symbols, positions=np.array(coords))
    mol_sA = geo.get_shape_measure(shape_code_A, central_atom=central_atom)
    mol_sB = geo.get_shape_measure(shape_code_B, central_atom=central_atom)

    ideal_A = get_ideal_coords(symbols, coords, shape_code_A, central_atom)
    ideal_B = get_ideal_coords(symbols, coords, shape_code_B, central_atom)

    if ideal_A is None or ideal_B is None or len(ideal_A) != len(ideal_B):
        return None, None, mol_sA, mol_sB

    path_sA, path_sB = [], []
    for t in np.linspace(0, 1, n_points):
        interp = (1 - t) * ideal_A + t * ideal_B
        try:
            geo_t = Geometry(symbols=symbols, positions=interp)
            path_sA.append(geo_t.get_shape_measure(shape_code_A, central_atom=central_atom))
            path_sB.append(geo_t.get_shape_measure(shape_code_B, central_atom=central_atom))
        except Exception:
            pass

    return path_sA, path_sB, mol_sA, mol_sB

def calculate_point_group(atom_positions_3d,atom_types):
    try:
        symbols   = atom_types
        positions = np.array(atom_positions_3d)
        geo       = Geometry(symbols=symbols, positions= positions)
        point_group = geo.get_pointgroup()
        return point_group
    except Exception as e:
        print(f'Error getting point group: {e}')

def _symmetry_eigenvector(matrix, target_eigenvalue):
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    idx = np.argmin(np.abs(np.real(eigenvalues) - target_eigenvalue))
    return np.real(eigenvectors[:, idx])

def get_symmetry_axis(matrix):
    # rotation axis: eigenvector whose eigenvalue is +1
    return _symmetry_eigenvector(matrix, 1.0)

def get_plane_normal(matrix):
    # mirror-plane normal: eigenvector whose eigenvalue is -1
    return _symmetry_eigenvector(matrix, -1.0)

def calc_distance(p1, p2):
    return np.linalg.norm(p2 - p1)

def calc_angle(p1, p2, p3):
    v1 = p1 - p2
    v2 = p3 - p2
    cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_a = np.clip(cos_a, -1.0, 1.0)
    return np.degrees(np.arccos(cos_a))

def calc_dihedral(p1, p2, p3, p4):
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return np.degrees(np.arctan2(y, x))

def create_dashed_line(p1, p2, n_dashes=20, color=(1, 1, 0, 1)):
    points = []
    for i in range(n_dashes):
        t_start = i / n_dashes
        t_end   = (i + 0.6) / n_dashes  
        points.append(p1 + t_start * (p2 - p1))
        points.append(p1 + t_end   * (p2 - p1))
        points.append([np.nan, np.nan, np.nan])  
    return GLLinePlotItem(
        pos=np.array(points),
        color=color,
        width=5,
        antialias=True,
        mode='lines'
    )

def calculate_shape(central_atom_idx, atom_positions_3d, atom_types):
    central_pos = atom_positions_3d[central_atom_idx]
    central_sym = atom_types[central_atom_idx]
    ligand_indices = get_ligands(central_atom_idx, atom_positions_3d, atom_types)

    print(f"Detected ligands: {len(ligand_indices)}")
    cn = len(ligand_indices)
    coords = [central_pos] + [atom_positions_3d[i] for i in ligand_indices]
    symbols = [central_sym] + [atom_types[i] for i in ligand_indices]

    #-----------GEOMETRÍA Y SHAPE MEASURES----------------
    geo = Geometry(symbols=symbols, positions=coords)
    results = {}
    shapes = SHAPE_REFERENCES.get(cn, {})

    for shape_code, shape_name in shapes.items():
        value = geo.get_shape_measure(shape_code, central_atom=1)
        results[shape_name] = value
        print(f"{shape_name}: {value:.4f}")

    return results,cn

#----------------------------------------------------------------------#

#FUNCTIONS FOR 3D VISUALIZATION
def create_sphere(radius=0.3, color=(1, 0, 0, 1)):
    md = MeshData.sphere(rows=20, cols=40, radius=radius)
    sphere = GLMeshItem(meshdata=md, smooth=True, color=color, shader='shaded', glOptions='opaque')
    return sphere

def covalent_threshold(type1, type2):
    r1 = atom_color_dict.get(type1, {}).get('cov_radius', 1.5)
    r2 = atom_color_dict.get(type2, {}).get('cov_radius', 1.5)
    return (r1 + r2) +0.3

def get_ligands(central_idx, positions, types):
    """Indices of atoms covalently bonded to the central atom.

    Single source of truth for the ligand-detection loop that was previously
    copy-pasted across the shape, polyhedron, symmetry, isolate and
    atom-info routines.
    """
    central_pos = np.asarray(positions[central_idx])
    central_sym = types[central_idx]
    return [
        i for i, (pos, sym) in enumerate(zip(positions, types))
        if i != central_idx
        and np.linalg.norm(np.asarray(pos) - central_pos) < covalent_threshold(central_sym, sym)
    ]






class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.export_data = {}
        self._bond_pairs = []
        self.isolate_mode = False
        self.setWindowTitle("Molektra")
        self.resize(600, 600)
        self.shape_results_window = ShapeResultsWindow(parent=self)
        self.shape_results_window.show_poly_button.clicked.connect(self.toggle_polyhedron)
        self.shape_results_window.distortion_path_button.clicked.connect(self.show_distortion_path)  
        self._polyhedron_items = []

        self.atom_positions_3d = []  
        self.atom_types = []          
        self.selected_atom_index = None  
        self.selected_atoms = []
        self.selection_mode = None
        self._drag_path = []
        
        self.atom_info_window = AtomInfoWindow(parent=self)

        #---------------Layout---------------------#
        container = QWidget()
        main_layout = QHBoxLayout()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

        #---------------3D View ---------------------#
        self.molecule_view = GLViewWidget()
        self.molecule_view.setMinimumSize(400, 600)
        self.molecule_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.molecule_view.opts['distance'] = 20
        self.molecule_view.setBackgroundColor('#f0f0f0')
        palette = self.molecule_view.palette()
        palette.setColor(QPalette.Window, QColor("#e0e0e0"))
        self.molecule_view.setPalette(palette)

        self.molecule_view.mouseMoveEvent = self.on_mouse_move
        self.molecule_view.mousePressEvent = self.on_mouse_press
        self.molecule_view.mouseReleaseEvent = self.on_mouse_release
        self.selection_overlay = SelectionOverlay(self.molecule_view)
        self.selection_overlay.setGeometry(self.molecule_view.rect())
        self.selection_overlay.hide()
        self.last_mouse_pos = None
        self.sym_label_overlay = QLabel("", self.molecule_view)
        self.sym_label_overlay.setStyleSheet("""
        QLabel {
            color: white;
            background-color: rgba(0, 0, 0, 140);
            border-radius: 6px;
            padding: 4px 10px;
            font-size: 14px;
            font-weight: bold;
            }
            """)
        self.sym_label_overlay.move(10, 10)
        self.sym_label_overlay.hide()


        #--------Visualization Toolbar----------------#
        toolbar = QWidget()
        toolbar.setFixedHeight(44)
        toolbar.setStyleSheet("background-color: #ffffff; border-bottom: 1px solid #e5e7eb;")
        toolbar_layout = QHBoxLayout()
        toolbar_layout.setContentsMargins(8, 5, 8, 5)
        toolbar_layout.setSpacing(6)
        toolbar.setLayout(toolbar_layout)

        self.delete_button = self.make_toolbar_button("delete.svg", "Delete selection", self.delete_selection, icon_size=20)
        self.boton_ejes     = self.make_toolbar_button("axis.svg",    "Show axes XYZ", self.toggle_ejes_xyz)
        self.hide_h_button  = self.make_toolbar_button("hydrogen.png","Show/Hide Hydrogens", self.toggle_hydrogens,icon_size=25)
        self.coord_env_button = self.make_toolbar_button("cn.png", "Isolate coordination environment", self.toggle_coord_env, icon_size=25)
        self.selection_button = self.make_toolbar_button("selection.png", "Select atoms", self.open_selection_menu, icon_size=25)
        self.reset_button = self.make_toolbar_button("reset.svg", "Show all", self.restore_original, icon_size=26)
        self.distance_button = self.make_toolbar_button("distance.png", "Measure distance", self.measure_distance, icon_size=25)
        self.angle_button = self.make_toolbar_button("angle.png", "Measure angle", self.measure_angle, icon_size=40)
        self.dihedral_button = self.make_toolbar_button("dihedral.png", "Measure dihedral angle", self.measure_dihedral, icon_size=28)



        
        toolbar_layout.addWidget(self.hide_h_button)
        toolbar_layout.addWidget(self.distance_button)
        toolbar_layout.addWidget(self.angle_button)
        toolbar_layout.addWidget(self.dihedral_button)
        toolbar_layout.addWidget(self.coord_env_button)
        toolbar_layout.addWidget(self.selection_button)
        toolbar_layout.addWidget(self.boton_ejes) 
        toolbar_layout.addWidget(self.delete_button)  
        toolbar_layout.addWidget(self.reset_button)
        toolbar_layout.addStretch()

        viewer_container = QWidget()
        viewer_layout = QVBoxLayout()
        viewer_layout.setContentsMargins(0, 0, 0, 0)
        viewer_layout.setSpacing(0)
        viewer_container.setLayout(viewer_layout)
        viewer_layout.addWidget(toolbar)
        viewer_layout.addWidget(self.molecule_view, stretch=1)

        #---------------Main Toolbar---------------------#
        button_panel = QWidget()
        button_layout = QVBoxLayout()
        button_panel.setLayout(button_layout)

        load_button = QPushButton("Load File")
        load_button.setFixedHeight(35)
        load_button.setStyleSheet("""
                   QPushButton {
                background-color: #5f6064;                       
                color: white;
                font-weight: bold;
                border-radius: 4px;
		border: 2px solid black;  }
            QPushButton:hover {background-color: #5f6064 }                                
                                  """)
        load_button.clicked.connect(self.load_file)
        button_layout.addWidget(load_button)

        #---------------BOTON ATOM---------------------#
        atom_info_button = QPushButton("Atom Information")
        atom_info_button.setFixedHeight(35)
        atom_info_button.setStyleSheet("""
            QPushButton {
                background-color: #383634;
                color:white;                       
                font-weight: bold;
                border-radius: 4px;
                border: 2px solid black; }
            QPushButton: hover {background-color: #383634}
            """   )
        atom_info_button.clicked.connect(self.show_atom_info)
        button_layout.addWidget(atom_info_button)


        self.btn_reset = QPushButton("Clear")
        self.btn_reset.setFixedHeight(35)
        self.btn_reset.setStyleSheet("""
            QPushButton {
            background-color: #383634;
            color: white;
            font-weight: bold;
            border-radius: 4px;
            border: 2px solid black; }
            QPushButton:hover { background-color: #383634 }
        """)
        self.btn_reset.clicked.connect(self.reset_measurements)
        self.btn_reset.setVisible(False)
        button_layout.addWidget(self.btn_reset)




        button_layout.addStretch()
        calc_shape_button = QPushButton("Calculate Shape")
        calc_shape_button.setFixedHeight(35)
        calc_shape_button.setStyleSheet(""" 
            QPushButton {
                background-color: #45af49;                       
                color: white;
                font-weight: bold;
                border-radius: 4px;
		border: 2px solid black;  }
            QPushButton:hover {background-color: #45a049 }
        """)
        calc_shape_button.clicked.connect(self.run_shape_calculation)
        button_layout.addWidget(calc_shape_button)

        #---------------BOTON Point---------------------#
        calc_pointgroup = QPushButton("Get Point Group")
        calc_pointgroup.setFixedHeight(35)
        calc_pointgroup.setStyleSheet(""" 
            QPushButton {
                background-color: #0354bb;                       
                color: white;
                font-weight: bold;
                border-radius: 4px;
		        border: 2px solid black;  }
            QPushButton:hover {background-color: #0354bb}

        """)
        calc_pointgroup.clicked.connect(self.run_point_group)
        button_layout.addWidget(calc_pointgroup)


        #------------------------TABLE----------------------------#
        # Pointgroup #
        self.pointgroup_table = QTableWidget()
        self.pointgroup_table.setColumnCount(1)
        self.pointgroup_table.setHorizontalHeaderLabels(['Point Group'])
        self.pointgroup_table.setMaximumHeight(140)
        self.pointgroup_table.verticalHeader().setVisible(False)
        self.pointgroup_table.setShowGrid(False)
        self.pointgroup_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.pointgroup_table.setAlternatingRowColors(False)
        self.pointgroup_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.pointgroup_table.setStyleSheet("""
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
        self.pointgroup_table.setVisible(False)
        button_layout.addWidget(self.pointgroup_table)

        self.close_pointgroup_button = QPushButton("Close")
        self.close_pointgroup_button.clicked.connect(lambda: (
        self.pointgroup_table.setVisible(False),
        self.close_pointgroup_button.setVisible(False)
        ))
        self.close_pointgroup_button.setStyleSheet(""" 
            QPushButton {
                background-color: #0354bb;                       
                color: white;
                font-weight: bold;
                border-radius: 4px;
		border: 2px solid black;  }
            QPushButton:hover {background-color: #0354bb}

        """)
        self.close_pointgroup_button.setVisible(False)
        button_layout.addWidget(self.close_pointgroup_button)



        #---------------BOTON Point---------------------#

        calc_sym_posym = QPushButton("Symmetry Measures")
        calc_sym_posym.setFixedHeight(35)
        calc_sym_posym.setStyleSheet(""" 
            QPushButton {
                background-color: #74034e;                       
                color: white;
                font-weight: bold;
                border-radius: 4px;
		        border: 2px solid black;  }
            QPushButton:hover {background-color: #74034e}

        """)
        calc_sym_posym.clicked.connect(self.run_symmetry_posym)
        button_layout.addWidget(calc_sym_posym)
        button_layout.addStretch()





        #BOTON RESULTS
        self.export_button = QPushButton("Export")
        self.export_button.setFixedHeight(35)
        self.export_button.setStyleSheet("""
            QPushButton {
                background-color: #5f6064;
                color: white;
                font-weight: bold;
                border-radius: 4px;
                border: 2px solid black; }
            QPushButton:hover { background-color: #5f6064 }
            """)
        self.export_button.clicked.connect(self.toggle_export)
        button_layout.addWidget(self.export_button)

        #BOTON Info 
        self.info_us = QPushButton("About us")
        self.info_us.setFixedHeight(35)
        self.info_us.setStyleSheet("""
            QPushButton {
                background-color: #5f6064;
                color: white;
                font-weight: bold;
                border-radius: 4px;
                border: 2px solid black; }
            QPushButton:hover { background-color: #5f6064 }
            """)
        self.info_us.clicked.connect(self.show_info)
        button_layout.addWidget(self.info_us)


        #BOTON EXIT
        self.exit_button = QPushButton("Exit")
        self.exit_button.setFixedHeight(35)
        self.exit_button.setStyleSheet("""
                   QPushButton {
                background-color: #5f6064;                       
                color: white;
                font-weight: bold;
                border-radius: 4px;
		border: 2px solid black;  }
            QPushButton:hover {background-color: #5f6064 }                                
                                  """)
        self.exit_button.clicked.connect(self.close)
        button_layout.addWidget(self.exit_button)
        button_layout.addStretch()

        main_layout.addWidget(viewer_container, stretch=1)
        main_layout.addWidget(button_panel, stretch=0)


        # ── Ejes XYZ ────────────────────────────────────────────────
        self.ejes_visibles = False
        self.eje_x, self.eje_y, self.eje_z = self.crear_ejes_xyz()
        

    # ────────────────────────────────────────────────────────────────
    # RATÓN
    # ────────────────────────────────────────────────────────────────

    def on_mouse_press(self, event):
        if self.selection_mode in ('rect', 'free'):
            self._drag_path = [np.array([event.pos().x(), event.pos().y()])]
            self.selection_overlay.setGeometry(self.molecule_view.rect())
            self.selection_overlay.show()
            return
        
        self._last_mouse_pos = event.pos()
        modifiers = event.modifiers()
        is_ctrl_pressed = modifiers == Qt.ControlModifier
        self.select_atom_closest(event.pos(), ctrl_pressed = is_ctrl_pressed)
        GLViewWidget.mousePressEvent(self.molecule_view, event)

    def on_mouse_release(self, event):
        if self.selection_mode in ('rect', 'free') and self._drag_path:
            self._finish_area_selection()
            self.selection_overlay.clear()
            self.selection_overlay.hide()
            self.selection_mode = None
            self._drag_path = []
            return
        self._last_mouse_pos = None
    
    def on_mouse_move(self, event):
        if self.selection_mode in ('rect', 'free') and self._drag_path:
            self._drag_path.append(np.array([event.pos().x(), event.pos().y()]))
            self.selection_overlay.update_path(self.selection_mode, self._drag_path)
            return
        GLViewWidget.mouseMoveEvent(self.molecule_view, event)


    #------------TOOLBAR----------#
    def make_toolbar_button(self, icon_name, tooltip, callback, icon_size=26):
        btn = QPushButton()
        icon_path = os.path.join(os.path.dirname(__file__), icon_name)
        btn.setIcon(QIcon(icon_path))
        btn.setIconSize(QSize(icon_size, icon_size))
        btn.setToolTip(tooltip)
        btn.setFixedSize(34, 34)
        btn.setStyleSheet("""
            QPushButton {
                background-color: transparent;
                border: none;
                border-radius: 6px;
            }
            QPushButton:hover { background-color: #f3f4f6; }
            """)
        btn.clicked.connect(callback)
        return btn
    
    def _icon(self, name):
        return QIcon(os.path.join(os.path.dirname(__file__), name))

    # ────────────────────────────────────────────────────────────────
    # PROYECCIÓN 3D → PANTALLA
    # ────────────────────────────────────────────────────────────────

    def project_to_screen(self, pos_3d):
        w = self.molecule_view.width()
        h = self.molecule_view.height()

        proj = self.molecule_view.projectionMatrix()
        mv   = self.molecule_view.viewMatrix()
        mvp  = proj * mv  # matriz completa

        v    = QVector4D(float(pos_3d[0]), float(pos_3d[1]), float(pos_3d[2]), 1.0)
        clip = mvp.map(v)

        if clip.w() == 0:
            return None 

        ndc_x = clip.x() / clip.w()
        ndc_y = clip.y() / clip.w()

        screen_x = (ndc_x + 1.0) * 0.5 * w
        screen_y = (1.0 - ndc_y) * 0.5 * h

        return np.array([screen_x, screen_y])

    # ────────────────────────────────────────────────────────────────
    # ATOM 
    # ────────────────────────────────────────────────────────────────

    def select_atom_closest(self, pos, threshold=20, ctrl_pressed=False):
        if not self.atom_positions_3d:
            return  

        click = np.array([pos.x(), pos.y()])
        min_dist = float('inf')
        closest_idx = None

        for i, atom_pos in enumerate(self.atom_positions_3d):
            screen_pos = self.project_to_screen(atom_pos)
            if screen_pos is None:
                continue
            dist = np.linalg.norm(screen_pos - click)
            if dist < min_dist:
                min_dist = dist
                closest_idx = i

        if closest_idx is not None and min_dist < threshold:
            if ctrl_pressed:
                if closest_idx in self.selected_atoms:
                    self.selected_atoms.remove(closest_idx)
                else:
                    self.selected_atoms.append(closest_idx)
            else:
                self.selected_atoms = [closest_idx]
                print(f"Selected atom: {self.atom_types[closest_idx]} (index {closest_idx})")
            self.resaltar_atomo(closest_idx)



    def _finish_area_selection(self):
        if len(self._drag_path) < 2:
            return
        path = np.array(self._drag_path)

        if self.selection_mode == 'rect':
            x0, y0 = path[0]
            x1, y1 = path[-1]
            xmin, xmax = min(x0, x1), max(x0, x1)
            ymin, ymax = min(y0, y1), max(y0, y1)
            def inside(sx, sy):
                return xmin <= sx <= xmax and ymin <= sy <= ymax
        else:
            from matplotlib.path import Path as MplPath
            poly = MplPath(path)
            def inside(sx, sy):
                return poly.contains_point((sx, sy))

        selected = []
        for i, atom_pos in enumerate(self.atom_positions_3d):
            screen = self.project_to_screen(atom_pos)
            if screen is not None and inside(screen[0], screen[1]):
                selected.append(i)

        self.selected_atoms = selected
        for idx in selected:
            self.resaltar_atomo(idx)
        print(f"Area selection: {len(selected)} atoms")
    #────────────────────────────────────────────────────────────────
    # SHAPE FORMAS
    #────────────────────────────────────────────────────────────────
    def run_shape_calculation(self):
        if not self.selected_atoms:
            print("No atoms selected")
            return

        shape_results = {}
        cn = len(self.selected_atoms)

        if len(self.selected_atoms) == 1:
            central_idx = self.selected_atoms[0]
            shape_results, cn = calculate_shape(
                central_idx,
                self.atom_positions_3d,
                self.atom_types
            )
        else:
            coords = [self.atom_positions_3d[i] for i in self.selected_atoms]
            symbols = [self.atom_types[i] for i in self.selected_atoms]

            print(f"Calculating shape for {cn} selected atoms: {symbols}")

            try:
                geo = Geometry(symbols=symbols, positions=np.array(coords))
                shape_results = {}
                shapes = SHAPE_REFERENCES.get(cn, {})
                for shape_code, shape_name in shapes.items():
                    value = geo.get_shape_measure(shape_code, central_atom=0)
                    shape_results[shape_name] = value
                    print(f"{shape_name}: {value:.4f}")
            except Exception as e:
                print(f"Error calculating shapes for {cn} selected atoms: {e}")
                shape_results = {}
        if shape_results:
            atom_label = f"Shape_{self.atom_types[self.selected_atoms[0]]}{self.selected_atoms[0]}"
            self.export_data[atom_label] = shape_results

        if not shape_results:
            QMessageBox.information(self, "Shape Measures", f"CN={cn} — no reference shapes defined")
        else:
            atom_symbol = self.atom_types[self.selected_atoms[0]]
            atom_idx    = self.selected_atoms[0]
            self.update_shape_table(shape_results, atom_symbol, atom_idx, cn)

    def show_distortion_path(self):
        if not self.selected_atoms:
            QMessageBox.warning(self, "No selection", "Please select an atom first.")
            return
        
        table = self.shape_results_window.table
        selected_rows = list({item.row() for item in table.selectedItems()})
        if len(selected_rows) != 2:
            QMessageBox.warning(self, "Selection", 
                "Select exactly 2 shapes in the table (Ctrl+click).")
            return

        name_A = table.item(selected_rows[0], 0).text()
        name_B = table.item(selected_rows[1], 0).text()

        # Recuperar los shape codes a partir de los nombres
        if len(self.selected_atoms) == 1:
            central_idx = self.selected_atoms[0]
            central_pos = self.atom_positions_3d[central_idx]
            central_sym = self.atom_types[central_idx]
            ligand_indices = get_ligands(central_idx, self.atom_positions_3d, self.atom_types)
            coords  = [central_pos] + [self.atom_positions_3d[i] for i in ligand_indices]
            symbols = [central_sym] + [self.atom_types[i] for i in ligand_indices]
            cn = len(ligand_indices)
            central_atom_param = 1
        else:
            coords  = [self.atom_positions_3d[i] for i in self.selected_atoms]
            symbols = [self.atom_types[i] for i in self.selected_atoms]
            cn = len(self.selected_atoms)
            central_atom_param = 0

        shapes = SHAPE_REFERENCES.get(cn, {})
        code_A = next((c for c, n in shapes.items() if n == name_A), None)
        code_B = next((c for c, n in shapes.items() if n == name_B), None)

        if not code_A or not code_B:
            QMessageBox.warning(self, "Error", "No shape codes were found.")
            return

        try:
            path_sA, path_sB, mol_sA, mol_sB = calculate_distortion_path(
                symbols, coords, code_A, code_B,
                central_atom=central_atom_param
            )
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Error calculating distortion path: {e}")
            return

        if path_sA is None:
            QMessageBox.warning(self, "Error", 
                "Ideal structures could not be obtained.")
            return

        # ── Plot ────────────────────────────────────────────────────────
        dialog = QDialog(self)
        dialog.setWindowTitle(f"Distortion Path: {name_A} ↔ {name_B}")
        dialog.resize(550, 480)
        layout = QVBoxLayout(dialog)

        fig, ax = plt.subplots(figsize=(5, 4))
        fig.patch.set_facecolor('#f8f9fa')
        ax.set_facecolor('#ffffff')

        ax.plot(path_sA, path_sB,
                color='#2563eb', linewidth=2.5, label='Min. distortion path')

        # Real molecule point
        ax.scatter([mol_sA], [mol_sB],
                   color='#dc2626', s=120, zorder=5,
                   label=f'Molecule  ({mol_sA:.2f}, {mol_sB:.2f})')

        # Labels
        ax.annotate(name_A, xy=(path_sA[0], path_sB[0]),
                    xytext=(6, 6), textcoords='offset points',
                    fontsize=9, color='#1d4ed8', fontweight='bold')
        ax.annotate(name_B, xy=(path_sA[-1], path_sB[-1]),
                    xytext=(6, -14), textcoords='offset points',
                    fontsize=9, color='#1d4ed8', fontweight='bold')

        ax.set_xlabel(f'S({name_A})', fontsize=11)
        ax.set_ylabel(f'S({name_B})', fontsize=11)
        ax.set_title('Minimum Distortion Path', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        canvas = FigureCanvas(fig)
        layout.addWidget(canvas)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #b45309; color: white;
                font-weight: bold; border-radius: 4px;
                border: 2px solid black; padding: 6px;
            }""")
        layout.addWidget(close_btn)
        dialog.show()

    def run_symmetry_calculation(self):
        if not self.selected_atoms:
            QMessageBox.warning(self, "No selection", "Select an atom")
            return

        shape_label = self.shape_results_window.get_selected_shape()
        if not shape_label:
            QMessageBox.warning(self, "No shape selected", 
                            "Select one shape in the table")
            return

        try:
            if len(self.selected_atoms) == 1:
                central_idx = self.selected_atoms[0]
                central_pos = self.atom_positions_3d[central_idx]
                central_sym = self.atom_types[central_idx]

                ligand_indices = get_ligands(central_idx, self.atom_positions_3d, self.atom_types)
                coords  = [central_pos] + [self.atom_positions_3d[i] for i in ligand_indices]
                symbols = [central_sym]  + [self.atom_types[i] for i in ligand_indices]
            else:
                coords  = [self.atom_positions_3d[i] for i in self.selected_atoms]
                symbols = [self.atom_types[i] for i in self.selected_atoms]

            geo   = Geometry(symbols=symbols, positions=np.array(coords))
            cn_atoms = len(ligand_indices) if len(self.selected_atoms) == 1 else len(self.selected_atoms)
            shape_code = next(
                (code for code, name in SHAPE_REFERENCES.get(cn_atoms, {}).items() 
                if name == shape_label),
                None
                )
            symmetries = SHAPE_SYMMETRIES.get(shape_code, {})
            if not symmetries:
                QMessageBox.warning(self, "No symmetries", 
                    f"There is no symmetry for '{shape_label}'")
                return
            results = {}
            for sym_code, sym_desc in symmetries.items():
                value = geo.get_symmetry_measure(sym_code,central_atom=1)
                results[sym_code] = value
                print(f"{sym_desc}: {value:.4f}")
            
            msg = f"Symmetry measures for {shape_label}:\n\n"
            for sym_code, value in results.items():
                msg += f"  {sym_code}: {value:.4f}\n"
            QMessageBox.information(self, "Symmetry Measures", msg)

            self.export_data[f"Symmetry_{shape_label}"] = results

        except Exception as e:
            QMessageBox.warning(self, "Error", f"Error calculating symmetry: {e}")

    def update_shape_table(self, shape_results, atom_symbol, atom_idx, cn): 
        self.shape_results_window.update_results(shape_results, atom_symbol, atom_idx, cn)
    
    def toggle_polyhedron(self):
        if self._polyhedron_items:
            for item in self._polyhedron_items:
                self.molecule_view.removeItem(item)
            self._polyhedron_items = []
            self.shape_results_window.show_poly_button.setText("Show Polyhedra")
            self.shape_results_window.poly_color_button.setVisible(False)
            return

        shape_code = self.shape_results_window.get_selected_shape()
        if not shape_code or not self.selected_atoms:
            return

        if not hasattr(self,'poly_color'):
            self.poly_color = (0.1, 0.2, 0.6, 1.0)
        self.draw_polyhedron(shape_code)
        self.shape_results_window.show_poly_button.setText("Hide Polyhedra")
        self.shape_results_window.poly_color_button.setVisible(True)
        try:
            self.shape_results_window.poly_color_button.clicked.disconnect()
        except TypeError:
            pass
        self.shape_results_window.poly_color_button.clicked.connect(self.change_poly_color)

    def toggle_coord_env(self):
        if getattr(self, 'isolated_mode', False):
            self.show_all()
        else:
            self.isolate_coordination_environment()

    def open_selection_menu(self):
        menu = QMenu(self)
        menu.setStyleSheet("""
            QMenu {
                background-color: #2b2b2b;color : white;
                border: 1px solid #ccc;
                border-radius: 4px; """)
        act_rect = menu.addAction("Rectangular selection")
        act_free = menu.addAction("Free selection")
        pos = self.selection_button.mapToGlobal(self.selection_button.rect().bottomLeft())
        chosen = menu.exec_(pos)
        if chosen == act_rect:
            self.selection_mode = 'rect'
        elif chosen == act_free:
            self.selection_mode = 'free'

    def draw_polyhedron(self, shape_code):
        for item in self._polyhedron_items:
            self.molecule_view.removeItem(item)
        self._polyhedron_items = []

        central_idx = self.selected_atoms[0]
        central_pos = self.atom_positions_3d[central_idx]
        central_sym = self.atom_types[central_idx]

        ligand_indices = get_ligands(central_idx, self.atom_positions_3d, self.atom_types)
        coords  = [central_pos] + [self.atom_positions_3d[i] for i in ligand_indices]
        symbols = [central_sym]  + [self.atom_types[i] for i in ligand_indices]

        ideal = get_ideal_coords(symbols, coords, shape_code, central_atom=1)
        print('IDEAL COORDINATES')
        print(ideal)
        if ideal is None or len(ideal) == 0:
            return

        ligand_pts = ideal[1:]
        n = len(ligand_pts)

        all_dists = [
            np.linalg.norm(ligand_pts[i] - ligand_pts[j])
            for i in range(n) for j in range(i+1, n)
            ]
        min_dist = min(all_dists)
        threshold = min_dist * 1.2 

        for i in range(n):
            for j in range(i+1, n):
                dist = np.linalg.norm(ligand_pts[i] - ligand_pts[j])
                if dist < threshold:
                    line = GLLinePlotItem(
                        pos=np.array([ligand_pts[i], ligand_pts[j]]),
                        color=self.poly_color,
                        width=6,
                        antialias=True
                    )
                    self.molecule_view.addItem(line)
                    self._polyhedron_items.append(line)

    def change_poly_color(self):
        initial = QColor(
            int(self.poly_color[0] * 255),
            int(self.poly_color[1] * 255),
            int(self.poly_color[2] * 255)
        )
        qcolor = QColorDialog.getColor(initial=initial, parent=self, title="Polyhedron Color")
        if not qcolor.isValid():
            return

        self.poly_color = (
            qcolor.red()   / 255.0,
            qcolor.green() / 255.0,
            qcolor.blue()  / 255.0,
            1.0
        )
        shape_code = self.shape_results_window.get_selected_shape()
        if shape_code:
            self.draw_polyhedron(shape_code)



        #----------------------Atom selection----------------------
    def run_point_group(self):
        if not self.atom_positions_3d or not self.atom_types:
            print('No molecule')
            return
        point_group = calculate_point_group(self.atom_positions_3d,self.atom_types)
        if point_group:
            self.current_point_group = str(point_group)
            print(f'Point group: {self.current_point_group}')
            self.pointgroup_table.setRowCount(0)
            self.pointgroup_table.setRowCount(1)
            item_label = QTableWidgetItem("Point Group")
            item_label.setFlags(item_label.flags() & ~Qt.ItemIsEditable)
            item_value = QTableWidgetItem(str(point_group))
            item_value.setFlags(item_value.flags() & ~Qt.ItemIsEditable)
            item_value.setTextAlignment(Qt.AlignCenter)
            self.pointgroup_table.setItem(0, 0, item_value)
            self.pointgroup_table.setVisible(True)
            self.close_pointgroup_button.setVisible(True)
            self.export_data["Point_Group"] = {"Point Group":str(point_group) }
        else:
            self.pointgroup_table.setVisible(False)


    def resaltar_atomo(self, idx):
        for i in range(len(self.spheres)):
            props = atom_color_dict.get(self.atom_types[i], {'color': 'gray', 'size': 0.3})
            if i in self.selected_atoms:
                self.spheres[i].setColor((0.969, 0.192, 0.027, 0.35)) 
                self.spheres[i].setGLOptions('translucent')
            else:
                self.spheres[i].setColor(colors.to_rgba(props['color']))
                self.spheres[i].setGLOptions('opaque')
        if self.selected_atoms and self.atom_info_window.isVisible():
            first_idx = self.selected_atoms[0]
            sym = self.atom_types[first_idx]
            pos = self.atom_positions_3d[first_idx]
            cn = len(get_ligands(first_idx, self.atom_positions_3d, self.atom_types))

            self.atom_info_window.update_atom(
                symbol=self.atom_types[first_idx],
                index=first_idx,
                x=pos[0],
                y=pos[1],
                z=pos[2],
                cn = cn,
            )
    
    def delete_selection(self):
        if not self.selected_atoms:
            QMessageBox.information(self, "Delete", "No atoms selected.")
            return

        to_delete = set(self.selected_atoms)
        remaining = [
            (self.atom_types[i], *self.atom_positions_3d[i])
            for i in range(len(self.atom_types))
            if i not in to_delete
        ]

        self.selected_atoms = []
        self.selected_atom_index = None
        self.display_molecule(remaining)

    def restore_original(self):
        if not getattr(self, '_original_atoms', None):
            QMessageBox.information(self, 'Reset', 'Load a molecule first.')
            return
        self.selected_atoms = []
        self.display_molecule(self._original_atoms)
    def clear_symmetry_elements(self):
        if hasattr(self, '_sym_items'):
            for item in self._sym_items:
                self.molecule_view.removeItem(item)
        self._sym_items = []
        if hasattr(self, 'sym_label_overlay'):
            self.sym_label_overlay.hide()

    def reset_measurements(self):
        if hasattr(self, '_measure_items'):
            for item in self._measure_items:
                self.molecule_view.removeItem(item)
            self._measure_items = []
        self.sym_label_overlay.hide()
        self.btn_reset.setVisible(False)


    def clear_all_analysis(self):
        self.reset_measurements()
        self.clear_symmetry_elements()
        for item in self._polyhedron_items:
            self.molecule_view.removeItem(item)
        self._polyhedron_items = []
        self.selected_atoms = []
        self.selected_atom_index = None
        self.export_data = {}
        self.current_point_group = None
        self.isolated_mode = False
        self._hydrogens_hidden = False
        self.pointgroup_table.setVisible(False)
        self.close_pointgroup_button.setVisible(False)
        if self.shape_results_window.isVisible():
            self.shape_results_window.hide()
        if self.atom_info_window.isVisible():
            self.atom_info_window.hide()


    def draw_symmetry_operation(self, operation, center,label_text=None):
        self.clear_symmetry_elements()
        label = operation.label
        matrix = np.real(operation.matrix_representation)
        length = 3.0  

        if label in ('E',):
            return

        elif label == 'i':
            #inversion center
            md = MeshData.sphere(rows=20, cols=40, radius=1.15)
            sphere = GLMeshItem(meshdata=md, smooth=True, color=(0.969, 0.192, 0.027, 0.2), glOptions='translucent')
            sphere.translate(*center)
            self.molecule_view.addItem(sphere)
            self._sym_items.append(sphere)

        elif label.startswith('C') or label.startswith('S'):
            # Rotation axis
            axis = get_symmetry_axis(matrix)
            axis = axis / np.linalg.norm(axis)
            p1 = center - axis * length
            p2 = center + axis * length
            line = GLLinePlotItem(
                pos=np.array([p1, p2]),
                color=(1, 0.3, 0, 1),
                width=3,
                antialias=True
            )
            self.molecule_view.addItem(line)
            self._sym_items.append(line)
            arrow_length = 0.4
            arrow_radius = 0.12
            n_sides = 16
            tip = p2+axis*arrow_length
            base_center = p2
            
            arbitrary = np.array([1, 0, 0]) if abs(axis[0]) < 0.9 else np.array([0, 1, 0])
            u = np.cross(axis, arbitrary)
            u /= np.linalg.norm(u)
            v = np.cross(axis, u)

            angles = np.linspace(0, 2 * np.pi, n_sides, endpoint=False)
            base_verts = np.array([base_center + arrow_radius * (np.cos(a) * u + np.sin(a) * v)
                                for a in angles])

            all_verts = np.vstack([base_verts, tip])  # n_sides + 1
            tip_idx = n_sides
            faces = np.array([[i, (i + 1) % n_sides, tip_idx] for i in range(n_sides)])

            md = MeshData(vertexes=all_verts, faces=faces)
            cone = GLMeshItem(meshdata=md, smooth=True, color=(1, 0.3, 0, 1),
                          shader='shaded', glOptions='opaque')
            self.molecule_view.addItem(cone)
            self._sym_items.append(cone)
        elif label.startswith('s') or label.startswith('σ'):
            
            normal = get_plane_normal(matrix)
            normal = normal / np.linalg.norm(normal)
            self._draw_plane(center, normal, color=(0.2, 0.5, 1.0, 0.25))
            
        display = label_text if label_text else label
        self.sym_label_overlay.setText(display)
        self.sym_label_overlay.adjustSize()
        self.sym_label_overlay.show()

    def _draw_plane(self, center, normal, size=3.0, color=(0.2, 0.5, 1.0, 0.3)):
        arbitrary = np.array([1, 0, 0]) if abs(normal[0]) < 0.9 else np.array([0, 1, 0])
        v1 = np.cross(normal, arbitrary)
        v1 /= np.linalg.norm(v1)
        v2 = np.cross(normal, v1)

        corners = np.array([
            center + size * (v1+v2),
            center + size * (v1-v2),
            center + size * (-v1-v2),
            center + size * (-v1+v2)
        ])

        faces = np.array([[0, 1, 2], [0, 2, 3]])


        md = MeshData(vertexes=corners, faces=faces)
        mesh = GLMeshItem(meshdata=md, smooth=False, color=color, glOptions='translucent')
        self.molecule_view.addItem(mesh)
        self._sym_items.append(mesh)
    #---------------------------ATOM INFO ------------------------------#
    def show_atom_info(self):
        if len(self.selected_atoms) == 0:
            QMessageBox.warning(self, "No selection", "There are no atom selected")
            return
        if len(self.selected_atoms) > 1:
            QMessageBox.warning(self, "Multiple selection","More than one atom has been selected." )
            return
        
        idx = self.selected_atoms[0]
        sym = self.atom_types[idx]
        central_pos = self.atom_positions_3d[idx]

        cn = len(get_ligands(idx, self.atom_positions_3d, self.atom_types))


        self.atom_info_window.update_atom(
        symbol=sym,
        index=idx,
        x=central_pos[0],
        y=central_pos[1],
        z=central_pos[2],
        cn = cn,
    )

    #---------------------------MEASURES------------------------------#
    def toggle_measurements(self):
        if not self.atom_positions_3d:
            QMessageBox.warning(self, "No molecule", "Load a molecule first.")
            return

        menu = QMenu(self)
        menu.setStyleSheet("""
            QMenu {
                background-color: #2b2b2b;
                color: white;
                border: 1px solid black;
                font-weight: bold;
            }
            QMenu::item:selected { background-color: #555555; }
            QMenu::item { padding: 6px 20px; }
        """)

        action_dist     = menu.addAction("Distance (2 atoms)")
        action_angle    = menu.addAction("Angle (3 atoms)")
        action_dihedral = menu.addAction("Dihedral Angle (4 atoms)")

        pos = self.measure_button.mapToGlobal(self.measure_button.rect().bottomLeft())
        chosen = menu.exec_(pos)

        if chosen == action_dist:
            self.measure_distance()
        elif chosen == action_angle:
            self.measure_angle()
        elif chosen == action_dihedral:
            self.measure_dihedral()

    def measure_distance(self):
        if len(self.selected_atoms) != 2:
            QMessageBox.warning(self, "Selection", "Select exactly 2 atoms.")
            return

        idx1, idx2 = self.selected_atoms[0], self.selected_atoms[1]
        p1 = self.atom_positions_3d[idx1]
        p2 = self.atom_positions_3d[idx2]
        d  = calc_distance(p1, p2)
        s1 = f"{self.atom_types[idx1]}{idx1}"
        s2 = f"{self.atom_types[idx2]}{idx2}"

        if hasattr(self, '_measure_items'):
            for item in self._measure_items:
                self.molecule_view.removeItem(item)
        self._measure_items = []

        # Línea de puntos
        dash = create_dashed_line(p1, p2, color=(0.969, 0.192, 0.027, 0.65))
        self.molecule_view.addItem(dash)
        self._measure_items.append(dash)

        # Label overlay
        self.sym_label_overlay.setText(f"{s1} — {s2}:  {d:.3f} Å")
        self.sym_label_overlay.adjustSize()
        self.sym_label_overlay.show()
        self.btn_reset.setVisible(True)

    def measure_angle(self):
        if len(self.selected_atoms) != 3:
            QMessageBox.warning(self, "Selection", "Select exactly 3 atoms.\nAngle is calculated at the middle atom.")
            return

        idx1, idx2, idx3 = self.selected_atoms[0], self.selected_atoms[1], self.selected_atoms[2]
        p1 = self.atom_positions_3d[idx1]
        p2 = self.atom_positions_3d[idx2] 
        p3 = self.atom_positions_3d[idx3]
        a  = calc_angle(p1, p2, p3)
        s1 = f"{self.atom_types[idx1]}{idx1}"
        s2 = f"{self.atom_types[idx2]}{idx2}"
        s3 = f"{self.atom_types[idx3]}{idx3}"

        
        if hasattr(self, '_measure_items'):
            for item in self._measure_items:
                self.molecule_view.removeItem(item)
        self._measure_items = []

        
        dash1 = create_dashed_line(p1, p2, color=(0.969, 0.192, 0.027, 0.65))
        dash2 = create_dashed_line(p3, p2, color=(0.969, 0.192, 0.027, 0.65))
        self.molecule_view.addItem(dash1)
        self.molecule_view.addItem(dash2)
        self._measure_items.extend([dash1, dash2])

        
        v1 = (p1 - p2) / np.linalg.norm(p1 - p2)
        v2 = (p3 - p2) / np.linalg.norm(p3 - p2)
        arc_radius = min(np.linalg.norm(p1 - p2), np.linalg.norm(p3 - p2)) * 0.3
        n_arc = 30
        arc_points = []
        for i in range(n_arc + 1):
            t = i / n_arc
            vec = (1 - t) * v1 + t * v2
            vec = vec / np.linalg.norm(vec)
            arc_points.append(p2 + arc_radius * vec)

        arc_line = GLLinePlotItem(
            pos=np.array(arc_points),
            color=(0.969, 0.192, 0.027, 0.65),
            width=5,
            antialias=True
        )
        self.molecule_view.addItem(arc_line)
        self._measure_items.append(arc_line)

        self.sym_label_overlay.setText(f"{s1} — {s2} — {s3}:  {a:.2f}°")
        self.sym_label_overlay.adjustSize()
        self.sym_label_overlay.show()
        self.btn_reset.setVisible(True)


    def measure_dihedral(self):
        if len(self.selected_atoms) != 4:
            QMessageBox.warning(self, "Selection", "Select exactly 4 atoms.")
            return

        idx1, idx2, idx3, idx4 = self.selected_atoms[0], self.selected_atoms[1], self.selected_atoms[2], self.selected_atoms[3]
        p1 = self.atom_positions_3d[idx1]
        p2 = self.atom_positions_3d[idx2]
        p3 = self.atom_positions_3d[idx3]
        p4 = self.atom_positions_3d[idx4]
        d = calc_dihedral(p1, p2, p3, p4)

        s1 = f"{self.atom_types[idx1]}{idx1}"
        s2 = f"{self.atom_types[idx2]}{idx2}"
        s3 = f"{self.atom_types[idx3]}{idx3}"
        s4 = f"{self.atom_types[idx4]}{idx4}"

        if hasattr(self, '_measure_items'):
            for item in self._measure_items:
                self.molecule_view.removeItem(item)
        self._measure_items = []

        for a, b in [(p1, p2), (p2, p3), (p3, p4)]:
            dash = create_dashed_line(a, b, color=(0.122, 0.467, 0.957, 0.65))
            self.molecule_view.addItem(dash)
            self._measure_items.append(dash)

        axis = p3 - p2
        axis_norm = np.linalg.norm(axis)

        if axis_norm > 1e-6:
            axis_u = axis / axis_norm

            u = (p1 - p2)
            u = u - np.dot(u, axis_u) * axis_u
            u_norm = np.linalg.norm(u)

            v = (p4 - p3)
            v = v - np.dot(v, axis_u) * axis_u
            v_norm = np.linalg.norm(v)

            if u_norm > 1e-6 and v_norm > 1e-6:
                u = u / u_norm
                v = v / v_norm

                arc_radius = min(np.linalg.norm(p1 - p2), np.linalg.norm(p3 - p2)) * 0.35
                arc_center = p2

                n_arc = 40
                arc_points = []
                for i in range(n_arc + 1):
                    t = i / n_arc
                    vec = (1 - t) * u + t * v
                    vec_norm = np.linalg.norm(vec)
                    if vec_norm > 1e-6:
                        vec = vec / vec_norm
                    arc_points.append(arc_center + arc_radius * vec)

                arc_line = GLLinePlotItem(
                    pos=np.array(arc_points),
                    color=(0.122, 0.467, 0.957, 0.65),
                    width=5,
                    antialias=True
                )
                self.molecule_view.addItem(arc_line)
                self._measure_items.append(arc_line)

        self.sym_label_overlay.setText(f"{s1} — {s2} — {s3} — {s4}:  {d:.2f}°")
        self.sym_label_overlay.adjustSize()
        self.sym_label_overlay.show()
        self.btn_reset.setVisible(True)


    #---------------------------ISOLATE------------------------------#
    def toggle_isolate(self):
        if not hasattr(self, 'isolated_mode'):
            self.isolated_mode = False


        if self.isolated_mode:
            self.show_all()
            return


        menu = QMenu(self)
        menu.setStyleSheet("""
            QMenu {
                background-color: #2b2b2b;
                color: white;
                border: 1px solid black;
                font-weight: bold;
            }
            QMenu::item:selected {
                background-color: #555555;
            }
            QMenu::item {
                padding: 6px 20px;
            }
        """)

        action_coord = menu.addAction("Coordination Environment")
        action_sel   = menu.addAction("Selection")


        pos = self.isolate_button.mapToGlobal(
            self.isolate_button.rect().bottomLeft()
        )
        chosen = menu.exec_(pos)

        if chosen == action_coord:
            self.isolate_coordination_environment()
        elif chosen == action_sel:
            self.isolate_selection()

    def isolate_atoms(self, visible_indices):
        visible_set = set(visible_indices)
        for i, sphere in enumerate(self.spheres):
            sphere.setVisible(i in visible_set)
        for bond_data, bond_item in zip(self._bond_pairs, self.bonds):
            i, j = bond_data
            bond_item.setVisible(i in visible_set and j in visible_set)
        self.isolated_mode = True

    def isolate_selection(self):
        if not self.selected_atoms:
            print("No atoms selected")
            return
        self.isolate_atoms(self.selected_atoms)
        self.isolate_button.setText("Show All")

    def isolate_coordination_environment(self):
        if not self.selected_atoms:
            print("No atoms selected")
            return

        central_idx = self.selected_atoms[0]
        central_pos = self.atom_positions_3d[central_idx]
        central_sym = self.atom_types[central_idx]

    
        ligand_indices = get_ligands(central_idx, self.atom_positions_3d, self.atom_types)

        visible_indices = [central_idx] + ligand_indices
        self.isolate_atoms(visible_indices)
        

    def show_all(self):
        for sphere in self.spheres:
            sphere.setVisible(True)
        for bond in self.bonds:
            bond.setVisible(True)
        self.isolated_mode = False
        

    def toggle_hydrogens(self):
        if not hasattr(self, 'spheres') or not self.spheres:
            return

        self._hydrogens_hidden = not getattr(self, '_hydrogens_hidden', False)

        for i, (sphere, atom_type) in enumerate(zip(self.spheres, self.atom_types)):
            if atom_type == 'H':
                sphere.setVisible(not self._hydrogens_hidden)

    
        for bond_data, bond_item in zip(self._bond_pairs, self.bonds):
            i, j = bond_data
            if self.atom_types[i] == 'H' or self.atom_types[j] == 'H':
                bond_item.setVisible(not self._hydrogens_hidden)

        self.hide_h_button.setToolTip(
            "Show Hydrogens" if self._hydrogens_hidden else "Hide Hydrogens"
        )


    #-----------------------XYZ AXIS------------------------ 

    def crear_ejes_xyz(self):
        eje_x = GLLinePlotItem(pos=np.array([[3, -3, -3], [4, -3, -3]]), color=(1, 0, 0, 1), width=3, antialias=True)
        eje_y = GLLinePlotItem(pos=np.array([[3, -3, -3], [3, -2, -3]]), color=(0, 1, 0, 1), width=3, antialias=True)
        eje_z = GLLinePlotItem(pos=np.array([[3, -3, -3], [3, -3, -2]]), color=(0, 0, 1, 1), width=3, antialias=True)
        return eje_x, eje_y, eje_z

    def toggle_ejes_xyz(self):
        if self.ejes_visibles:
            for eje in [self.eje_x, self.eje_y, self.eje_z]:
                if eje is not None and eje in self.molecule_view.items:
                    self.molecule_view.removeItem(eje)
            self.eje_x = self.eje_y = self.eje_z = None
            self.ejes_visibles = False
        else:
            self.eje_x, self.eje_y, self.eje_z = self.crear_ejes_xyz()
            for eje in [self.eje_x, self.eje_y, self.eje_z]:
                self.molecule_view.addItem(eje)
            self.ejes_visibles = True

    #----------------------Load and visualization--------------------
    def display_molecule(self, atoms):
        if hasattr(self, 'spheres'):
            for s in self.spheres:
                self.molecule_view.removeItem(s)
        if hasattr(self, 'bonds'):
            for b in self.bonds:
                self.molecule_view.removeItem(b)

        self.spheres = []
        self.bonds = []
        self._bond_pairs = []   
        self.atom_positions_3d = []
        self.atom_types = []
        self.selected_atom_index = None

        for atom_type, x, y, z in atoms:
            props = atom_color_dict.get(atom_type, {'color': 'gray', 'size': 0.15})
            color = colors.to_rgba(props['color'])
            sphere = create_sphere(radius=props['size']*2, color=color)
            sphere.translate(x, y, z)
            self.molecule_view.addItem(sphere)
            self.spheres.append(sphere)
            self.atom_positions_3d.append(np.array([x, y, z]))
            self.atom_types.append(atom_type)

        bond_count = 0
        for i, (type1, x1, y1, z1) in enumerate(atoms):
            for j, (type2, x2, y2, z2) in enumerate(atoms):
                if i < j:
                    pos1 = np.array([x1, y1, z1])
                    pos2 = np.array([x2, y2, z2])
                    dist = np.linalg.norm(pos1 - pos2)
                    threshold = covalent_threshold(type1, type2)
                    
                    if dist < threshold:
                        bond_count += 1
                        cylinder = create_cylinder(pos1, pos2, radius=0.06, color=(0.3, 0.3, 0.3, 1))
                        self._bond_pairs.append((i, j))  
                        self.bonds.append(cylinder)      
                        self.molecule_view.addItem(cylinder)
        if self.atom_positions_3d:
            centroid = np.mean(self.atom_positions_3d,axis=0)
            self.molecule_view.opts['center'] = pg.Vector(centroid[0], centroid[1], centroid[2])
            max_extent = max(np.linalg.norm(pos-centroid) for pos in self.atom_positions_3d)
            self.molecule_view.opts['distance'] = max(max_extent * 3,10)

        

    def load_file(self):
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Molecule File", "",
            "Structure Files (*.xyz *.pdb *.cif);;All Files (*)",
            options=options
        )
        if not file_path:
            return
        try:
            atoms,cell = load_structure(file_path)
        except Exception as e:
            QMessageBox.warning(self, "Load error", f"Could not read file:\n{e}")
            return
        self.cell = cell
        self._original_atoms = list(atoms)
        self.clear_all_analysis()
        self.display_molecule(atoms)
        self.btn_reset.setVisible(False)


    def close_application(self):
        self.close()

    def export_image(self, fmt="png"):
        from PyQt5.QtWidgets import QComboBox

        fmt_labels = {
            "png":  "PNG Files (*.png)",
            "tiff": "TIFF Files (*.tiff)",
            "jpg":  "JPEG Files (*.jpg)",
        }

        dialog = QDialog(self)
        dialog.setWindowTitle("Export Image")
        dialog.setFixedSize(300, 200)
        layout = QVBoxLayout(dialog)

        layout.addWidget(QLabel("Resolution:"))
        combo = QComboBox()
        combo.addItems([
            "1x  (screen size)",
            "2x  (publication)",
            "4x  (high resolution)",
        ])
        combo.setCurrentIndex(1)
        combo.setStyleSheet("padding: 4px; font-size: 13px;")
        layout.addWidget(combo)

        layout.addWidget(QLabel("Background:"))
        bg_combo = QComboBox()
        bg_combo.addItems(["Current (grey)", "White", "Black"])
        bg_combo.setStyleSheet("padding: 4px; font-size: 13px;")
        layout.addWidget(bg_combo)

        btn_row = QHBoxLayout()
        ok_btn = QPushButton("Export")
        ok_btn.setStyleSheet("""QPushButton {
            background-color: #5f6064; color: white;
            font-weight: bold; border-radius: 4px;
            border: 2px solid black; padding: 10px;}""")
        cancel_btn = QPushButton("Cancel")
        cancel_btn.setStyleSheet(ok_btn.styleSheet())
        ok_btn.clicked.connect(dialog.accept)
        cancel_btn.clicked.connect(dialog.reject)
        btn_row.addWidget(ok_btn)
        btn_row.addWidget(cancel_btn)
        layout.addLayout(btn_row)

        if dialog.exec_() != QDialog.Accepted:
            return

        factor = [1, 2, 4][combo.currentIndex()]
        bg_choice = bg_combo.currentIndex()

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Image", f"molecule.{fmt}",
            f"{fmt_labels[fmt]};;All Files (*)"
        )
        if not file_path:
            return

        if bg_choice == 1:
            self.molecule_view.setBackgroundColor('w')
        elif bg_choice == 2:
            self.molecule_view.setBackgroundColor('k')

        original_size = self.molecule_view.size()
        w = original_size.width()  * factor
        h = original_size.height() * factor
        self.molecule_view.resize(w, h)
        self.molecule_view.repaint()


        img = self.molecule_view.grabFramebuffer()

        self.molecule_view.resize(original_size)
        if bg_choice in (1, 2):
            self.molecule_view.setBackgroundColor('#f8f9fa')
        self.molecule_view.repaint()

        quality = 95 if fmt == "jpg" else -1
        if img.save(file_path, fmt.upper(), quality):
            print(f"Image saved: {file_path}  ({img.width()}x{img.height()} px)")
        else:
            QMessageBox.warning(self, "Error", f"Could not save image to:\n{file_path}")


    def export_f(self):
        if not self.export_data:
            QMessageBox.information(self, "Export", "No results to export yet.")
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Results", "results.csv", "CSV Files (*.csv)"
        )
        if not file_path:
            return
        try:
            with open(file_path, "w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh, delimiter=';')
                for section, entries in self.export_data.items():
                    writer.writerow([section])
                    for prop, value in entries.items():
                        writer.writerow([prop, value])
                    writer.writerow([])
            print(f"Results exported: {file_path}")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not export results:\n{e}")

    def toggle_export(self):
        menu = QMenu(self)
        menu.setStyleSheet("""
            QMenu {
                background-color: #2b2b2b;
                color: white;
                border: 1px solid black;
                font-weight: bold;
            }
            QMenu::item:selected { background-color: #555555; }
            QMenu::item { padding: 10px 20px; }
        """)

        action_results = menu.addAction("Export Results (.csv)")
        menu.addSeparator()
        action_png  = menu.addAction("Export Image as PNG")
        action_tiff = menu.addAction("Export Image as TIFF")
        action_jpg  = menu.addAction("Export Image as JPEG")

        pos = self.export_button.mapToGlobal(self.export_button.rect().bottomLeft())
        chosen = menu.exec_(pos)

        if chosen == action_results:
            self.export_f()
        elif chosen == action_png:
            self.export_image(fmt="png")
        elif chosen == action_tiff:
            self.export_image(fmt="tiff")
        elif chosen == action_jpg:
            self.export_image(fmt="jpg")

    def run_symmetry_posym(self):
        if not self.selected_atoms:
            QMessageBox.warning(self, "No selection", "Please select an atom first.")
            return
        if not hasattr(self, 'current_point_group') or not self.current_point_group:
            QMessageBox.warning(self, "Missing step", "Ep! First get the point group!")
            return
        try:
            if len(self.selected_atoms) == 1:
                central_idx = self.selected_atoms[0]
                central_pos = self.atom_positions_3d[central_idx]
                central_sym = self.atom_types[central_idx]

                ligand_indices = get_ligands(central_idx, self.atom_positions_3d, self.atom_types)
                coords  = [central_pos] + [self.atom_positions_3d[i] for i in ligand_indices]
                symbols = [central_sym]  + [self.atom_types[i] for i in ligand_indices]
            else:
                coords  = [self.atom_positions_3d[i] for i in self.selected_atoms]
                symbols = [self.atom_types[i] for i in self.selected_atoms]

            if not self.current_point_group:
                QMessageBox.warning(self, "No point group", "Could not determine the molecule's point group.")
                return

            results, operations,center = calculate_symmetry_posym(coords, symbols, groups=[self.current_point_group])

            if not results:
                QMessageBox.information(self, "Error",'No compatible point group found')
                return

            self.show_posym_results(results,operations,center)

            central_label = f"{self.atom_types[self.selected_atoms[0]]}{self.selected_atoms[0]}"
            self.export_data[f"CSM_{central_label}"] = results

            if operations:
                ops_data = {}
                counter = {}
                for op in operations:
                    counter[op.label] = counter.get(op.label, 0) + 1
                    n_total = sum(1 for o in operations if o.label == op.label)
                    name = f"{op.label}_{counter[op.label]}" if n_total > 1 else op.label
                    matrix = np.round(np.real(op.matrix_representation), 4)
                    ops_data[name] = " | ".join(
                        "[" + ", ".join(f"{v:.4f}" for v in row) + "]" for row in matrix
                    )
                self.export_data[f"SymmetryOps_{central_label}"] = ops_data

        except Exception as e:
            QMessageBox.warning(self, "Error posym", f"Error calculating CSM: {e}")

    def show_posym_results(self, results, operations,center):
        dialog = QDialog(self)
        dialog.setWindowTitle("Continuous Symmetry Measures")
        dialog.resize(380, 550)

        layout = QVBoxLayout()
        dialog.setLayout(layout)

        label = QLabel("CSM (0 = perfectly symmetric)")
        label.setStyleSheet("font-size: 11px; color: #555;")
        layout.addWidget(label)

        table = QTableWidget()
        table.setColumnCount(2)
        table.setFixedHeight(80)
        table.setHorizontalHeaderLabels(["Point Group", "CSM"])
        table.setRowCount(len(results))
        table.verticalHeader().setVisible(False)
        table.setShowGrid(False)
        table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        table.setStyleSheet("""
            QTableWidget { background-color: #ffffff; border: 1px solid #e5e7eb;
                            border-radius: 8px; font-size: 13px; }
            QTableWidget::item { padding: 6px 12px; border-bottom: 1px solid #f3f4f6; color: #1f2937; }
            QTableWidget::item:selected { background-color: #f0f4ff; color: #1f2937; }
            QHeaderView::section { background-color: #f9fafb; color: #6b7280; font-weight: 600;
                                    font-size: 12px; padding: 6px 12px; border-bottom: 1px solid #e5e7eb; }
        """)

        for row, (group, value) in enumerate(results.items()):
            item_group = QTableWidgetItem(group)
            item_value = QTableWidgetItem(f"{value:.4f}")
            item_group.setFlags(item_group.flags() & ~Qt.ItemIsEditable)
            item_value.setFlags(item_value.flags() & ~Qt.ItemIsEditable)
            item_value.setTextAlignment(Qt.AlignCenter)
            if value < 0.5:
                item_group.setForeground(QColor("#16a34a"))
                item_value.setForeground(QColor("#16a34a"))
            elif value < 3.0:
                item_group.setForeground(QColor("#d97706"))
                item_value.setForeground(QColor("#d97706"))
            table.setItem(row, 0, item_group)
            table.setItem(row, 1, item_value)

        layout.addWidget(table)

        # ── Operaciones de simetría ──────────────────────────────────
        if operations:
            ops_label = QLabel("Symmetry operations")
            ops_label.setStyleSheet("font-size: 12px; font-weight: bold; color: #1f2937; margin-top: 8px;")
            layout.addWidget(ops_label)

            # Agrupar operaciones por label
            from collections import defaultdict
            grouped = defaultdict(list)
            for op in operations:
                grouped[op.label].append(op)

            ops_table = QTableWidget()
            ops_table.setColumnCount(2)
            ops_table.setHorizontalHeaderLabels(["Operation", "Matrix"])
            ops_table.verticalHeader().setVisible(False)
            ops_table.setShowGrid(False)
            ops_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
            ops_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
            ops_table.setStyleSheet("""
                QTableWidget { background-color: #ffffff; border: 1px solid #e5e7eb;
                                border-radius: 8px; font-size: 12px; }
                QTableWidget::item { padding: 4px 10px; border-bottom: 1px solid #f3f4f6; color: #1f2937; }
                QTableWidget::item:selected { background-color: #f0f4ff; color: #1f2937; }
                QHeaderView::section { background-color: #f9fafb; color: #6b7280; font-weight: 600;
                                        font-size: 12px; padding: 6px 10px; border-bottom: 1px solid #e5e7eb; }
            """)

            rows = []  # lista de (label, matrix_str, is_header, parent_label, nth)
            for label, ops in grouped.items():
                if len(ops) == 1:
                    matrix_str = "; ".join(
                        "[" + ", ".join(f"{v:.2f}" for v in r) + "]"
                        for r in np.round(ops[0].matrix_representation, decimals=2)
                    )
                    rows.append((label, matrix_str, False, None, 0))
                else:
                    rows.append((f"▶  {label}  ({len(ops)})", "", True, label, 0))
                    for n, op in enumerate(ops):
                        matrix_str = "; ".join(
                            "[" + ", ".join(f"{v:.2f}" for v in r) + "]"
                            for r in np.round(op.matrix_representation, decimals=2)
                        )
                        rows.append((f"   {label}_{n+1}", matrix_str, False, label, n+1))

            ops_table.setRowCount(len(rows))

            hidden_rows = {}  
            expanded = {}     

            for row_idx, (lbl, matrix_str, is_header, parent, nth) in enumerate(rows):
                item_label = QTableWidgetItem(lbl)
                item_matrix = QTableWidgetItem(matrix_str)
                item_label.setFlags(item_label.flags() & ~Qt.ItemIsEditable)
                item_matrix.setFlags(item_matrix.flags() & ~Qt.ItemIsEditable)

                if is_header:
                    item_label.setForeground(QColor("#0354bb"))
                    item_label.setFont(item_label.font())
                    expanded[parent] = False
                    hidden_rows[parent] = []
                elif parent and nth > 0:
                    hidden_rows[parent].append(row_idx)
                    ops_table.setRowHidden(row_idx, True)  

                ops_table.setItem(row_idx, 0, item_label)
                ops_table.setItem(row_idx, 1, item_matrix)

            self._current_operations = operations
            self._current_sym_center = center

            def on_cell_clicked(row, col, rows=rows, hidden_rows=hidden_rows, expanded=expanded):
                lbl, _, is_header, parent, nth = rows[row]

                if is_header:
                    label_key = parent
                    expanded[label_key] = not expanded[label_key]
                    for child_row in hidden_rows[label_key]:
                        ops_table.setRowHidden(child_row, not expanded[label_key])
                    arrow = "▼" if expanded[label_key] else "▶"
                    count = len(hidden_rows[label_key])
                    ops_table.item(row, 0).setText(f"{arrow}  {label_key}  ({count})")

                else:
                    label_clean = lbl.strip()
                    count = {}
                    for op in self._current_operations:
                        count[op.label] = count.get(op.label, 0) + 1
                        n_total = len([o for o in self._current_operations if o.label == op.label])
                        suffix = f"_{count[op.label]}" if n_total > 1 else ""
                        full = f"{op.label}{suffix}"
                        if full == label_clean or op.label == label_clean:
                            self.draw_symmetry_operation(op, self._current_sym_center, label_text=label_clean)
                            break

            ops_table.cellClicked.connect(on_cell_clicked)
            layout.addWidget(ops_table)
        self.clear_sym_button = QPushButton("Clear Symmetry Elements")
        self.clear_sym_button.setFixedHeight(30)
        self.clear_sym_button.setStyleSheet("""
        QPushButton {
            background-color: #74034e;
            color: white;
            font-weight: bold;
            border-radius: 4px;
            border: 2px solid black; }
            QPushButton:hover { background-color: #74034e }
        """)
        self.clear_sym_button.clicked.connect(self.clear_symmetry_elements)
        layout.addWidget(self.clear_sym_button)
        dialog.show()


    def show_info(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("About Us")
        dialog.setFixedSize(420, 320)
        dialog.setStyleSheet("background-color: #1e1e2e;")
        layout = QVBoxLayout(dialog)
        layout.setSpacing(10)
        layout.setContentsMargins(24, 20, 24, 20)
        #-------LOGO------
        logo_label = QLabel()
        logo_path = os.path.join(os.path.dirname(__file__), "molektra.png")  
        if os.path.exists(logo_path):
            pix = QPixmap(logo_path).scaled(110, 110, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            logo_label.setPixmap(pix)
        else:
            logo_label.setText("X")
            logo_label.setStyleSheet("font-size: 48px;")
        logo_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(logo_label)

        #---NAME---
        title = QLabel("Molektra")
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("color: white; font-size: 16px; font-weight: bold;")
        layout.addWidget(title)

        subtitle = QLabel("A visual tool for shape and symmetry analysis of molecules")
        subtitle.setAlignment(Qt.AlignCenter)
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color: #9ca3af; font-size: 11px;")
        layout.addWidget(subtitle)

        sep = QLabel()
        sep.setFixedHeight(1)
        sep.setStyleSheet("background-color: #374151;")
        layout.addWidget(sep)

    #------------LINKS------------------#
        links = [
            ("GitHub",  "https://github.com/",  "#6366f1"),
            ("Zenodo",  "https://zenodo.org/",        "#0ea5e9"),
            ("Paper",   "https://doi.org/",                 "#10b981"),
        ]

        btn_row = QHBoxLayout()
        btn_row.setSpacing(10)
        for label, url, color in links:
            btn = QPushButton(label)
            btn.setFixedHeight(34)
            btn.setStyleSheet(f"""
                QPushButton {{
                    background-color: {color};
                    color: white;
                    font-weight: bold;
                    border-radius: 6px;
                    border: none;
                    font-size: 12px;
                }}
                QPushButton:hover {{ opacity: 0.85; }}
            """)
            btn.setCursor(Qt.PointingHandCursor)
            btn.clicked.connect(lambda _, u=url: QDesktopServices.openUrl(QUrl(u)))
            btn_row.addWidget(btn)
        layout.addLayout(btn_row)

        layout.addStretch()

        # ── Créditos + Close ─────────────────────────────────────────────
        credits = QLabel("© 2026 Electronic Structure & Symmetry Group · Version 1.0")
        credits.setAlignment(Qt.AlignCenter)
        credits.setStyleSheet("color: #6b7280; font-size: 10px;")
        layout.addWidget(credits)

        close_btn = QPushButton("Close")
        close_btn.setFixedHeight(32)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #374151; color: white;
                font-weight: bold; border-radius: 4px;
                border: 2px solid #4b5563;
            }
            QPushButton:hover { background-color: #4b5563; }
        """)
        close_btn.clicked.connect(dialog.close)
        layout.addWidget(close_btn)

        dialog.exec_()