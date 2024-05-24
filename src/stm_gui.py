"""
This module contains functions to build the truss/strut-and-tie model and to plot the corresponding stress fields.
Version 0.3: Now contains a simple GUI.
----------------
Older versions:
Version 0.1: Initial release, only includes validation of strut-and-tie model without nodal zones.
Version 0.2: Includes validation of hydrostatic nodes for corresponding stress fields with concentrated struts and ties.
"""

# Copyright (C) 2024
# Karin Yu
# ETH Zurich

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

__author__ = 'Karin Yu'
__email__ = 'karin.yu@ibk.baug.ethz.ch'
__copyright__ = 'Copyright 2024, Karin Yu'
__license__ = 'Apache 2.0'
__version__ = '0.3'
__maintainer__ = 'Karin Yu'

# Tutorial
from PyQt6.QtGui import *
from PyQt6.QtWidgets import *
from PyQt6.QtCore import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import matplotlib
import math

import sys
import stm_trusssystem as TS
from stm_methods import *

### TODOs, Karin Yu, 24.05.2024
# TODO: Add remove node
# TODO: Add remove edge
# TODO: Add modify node
# TODO: Add modify edge connection
# TODO: Save model or combine with jupyter notebook
# TODO: Add functions to modify material properties

class MainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setWindowTitle("StrutandTieModelling")

        self.resize(QSize(800, 500))

        # Central Widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        layout_all = QHBoxLayout(central_widget)  # layout

        layoutButtons = QVBoxLayout()  # layout
        self.layoutDraw = QVBoxLayout()  # layout

        # Variables specific to STM
        self.Truss = TS.TrussSystem()
        self.polygon_in = []
        self.polygon_out = []

        # Default material
        # Reinforcement material values
        Es = 205000  # MPa
        fsy = 435  # MPa
        eps_sy = fsy / Es
        ds = 1
        # Concrete material
        Ec = 30000  # MPa
        fck = 20  # MPa, fcd
        eps_c = 0.003
        dc = -1
        self.Truss.addMaterial(TS.Material(Ec,fck,eps_c,dc))
        self.Truss.addMaterial(TS.Material(Es,fsy,eps_sy,ds))

        # Label with material properties
        self.text_material = QLabel(self)
        curr_text = "Material properties:\nConcrete: Ec = {:.0f} GPa, fck = ".format(Ec/1e3) + str(fck) + " MPa\nSteel: Es = {:.0f} GPa, fsy = ".format(Es/1e3) + str(fsy) + " MPa"
        self.text_material.setText(curr_text)
        layoutButtons.addWidget(self.text_material, stretch = 1)

        # Button to add thickness
        self.button_add_thickness = QPushButton("Add thickness")
        self.button_add_thickness.clicked.connect(self.add_thickness)
        layoutButtons.addWidget(self.button_add_thickness, stretch = 3)

        # Button to add polygon in
        self.button_add_polygon_in = QPushButton("Add polygon")
        self.button_add_polygon_in.clicked.connect(lambda: self.add_polygon_clicked(True))
        layoutButtons.addWidget(self.button_add_polygon_in, stretch = 3)

        # Button to add polygon out
        """self.button_add_polygon_out = QPushButton("Add polygon outside")
        self.button_add_polygon_out.clicked.connect(lambda: self.add_polygon_clicked(False))
        layoutButtons.addWidget(self.button_add_polygon_out, stretch = 3)"""

        # Button to add nodes
        self.button_add_node = QPushButton("Add node")
        self.button_add_node.clicked.connect(self.add_node)
        layoutButtons.addWidget(self.button_add_node, stretch = 3)

        # Button to add external force
        self.button_add_support = QPushButton("Add support")
        self.button_add_support.clicked.connect(self.add_support)
        layoutButtons.addWidget(self.button_add_support, stretch=3)

        # Button to add external force
        self.button_add_force = QPushButton("Add force")
        self.button_add_force.clicked.connect(self.add_force)
        layoutButtons.addWidget(self.button_add_force, stretch=3)

        # Button to add edges
        self.button_add_edges = QPushButton("Add edge")
        self.button_add_edges.clicked.connect(self.add_edge)
        layoutButtons.addWidget(self.button_add_edges, stretch = 3)

        # Button to draw truss
        self.button_draw_truss = QPushButton("Draw truss")
        self.button_draw_truss.clicked.connect(self.draw_truss)
        layoutButtons.addWidget(self.button_draw_truss, stretch = 3)

        # Button to calculate truss
        self.button_solve_truss = QPushButton("Solve truss")
        self.button_solve_truss.clicked.connect(self.solve_truss)
        layoutButtons.addWidget(self.button_solve_truss, stretch=3)

        # Button to calculate truss
        self.button_stress_fields = QPushButton("Show stress fields")
        self.button_stress_fields.clicked.connect(self.stress_fields)
        layoutButtons.addWidget(self.button_stress_fields, stretch=3)

        # Button to clear all
        self.button_clear = QPushButton("Clear all")
        self.button_clear.clicked.connect(self.clear_all)
        layoutButtons.addWidget(self.button_clear, stretch = 3)

        # Label with copyright
        self.text_copyright = QLabel(self)
        curr_text = "Apache License, 2024, Karin Yu, ETH Zurich"
        self.text_copyright.setText(curr_text)
        layoutButtons.addWidget(self.text_copyright, stretch=1)

        self.canvas = pltCanvas(plt.figure())
        self.layoutDraw.addWidget(self.canvas)

        layout_all.addLayout(layoutButtons)
        layout_all.addLayout(self.layoutDraw)

    def solve_truss(self):
        for e in self.Truss.EdgeList:
            e.force = 0
            e.area = 1
            e.stress = 0
        self.Truss.solveTruss(update=True)
        if not self.Truss.checkEquilibrium():
            self.display_message("Truss is not in equilibrium", "Error", QMessageBox.Icon.Warning)
        else:
            self.draw_truss(wo_forces=False)

    def add_thickness(self):
        dial_thickness = ThicknessDialog(self)
        dial_thickness.exec()
        if dial_thickness.getThickness() <= 0:
            self.display_message("Thickness must be positive", "Error", QMessageBox.Icon.Warning)
        else:
            self.thickness = dial_thickness.getThickness()

    def display_message(self, message, title_message="Message", icon_type = QMessageBox.Icon.Information):
        app.setAttribute(Qt.ApplicationAttribute.AA_DontUseNativeDialogs, True)
        msg = QMessageBox()
        msg.setWindowTitle(title_message)
        msg.setIcon(icon_type)
        msg.setText(message)
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()
        app.setAttribute(Qt.ApplicationAttribute.AA_DontUseNativeDialogs, False)

    def on_click_return_node(self, event):
        x = event.xdata
        y = event.ydata
        curr_point = TS.Point(x, y, 0)
        print('Point', curr_point.printPoint())
        for n in self.Truss.NodeList:
            if n.Point == curr_point:
                self.curr_node = n
        self.curr_node = None

    def clear_all(self):
        # Variables specific to STM
        self.Truss = TS.TrussSystem()
        self.polygon_in = []
        self.polygon_out = []

        # Default material
        # Reinforcement material values
        Es = 205000  # MPa
        fsy = 435  # MPa
        eps_sy = fsy / Es
        ds = 1
        # Concrete material
        Ec = 30000  # MPa
        fck = 20  # MPa, fcd
        eps_c = 0.003
        dc = -1
        self.Truss.addMaterial(TS.Material(Ec, fck, eps_c, dc))
        self.Truss.addMaterial(TS.Material(Es, fsy, eps_sy, ds))

        # clear canvas
        fig = plt.figure(figsize=(10, 5))
        canvas = pltCanvas(fig)
        canvas.draw()
        self.canvas.figure.clf()  # Clear the figure
        self.layoutDraw.replaceWidget(self.canvas, canvas)
        self.canvas = canvas

    def stress_fields(self):
        fig = self.Truss.plotStressField(fig_size=(10, 5), t=self.thickness, tol = 0.5)
        canvas = pltCanvas(fig)
        canvas.draw()
        self.canvas.figure.clf()  # Clear the figure
        self.layoutDraw.replaceWidget(self.canvas, canvas)
        self.canvas = canvas

    def draw_truss(self,  wo_forces= True):
        fig = self.Truss.plotSTM(fig_size=(10, 5), polygon_in=self.polygon_in, polygon_out=self.polygon_out,
                           plot_legend=False, without_forces= wo_forces, plot_scale=1e5, label_edges = True)
        canvas = pltCanvas(fig)
        canvas.draw()
        self.canvas.figure.clf()  # Clear the figure
        self.layoutDraw.replaceWidget(self.canvas, canvas)
        self.canvas = canvas

    def add_edge(self):
        # State variables to store the clicked nodes
        self.node1 = None
        self.node2 = None

        def on_click(event):
            if event.xdata is not None and event.ydata is not None:
                curr_point = TS.Point(event.xdata, event.ydata,0)
                curr_point.change_abs_tol(150)
                print('Current Point', curr_point.printPoint())
                for n in self.Truss.NodeList:
                    if curr_point == n.point:
                        clicked_node = n
                        print('Clicked Node', n.point.printPoint())
                if self.node1 is None:
                    self.node1 = clicked_node
                elif self.node2 is None:
                    self.node2 = clicked_node
                    self.canvas.mpl_disconnect(cid)
                    # Process the nodes (e.g., add an edge)
                    self.Truss.addEdge([TS.Edge(self.node1, self.node2,mate=self.Truss.MaterialList[0])])
                    self.draw_truss()

        # Connect the event handler
        cid = self.canvas.mpl_connect("button_press_event", on_click)

    def add_force(self):
        def on_click(event):
            if event.xdata is not None and event.ydata is not None:
                curr_point = TS.Point(event.xdata, event.ydata,0)
                curr_point.change_abs_tol(150)
                print('Current Point', curr_point.printPoint())
                for n in self.Truss.NodeList:
                    if curr_point == n.point:
                        clicked_node = n
                        print('Clicked Node', n.point.printPoint())
                self.canvas.mpl_disconnect(cid)
                # get DOFs
                DOF_force = TS.DOF(0,0,0,0,0,0)
                dial_force = setForceDialog(self)
                dial_force.exec()
                Force = dial_force.getForce()
                Fx, Fy = Force[0], Force[1]
                if not math.isclose(Fx, 0):
                    DOF_force.Dx = 1
                if not math.isclose(Fy, 0):
                    DOF_force.Dy = 1
                self.Truss.addForce(TS.Force_ext(clicked_node, [Fx,Fy,0,0,0,0], direc_Force = DOF_force))
                self.draw_truss()

        # Connect the event handler
        cid = self.canvas.mpl_connect("button_press_event", on_click)

    def add_support(self):
        def on_click(event):
            if event.xdata is not None and event.ydata is not None:
                curr_point = TS.Point(event.xdata, event.ydata,0)
                curr_point.change_abs_tol(150)
                print('Current Point', curr_point.printPoint())
                for n in self.Truss.NodeList:
                    if curr_point == n.point:
                        clicked_node = n
                        print('Clicked Node', n.point.printPoint())
                self.canvas.mpl_disconnect(cid)
                # get DOFs
                dial_support = SetSupportDialog(self)
                dial_support.exec()
                DOFSupport = dial_support.get_support_conditions()
                if DOFSupport != TS.DOF(0,0,0,0,0,0):
                    if len(self.Truss.NodeList) == 0:
                        DOFSupport.Dz = 1
                    self.Truss.addSupport(TS.Support(clicked_node, DOFSupport))
                self.draw_truss()

        # Connect the event handler
        cid = self.canvas.mpl_connect("button_press_event", on_click)

    def add_node(self):
        try:
            dial_Point = AddPointDialog(self)
            dial_Point.setWindowTitle("Define coordinate position")
            dial_Point.resize(QSize(250, 150))
            if dial_Point.exec():
                x_coord, y_coord = dial_Point.get_coordinates()
                Point = TS.Point(float(x_coord), float(y_coord), 0)
            DOF0 = TS.DOF(0,0,0,0,0,0)
            n = TS.Node(Point, DOF0)
            self.Truss.addNode(n)
            self.draw_truss()
        except:
            self.display_message("Point is not inside the search space.", "Error", QMessageBox.Icon.Warning)

    def add_polygon_clicked(self, poly_in=True):
        dial_Polygon = AddPolygonDialog(self)
        dial_Polygon.exec()
        polygon = dial_Polygon.get_polygon()
        if poly_in:
            self.polygon_in.append(polygon)
            self.Truss.addPolygontoSS(polygon)
        else:
            self.polygon_out.append(polygon)
        self.draw_truss()

class EdgeMaterialDialog(QDialog):
    def __init__(self, *args, **kwargs):
        super(EdgeMaterialDialog, self).__init__(*args, **kwargs)
        self.setWindowTitle("Choose material of edge")
        self.resize(QSize(250, 150))

        layout = QVBoxLayout()

        self.material_index = QLineEdit()
        self.material_index.setMaxLength(1)
        self.material_index.setPlaceholderText("Material index, starting with 0")
        layout.addWidget(self.material_index)

        button_close = QPushButton("Done")
        button_close.clicked.connect(self.accept)
        layout.addWidget(button_close)

        self.setLayout(layout)

    def getMaterial(self):
        # index starts with 0
        return int(self.material_index.text())

class ThicknessDialog(QDialog):
    def __init__(self, *args, **kwargs):
        super(ThicknessDialog, self).__init__(*args, **kwargs)
        self.setWindowTitle("Define thickness")
        self.resize(QSize(250, 150))

        layout = QVBoxLayout()

        self.thickness = QLineEdit()
        self.thickness.setMaxLength(10)
        self.thickness.setPlaceholderText("thickness")
        layout.addWidget(self.thickness)

        button_close = QPushButton("Done")
        button_close.clicked.connect(self.accept)
        layout.addWidget(button_close)

        self.setLayout(layout)

    def getThickness(self):
        return float(self.thickness.text())

class setForceDialog(QDialog):
    def __init__(self,*args, **kwargs):
        super(setForceDialog, self).__init__(*args, **kwargs)
        self.setWindowTitle("Choose Forces")
        self.resize(QSize(250, 150))

        layout = QVBoxLayout()

        layout_force = QHBoxLayout()  # layout

        self.force_x = QLineEdit()
        self.force_x.setMaxLength(10)
        self.force_x.setPlaceholderText("Fx in x-direction")
        layout_force.addWidget(self.force_x)

        self.force_y = QLineEdit()
        self.force_y.setMaxLength(10)
        self.force_y.setPlaceholderText("Fy in y-direction")
        layout_force.addWidget(self.force_y)

        button_close = QPushButton("Done")
        button_close.clicked.connect(self.accept)

        layout.addLayout(layout_force)
        layout.addWidget(button_close)

        self.setLayout(layout)

    def getForce(self):
        return [float(self.force_x.text()), float(self.force_y.text())]

class SetSupportDialog(QDialog):
    def __init__(self, *args, **kwargs):
        super(SetSupportDialog, self).__init__(*args, **kwargs)
        self.setWindowTitle("Set Support Conditions")
        self.resize(QSize(250, 150))

        layout = QVBoxLayout()

        layout_support = QHBoxLayout() # layout

        self.support_dx = QCheckBox("Dx")
        layout_support.addWidget(self.support_dx)

        self.support_dy = QCheckBox("Dy")
        layout_support.addWidget(self.support_dy)

        button_close = QPushButton("Done")
        button_close.clicked.connect(self.accept)

        layout.addLayout(layout_support)
        layout.addWidget(button_close)

        self.setLayout(layout)

    def get_support_conditions(self):
        dx, dy = 0, 0
        if self.support_dx.isChecked():
            dx = 1
        if self.support_dy.isChecked():
            dy = 1
        return TS.DOF(dx, dy, 0, 0, 0, 0)


class AddPointDialog(QDialog):
    def __init__(self, *args, **kwargs):
        super(AddPointDialog, self).__init__(*args, **kwargs)
        self.setWindowTitle("Add Point")
        self.resize(QSize(250, 150))

        layout = QVBoxLayout() # layout

        # add text fields for x and y coordinates
        self.widget_x_coord = QLineEdit()
        self.widget_x_coord.setMaxLength(10)
        self.widget_x_coord.setPlaceholderText("X coordinate")
        self.widget_y_coord = QLineEdit()
        self.widget_y_coord.setMaxLength(10)
        self.widget_y_coord.setPlaceholderText("Y coordinate")

        layout.addWidget(self.widget_x_coord)
        layout.addWidget(self.widget_y_coord)

        button_close = QPushButton("Done")
        button_close.clicked.connect(self.accept)
        layout.addWidget(button_close)

        self.setLayout(layout)

    def get_coordinates(self):
        x_coord = self.widget_x_coord.text()
        y_coord = self.widget_y_coord.text()
        return x_coord, y_coord

    def add_point_clicked(self):
        print("Add point clicked")

class AddPolygonDialog(QDialog):
    def __init__(self, *args, **kwargs):
        super(AddPolygonDialog, self).__init__(*args, **kwargs)
        self.setWindowTitle("Add Polygon")
        self.resize(QSize(400, 200))

        layout = QVBoxLayout()  # layout

        self.Points = []

        button_add_point = QPushButton("Add point, clockwise")
        button_add_point.clicked.connect(self.add_point_clicked)
        layout.addWidget(button_add_point)

        button_close = QPushButton("Close")
        button_close.clicked.connect(self.accept)
        layout.addWidget(button_close)

        self.setLayout(layout)

    def add_point_clicked(self):
        dial_Point = AddPointDialog(self)
        dial_Point.setWindowTitle("Add Point")
        dial_Point.resize(QSize(250, 150))
        if dial_Point.exec():
            x_coord, y_coord = dial_Point.get_coordinates()
            self.Points.append(TS.Point(float(x_coord), float(y_coord), 0))

    def get_polygon(self):
        Polygon = TS.Polygon()
        Polygon.addPoint(self.Points)
        return Polygon

class pltCanvas(FigureCanvas):
    def __init__(self, fig = None):
        super(pltCanvas, self).__init__(fig)



app = QApplication(sys.argv) # no arguments then [] instead of sys.argv

window = MainWindow()
window.show() # Show

app.exec() # start the event loop

# need one main window
