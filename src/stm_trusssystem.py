"""
This module contains functions to build the truss/strut-and-tie model.
Version 0.1: Initial release, only includes validation of strut-and-tie models without nodal zones.
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
__version__ = '0.1'
__maintainer__ = 'Karin Yu'

import numpy as np
import math
import matplotlib.pyplot as plt
from stm_methods import *

# Define colors
color_map = {
    'g': 'tab:green',
    'b': 'tab:blue',
    'r': 'tab:red',
    'orange': 'tab:orange',
    'grey': 'tab:gray',
}

def BoundAngles(angle, tol = 1e-4):
    while angle >= 2*math.pi-tol:
        angle = angle - 2*math.pi
    while angle < 0:
        angle =  angle + 2 * math.pi
    return angle

def getAngle(dx: float, dy: float):
    if dx == 0:
        if dy > 0:
            return math.pi / 2
        elif dy < 0:
            return math.pi / 2 * 3
        else:
            return 0
    elif dx > 0:
        return math.atan(dy / dx)
    else:
        return math.atan(dy / dx) + math.pi

# Geometric Objects
class Point:
    x: float
    y: float
    z: float

    def __init__(self, x: float, y: float, z: float):
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, other):
        return all([abs(self.x - other.x) < 1e-5, abs(self.y - other.y) < 1e-5, abs(self.z - other.z) < 1e-5])

    def printPoint(self, verbose = False):
        if verbose: print('x: {:.2f}, y: {:.2f}, z: {:.2f}'.format(self.x, self.y, self.z))
        return [self.x, self.y, self.z]

class Line:
    def __init__(self, start: Point, end: Point):
        self.start = start
        self.end = end
        self.length = ((self.start.x-self.end.x)**2+(self.start.y-self.end.y)**2+(self.start.z-self.end.z)**2)**(1/2)
        self.updateAngle()

    def mirror(self, start_Point = True):
        direc = [self.end.x - self.start.x, self.end.y - self.start.y]
        # start_Point, mirror end point around start point
        if start_Point:
            self.end.x = self.end.x - 2*direc[0]
            self.end.y = self.end.y - 2 * direc[1]
            self.angle = BoundAngles(self.angle + math.pi)
        else:
            # not start_Point, mirror start point around end point
            self.start.x = self.start.x + 2*direc[0]
            self.start.y = self.start.y + 2 * direc[1]
            self.angle = BoundAngles(self.angle + math.pi)

    def updateAngle(self):
        UV = UnitVector(self)
        self.angle = UV.angle

class UnitVector:
    def __init__(self, L: Line):
        if L.length == 0:
            self.ux = 0
            self.uy = 0
            self.uz = 0
            self.angle = 0
        else:
            self.ux = (L.end.x-L.start.x)/L.length
            self.uy = (L.end.y-L.start.y)/L.length
            self.uz = (L.end.z-L.start.z)/L.length
            self.angle = getAngle(self.ux, self.uy)
            self.angle = BoundAngles(self.angle)

    def opposite(self):
        self.ux, self.uy, self.uz = -1*self.ux, -1*self.uy, -1*self.uz

    def printVec(self):
        return [self.ux, self.uy, self.uz]

#To check geometry
def onLine(p: Point, L: Line):
    #checks whether point p is on line L
    if (p.x <= max(L.start.x, L.end.x)) and (p.x >= min(L.start.x, L.end.x)) and (p.y >= min(L.start.y, L.end.y)) and (p.y <= max(L.start.y, L.end.y)) and (p.z >= min(L.start.z, L.end.z)) and (p.z <= max(L.start.z, L.end.z)):
        temp = Line(p, L.end)
        Oppo = UnitVector(temp)
        Oppo.opposite()
        UV_L = UnitVector(L)
        return all((math.isclose(n1, n2, abs_tol = 1e-5) for n1, n2 in zip(UV_L.printVec(),  Oppo.printVec()))) or all((math.isclose(n1, n2, abs_tol = 1e-5) for n1, n2 in zip(UV_L.printVec(),  UnitVector(temp).printVec())))
    else:
        return False

def Intersect(L1: Line, L2: Line):
    # checks whether two lines L1 and L2 intersect, if they are on the same plane (only differences in x and y)
    a = L1.end.x-L1.start.x
    b = -(L2.end.x-L2.start.x)
    c = L1.end.y-L1.start.y
    d = -(L2.end.y-L2.start.y)
    bx = L2.start.x - L1.start.x
    by = L2.start.y - L1.start.y
    if a*d-b*c == 0: return False
    x = 1/(a*d-b*c)*np.array([[d, -b], [-c, a]]).dot(np.array([bx, by]).T)
    return all([True if 0 <= i <= 1 else False for i in x])

# Material
class Material:
    E: float  # Young's modulus in MPa
    fy: float  # yield strength in MPa
    eps_u: float  # ultimate strain
    eps_y: float  # yield strain
    isotropic: bool  # is it isotropic?
    direc: int  # -1: only compression, 0: both, 1: only tension


    def __init__(self, E: float, fy: float, eu: float, d: int):
        self.isotropic = d == 0
        self.E = E
        self.fy = fy
        self.eps_u = eu
        self.eps_y = fy / E
        self.direc = d


# Element Stuff
class DOF:  # mkr
    Dx: float
    Dy: float
    Dz: float
    Rx: float
    Ry: float
    Rz: float

    def __init__(self, Dx: float, Dy: float, Dz: float, Rx: float, Ry: float, Rz: float):
        # 1 if direction of axis is in tension, -1 if direction of axis is in compression
        self.Dx = Dx
        self.Dy = Dy
        self.Dz = Dz
        self.Rx = Rx
        self.Ry = Ry
        self.Rz = Rz

    def __eq__(self, other):
        return all(
            [self.Dx == other.Dx, self.Dy == other.Dy, self.Dz == other.Dz, self.Rx == other.Rx, self.Ry == other.Ry,
             self.Rz == other.Rz])


class Node:
    point: Point
    dof: DOF
    stress: list
    force: list  # Fx, Fy, Fz
    isSupport: bool

    # attributes

    def __init__(self, p: Point, DOF_node: DOF = DOF(0, 0, 0, 0, 0, 0)):
        self.point = p
        self.dof = DOF_node
        self.isSupport = False

    def __eq__(self, other):
        if isinstance(other, Node):
            return self.point == other.point
        else:
            return False


class Edge:
    start_node: Node
    end_node: Node
    mat: Material
    line: Line
    # specific to reinforced concrete
    type: int  # -1: strut, 0: 0 force,  1: tie
    # attributes
    area: float
    width: float
    force: float
    mom_inertia: float
    stress: float
    start_id: int
    end_id: int
    mat_id: int
    vol: float

    def __init__(self, start: Node, end: Node, mate: Material, A: float = 1, F: float = 0, I: float = 1e-6,
                 b: float = 0.5):
        self.start_node = start
        self.end_node = end
        self.mat = mate
        self.line = Line(start.point, end.point)
        self.area = A
        self.mom_inertia = I
        self.force = F
        self.stress = F / A
        self.type = np.sign(F)
        self.width = b
        self.end_id = None
        self.start_id = None

    def __eq__(self, other):
        if (self.start_node == other.start_node) and (self.end_node == other.end_node):
            return True
        elif (self.start_node == other.end_node) and (self.end_node == other.start_node):
            return True
        else:
            return False

    def CheckStress(self):
        if self.stress < 0:
            if self.mat.direc <= 0:
                return abs(self.stress) <= self.mat.fy
            else:
                raise ValueError('This truss is being compressed but the material cannot take compression.')
        else:
            if self.mat.direc >= 0:
                return self.stress <= self.mat.fy
            else:
                raise ValueError('This truss is a tie but the material cannot be under tension.')

    def UpdateForce(self, F: float = 0):
        self.force = F
        self.stress = F / self.area

    def UpdateAreaAuto(self, A_min=0):
        t = self.area / self.width
        self.area = max(abs(self.force / self.mat.fy), A_min)
        self.mom_inertia = self.area ** 2 / (4 * np.pi)
        self.stress = self.force / self.area
        self.width = self.area / t
        self.type = np.sign(self.force)

    @property
    def vol(self):
        return self.area * self.line.length

    def getprojectedForce(self):
        unit_vec = UnitVector(self.line)
        self.fx, self.fy, self.fz = unit_vec.ux * self.force, unit_vec.uy * self.force, unit_vec.uz * self.force

    def UpdateLine(self):
        if self.line.start != self.start_node.point:
            self.line.end, self.line.start = self.line.start, self.line.end
            self.line.updateAngle()


# Boundary Conditions
class Support:
    def __init__(self, n: Node, DOF_support: DOF):
        self.node = n
        self.dof = DOF_support

    def __eq__(self, other):
        if type(other) is Node:
            return self.node == other
        else:
            return (self.node == other.node) and (self.dof == other.dof)


class Force_ext:
    def __init__(self, n: Node, F_mag, direc_Force: DOF = DOF(0, 0, 0, 0, 0, 0)):
        self.node = n
        self.Force_magnitude = F_mag  # corresponds to Fx, Fy, Fz, Mx, My, Mz
        self.dof = direc_Force

    def __eq__(self, other):
        if type(other) is Node:
            return self.node == other
        else:
            return (self.node == other.node) and (self.Force_magnitude == other.Force_magnitude)

# Polygon
class Polygon:  # assume all points are on one plane
    Points: list  # need to be in clockwise direction
    numEdges: int
    isPolygon: bool

    def __init__(self):
        self.Points = list()
        self.numEdges = 0

    def addPoint(self, pts):
        if type(pts) != list: pts = [pts]
        for p in pts:
            if p not in self.Points:
                self.Points.append(p)
            else:
                raise ValueError('Point p is already part of search space.')
            if len(self.Points) == 2:
                self.numEdges += 1
            elif len(self.Points) > 2:
                self.numEdges = len(self.Points)
            else:
                self.numEdges = 0

    def removePoint(self, p: Point):
        if p not in self.Points:
            raise ValueError('Point p is not part of search space.')
        else:
            self.Points.remove(p)
            if len(self.Points) >= 3:
                self.numEdges -= 1
            elif len(self.Points) == 2:
                self.numEdges -= 2
            elif len(self.Points) <= 1:
                self.numEdges = 0

    @property
    def isPolygon(self):
        return self.numEdges >= 3

    def RotationMatrixforCoordinateTransform(self):
        # so that z coordinate is equal
        if self.isPolygon:
            L1 = Line(self.Points[0], self.Points[1])
            L2 = Line(self.Points[1], self.Points[2])
            x, y, z = L1.end.x - L1.start.x, L1.end.y - L1.start.y, L1.end.z - L1.start.z
            a, b, c = L2.end.x - L2.start.x, L2.end.y - L2.start.y, L2.end.z - L2.start.z
            UV_Vec = UnitVector(Line(Point(0, 0, 0), Point(y * c - z * b, z * a - x * c, x * b - y * a)))
            UV_Vec.opposite()
            normVec = np.array([UV_Vec.printVec()]).T
            Vec_z = np.array([[0, 0, 1]]).T
            n_n = np.cross(np.squeeze(normVec), np.squeeze(Vec_z))
            alpha = np.arccos(normVec.T.dot(Vec_z))
            self.RotationMatrix = n_n.dot(n_n.T) + np.cos(alpha) * (np.eye(3) - n_n.dot(n_n.T)) + np.sin(
                alpha) * np.cross(n_n, np.identity(3) * -1)
            # Check: https://en.wikipedia.org/wiki/Rotation_matrix: Rotation matrix from axis and angle
        else:
            raise ValueError('It is not a polygon.')

    def LineCoordinateTransform(self, L: Line):
        P1 = np.array([L.start.printPoint()]).T
        P2 = np.array([L.end.printPoint()]).T
        temp1 = np.squeeze(self.RotationMatrix.dot(P1))
        temp2 = np.squeeze(self.RotationMatrix.dot(P2))
        return Line(Point(temp1[0], temp1[1], temp1[2]), Point(temp2[0], temp2[1], temp2[2]))

    #checks if point is inside polygon
    def containsPoint(self, p: Point, verbose=False):
        if self.isPolygon:
            self.RotationMatrixforCoordinateTransform()
            # Ray-casting algorithm
            Ray = Line(p, Point(self.Points[0].x, self.Points[0].y - 100, self.Points[0].z))
            Ray = self.LineCoordinateTransform(Ray)
            currLine = self.LineCoordinateTransform(Line(self.Points[-1], self.Points[0]))
            count = 0
            for i in range(len(self.Points)):
                if onLine(p, currLine): return True
                if Intersect(currLine, Ray):
                    count += 1
                if i < len(self.Points) - 1: currLine = self.LineCoordinateTransform(
                    Line(self.Points[i], self.Points[i + 1]))
            if verbose: print('Number of Counts: ', count)
            return count % 2 == 1
        else:
            raise ValueError('Search space is not a polygon.')


class TrussSystem:
    SearchSpace: list  # list of polygons
    MaterialList: list  # list of materials
    NodeList: list  # list of nodes
    EdgeList: list  # list of edges / trusses
    ForceList: list  # list of forces
    SupportList: list  # list of supports
    issolved: bool

    def __init__(self):
        self.SearchSpace = list()
        self.MaterialList = list()
        self.NodeList = list()
        self.EdgeList = list()
        self.ForceList = list()
        self.SupportList = list()
        self.solvermethod = DSM()
        self.ndim = 3  # 2 dofs per node
        self.numDofs = 0
        self.issolved = False

    def CheckNodeInPolygon(self, n: Node):
        for poly in self.SearchSpace:
            if poly.containsPoint(n.point, False): return True
        return False

    def addPolygontoSS(self, Poly: Polygon):  # , mat: Material
        if Poly not in self.SearchSpace:
            self.SearchSpace.append(Poly)

    def addMaterial(self, mat: Material):
        if mat not in self.MaterialList:
            self.MaterialList.append(mat)

    def autoAllocateMat(self):
        # 1. Material is for compression & second Material is for tension
        if not self.issolved: raise ValueError('System has not been solved.')
        for i in range(len(self.EdgeList)):
            if self.EdgeList[i].force <= 0:
                self.EdgeList[i].mat_id = 0
            else:
                self.EdgeList[i].mat_id = 1

    def addNode(self, nodes):
        if type(nodes) != list: nodes = [nodes]
        for n in nodes:
            if n not in self.NodeList:
                if self.CheckNodeInPolygon(n):
                    self.NodeList.append(n)
                    self.numDofs += self.ndim
                else:
                    raise ValueError('Node ', n.point.printPoint(), ' is not inside the defined search space.')

    def removeNode(self, nodes):
        if type(nodes) != list: nodes = [nodes]
        for n in nodes:
            if n in self.NodeList:
                self.NodeList.remove(n)
            else:
                raise ValueError('Node ', n.point.printPoint(), ' is not in NodeList.')
        self.UpdateEdgeIDs()

    def changeNodePosition(self, old_loc: Point, new_loc: Point):
        NewNode = Node(old_loc)
        if not NewNode in self.NodeList: raise ValueError('The given node is not in the node list.')
        node_id = self.NodeList.index(NewNode)
        self.NodeList[node_id].point = new_loc

    def addEdge(self, edges):
        for e in edges:
            if e not in self.EdgeList:
                if e.start_node == e.end_node:
                    raise ValueError('Start node coincides with end node.')
                if (e.start_node in self.NodeList) and (e.end_node in self.NodeList):
                    if e.mat in self.MaterialList:
                        e.start_id = self.NodeList.index(e.start_node)
                        e.end_id = self.NodeList.index(e.end_node)
                        e.mat_id = self.MaterialList.index(e.mat)
                        self.EdgeList.append(e)
                    else:
                        raise ValueError('Material has not been defined in MaterialList.')
                else:
                    raise ValueError('Nodes of edges are not in NodeList.')

    def removeEdge(self, e: Edge):
        if e in self.EdgeList:
            self.EdgeList.remove(e)
        else:
            raise ValueError('Edge is not in EdgeList.')

    def UpdateEdgeIDs(self):
        for e in self.EdgeList:
            e.start_id = self.NodeList.index(e.start_node)
            e.end_id = self.NodeList.index(e.end_node)
            e.mat_id = self.MaterialList.index(e.mat)

    def sortEdgeList(self):
        self.EdgeList.sort(key=lambda obj: obj.start_id)

    def addForce(self, f: Force_ext):
        if f not in self.ForceList:
            self.ForceList.append(f)
        if f.node not in self.NodeList:
            for i in range(len(self.EdgeList)):
                if onLine(f.node.point, self.EdgeList[i].line):
                    self.NodeList = self.NodeList[:self.EdgeList[i].start_id + 1] + [f.node] + self.NodeList[
                                                                                               self.EdgeList[i].end_id:]
                    self.numDofs += self.ndim
                    e1 = Edge(self.EdgeList[i].start_node, f.node, mate=self.EdgeList[i].mat)
                    e2 = Edge(f.node, self.EdgeList[i].end_node, mate=self.EdgeList[i].mat)
                    self.EdgeList = self.EdgeList[:i] + [e1, e2] + self.EdgeList[i + 1:]
                    break
        self.UpdateEdgeIDs()

    def addSupport(self, sup: Support):
        if sup not in self.SupportList:
            if sup.node in self.NodeList:
                sup.id = self.NodeList.index(sup.node)
                self.NodeList[sup.id].isSupport = True
                self.NodeList[sup.id].dof = sup.dof
                self.SupportList.append(sup)
            else:
                raise ValueError('Support node is not in NodeList.')

    def solveTruss(self, update=False, A_min=0, d_max=None):
        self.EdgeForces, self.nodal_disp, self.restrainedForces, self.ShearForces, self.BendingMoments = self.solvermethod.solveSystem(
            self.NodeList, self.EdgeList, self.numDofs, self.ForceList, self.SupportList)
        if len(self.nodal_disp[0]) == 1:
            self.nodal_disp = np.reshape(self.nodal_disp, (len(self.nodal_disp) // 3, 3))
        if update:
            for i in range(len(self.EdgeList)):
                self.EdgeList[i].UpdateForce(self.EdgeForces[i])
                self.EdgeList[i].UpdateAreaAuto(A_min)
        self.issolved = True

    def ValidTruss(self):
        if self.getStaticDeterminancy() < 0:
            raise ValueError(
                'The truss is statically overdetermined (Static determinancy = {:.0f}). Make sure that the structure has a static determinancy of at least 0.'.format(
                    self.statdet))
        return self.solvermethod.checkvalidTruss()

    def checkEquilibrium(self):
        if not self.issolved: raise ValueError('System has not been solved.')
        Nod = [[0, 0, 0] for _ in self.NodeList]
        for e in self.EdgeList:
            e.getprojectedForce()
            Nod[e.start_id][0] -= e.fx
            Nod[e.start_id][1] -= e.fy
            Nod[e.start_id][2] -= e.fz
            Nod[e.end_id][0] += e.fx
            Nod[e.end_id][1] += e.fy
            Nod[e.end_id][2] += e.fz
        res = 0
        for i in range(len(self.NodeList)):
            if not self.NodeList[i].isSupport:
                res += abs(self.NodeList[i].isSupport)
        return math.isclose(res, 0)

    def getStaticDeterminancy(self):
        # TODO: Implement static determinancy for 3D
        # implemented for 2D
        n = len(self.NodeList)  # number of nodes
        r = 0  # initialize --> number of reaction forces
        for support in self.SupportList:  # currently in 2D
            if support.dof.Dx != 0: r += 1
            if support.dof.Dy != 0: r += 1
            # if support.dof.Dz != 0: r += 1
            # if support.dof.Rx != 0: r += 1
            # if support.dof.Ry != 0: r += 1
            if support.dof.Rz != 0: r += 1
        p = len(self.EdgeList)  # number of edges
        self.statdet = p + r - 2 * n
        return self.statdet

    def plotSTM(self, fig_size=(12, 6), savefig=None, ForceList=None, plot_scale=1e4, draw_geometry=True,
                without_forces=False, label_edges = False):
        plt.figure(figsize=fig_size)
        if self.SearchSpace:
            ymin, ymax = -500, 500
            for poly in self.SearchSpace:
                xval = []
                yval = []
                for p in poly.Points:
                    xval.append(p.x)
                    yval.append(p.y)
                xval.append(xval[0])
                yval.append(yval[0])
                if draw_geometry:
                    plt.plot(xval, yval, 'k-', linewidth=1)
                ymin = min([0.9 * min(yval), ymin, 1.1 * min(yval)])
                ymax = max([1.1 * max(yval), ymax])
            plt.ylim([ymin, ymax])
        for n in self.NodeList:
            plt.plot(n.point.x, n.point.y, '.', color='#bbbbbb', markersize=10)
        for i in range(len(self.EdgeList)):
            if ForceList is None:
                f = self.EdgeList[i].force
            else:
                if len(ForceList[i]) > 1:
                    f = ForceList[i][0]
                else:
                    f = ForceList[i]
            if without_forces:
                c = 'k'
                line_t = 1.5
            else:
                if f < 0:
                    c = color_map['g']
                    line_t = -1 * f / plot_scale
                elif f == 0:
                    c = 'k'
                    line_t = 1.5
                else:
                    c = color_map['b']
                    line_t = f / plot_scale
            plt.plot([self.EdgeList[i].start_node.point.x, self.EdgeList[i].end_node.point.x],
                     [self.EdgeList[i].start_node.point.y, self.EdgeList[i].end_node.point.y], '-', color=c,
                     linewidth=line_t)
            if label_edges:
                x_text = (self.EdgeList[i].start_node.point.x+ self.EdgeList[i].end_node.point.x)/2
                y_text = (self.EdgeList[i].start_node.point.y+ self.EdgeList[i].end_node.point.y)/2
                plt.text(x_text, y_text, '{:.0f}'.format(f), fontsize=9, verticalalignment='bottom')
        for i in range(len(self.ForceList)):
            f = self.ForceList[i]
            marker_size = 15
            if f.node in self.SupportList:
                marker_size = 0
            if i == 0 and marker_size > 0:
                plt.plot(f.node.point.x, f.node.point.y, '.', color=color_map['orange'], markersize=marker_size, label = 'Force')
            else:
                plt.plot(f.node.point.x, f.node.point.y, '.', color=color_map['orange'], markersize=marker_size)
        for j in range(len(self.SupportList)):
            s = self.SupportList[j]
            if j == 0:
                plt.plot(s.node.point.x, s.node.point.y, 'x', color=color_map['grey'], markersize=10, label = 'Support')
            else:
                plt.plot(s.node.point.x, s.node.point.y, 'x', color=color_map['grey'], markersize=10)
        plt.legend(loc='upper right')
        if savefig is not None:
            plt.savefig(savefig, dpi=600)
