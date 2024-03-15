"""
This module contains the methods to validate strut-and-tie models and to evaluate the hydrostatic nodal zones
Version 0.2: Includes validation of hydrostatic nodes for corresponding stress fields with concentrated struts and ties.
----------------
Older versions:
Version 0.1: Initial release, only includes validation of strut-and-tie model without nodal zones.
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
__version__ = '0.2'
__maintainer__ = 'Karin Yu'

import copy
import math
import numpy as np
from abc import ABC, abstractmethod
import stm_trusssystem as TS

# functions needed to solve the truss
def LocalStiffnessTruss(E, A, I, Li):
    #k = E*A/L*np.array([[1, 0,-1, 0], [0, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 0]])
    k = np.array([
            [E*A/Li, 0, 0, -E*A/Li, 0, 0],
            [0, 12*E*I/Li**3, 6*E*I/Li**2, 0, -12*E*I/Li**3, 6*E*I/Li**2],
            [0, 6*E*I/Li**2, 4*E*I/Li, 0, -6*E*I/Li**2, 2*E*I/Li],
            [-E*A/Li, 0, 0, E*A/Li, 0, 0],
            [0, -12*E*I/Li**3, -6*E*I/Li**2, 0, 12*E*I/Li**3, -6*E*I/Li**2],
            [0, 6*E*I/Li**2, 2*E*I/Li, 0, -6*E*I/Li**2, 4*E*I/Li]
        ])
    return k


def getRotationMat(cos_theta, sin_theta):
    # rotation in x-y plane
    R = np.array([[cos_theta, sin_theta, 0, 0, 0, 0], [-sin_theta, cos_theta, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
                  [0, 0, 0, cos_theta, sin_theta, 0], [0, 0, 0, -sin_theta, cos_theta, 0], [0, 0, 0, 0, 0, 1]])
    return R

#Abstract Base Class for all methods
class Truss_Solve_Method(ABC):
    validTruss: bool
    isinEQ: bool

    def __init__(self):
        self.isinEQ = None
        self.validTruss = None

    @abstractmethod
    def solveSystem(self, nodes, edges, numberOfDofs, forces, supports):
        pass

    @abstractmethod
    def checkvalidTruss(self):
        pass

# Direct Stiffness Method
# adapted from 'Finite Element Method I' course at ETH Zürich
# by Prof. Dr. Eleni Chatzi supervised by Dr. Konstantinos Tatsis from spring semester 2021
class DSM(Truss_Solve_Method):
    def solveSystem(self, nodes, edges, numberOfDofs, forces, supports):
        def get_local_matrices(edge):
            # 1. Calculate length
            xi, yi = edge.start_node.point.x, edge.start_node.point.y
            xj, yj = edge.end_node.point.x, edge.end_node.point.y
            L = np.sqrt((xj - xi) ** 2 + (yj - yi) ** 2)

            # 2. Calculate element local stiffness matrix
            ke_loc = LocalStiffnessTruss(edge.mat.E, edge.area, edge.mom_inertia, L)

            # 3. Evaluate the rotation matrix
            cos_theta, sin_theta = (xj - xi) / L, (yj - yi) / L

            R_mat = getRotationMat(cos_theta, sin_theta)
            return ke_loc, R_mat

        # initialize global stiffness matrix and global displacement vector
        K = np.zeros((numberOfDofs, numberOfDofs))
        U = np.zeros((numberOfDofs, 1))

        # iterate over all edges
        for e in edges:
            # 1.-3. get local stiffness matrix and rotation matrix
            ke_local, R = get_local_matrices(e)
            # 4. Rotate element stiffness matrix to global coordinate system
            ke_global = R.T.dot(ke_local).dot(R)
            i, j = e.start_id, e.end_id
            # 5. Assemble the element global stiffness to the system global stiffness
            K[3*i:3*i+3, 3*i:3*i+3] += ke_global[0:3, 0:3]
            K[3*i:3*i+3, 3*j:3*j+3] += ke_global[0:3, 3:6]
            K[3*j:3*j+3, 3*i:3*i+3] += ke_global[3:6, 0:3]
            K[3*j:3*j+3, 3*j:3*j+3] += ke_global[3:6, 3:6]

        #Get global force vector
        F = np.zeros((numberOfDofs, 1))
        for f in forces:
            ind = nodes.index(f.node)
            for i in range(3):
                F[3*ind+i] += f.Force_magnitude[i]

        #Apply boundary conditions
        allDofs = np.arange(0, numberOfDofs)
        restrainedDofs = list()
        for bc in supports:
            ind = nodes.index(bc.node)
            if bc.dof.Dx != 0: restrainedDofs.append(3*ind)
            if bc.dof.Dy != 0: restrainedDofs.append(3*ind+1)
            if bc.dof.Rz != 0: restrainedDofs.append(3*ind+2)
        freeDofs = np.setdiff1d(allDofs, restrainedDofs)

        #Generate stiffness matrices
        Kff = K[freeDofs, :][:, freeDofs]
        Kfr = K[freeDofs, :][:, restrainedDofs]
        Krf = K[restrainedDofs, :][:, freeDofs]
        Krr = K[restrainedDofs, :][:, restrainedDofs]
        if np.linalg.det(Kff) == 0:
            self.validTruss = False
        else:
            self.validTruss = True
        Ff = F[freeDofs, :]
        #Solve for displacements
        Uf = np.linalg.solve(Kff, Ff)
        Fr = Krf.dot(Uf)

        U[freeDofs] = Uf

        F[restrainedDofs] = Fr

        EdgeForces = np.zeros(len(edges))
        ShearForces = np.zeros((len(edges),2))
        Moments = np.zeros((len(edges),2))

        #Calculate element forces
        for c,e in enumerate(edges):
            i, j = e.start_id, e.end_id
            elementDofs = [3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2]
            Ue_global = U[elementDofs, :]
            # 1.-3. get local stiffness matrix and rotation matrix
            ke_local, R = get_local_matrices(e)

            # Rotate element global displacements to the local coordinate system
            Ue_local = R.dot(Ue_global)

            Fe_local = ke_local.dot(Ue_local)

            EdgeForces[c] = -1*Fe_local[0]
            ShearForces[c,0] = Fe_local[1]
            ShearForces[c, 1] = -1*Fe_local[4]
            Moments[c,0] = -1*Fe_local[2]
            Moments[c,1] = Fe_local[5]
        return EdgeForces, U, Fr, ShearForces, Moments

    def checkvalidTruss(self):
        return self.validTruss

# class to check nodal zones
class CheckNodalZone():
    def __init__(self, tol = 1e-4, print_disc_points = False, check_disc_points = False, polygon_in = None, polygon_out = None):
        self.tol = tol # tolerance for comparison
        self.print_discontinuity_points = print_disc_points
        self.check_discontinuity_points = check_disc_points
        if self.check_discontinuity_points:
            if polygon_in is None:
                raise ValueError("Polygons are not defined.")
            self.polygon_in = polygon_in
            self.polygon_out = polygon_out

    def applybtoEdge(self, d, positiveside=True):
        # d = [endpoint_x, endpoint_y, ux, uy, width]
        # moves edge to the positive or negative side --> show strut thickness
        # positive side is right side
        # for each strut it takes into account the width of the strut and adds it
        new = copy.deepcopy(d)
        if positiveside:
            move = d[4] / 2  # width of strut/2 and direction +
        else:
            move = -1 * d[4] / 2  # width of strut/2 and direction -
        perpendicular_direction = [-1 * d[3], d[2]]  # [-uy, ux]
        new[0] = d[0] + move * perpendicular_direction[0]
        new[1] = d[1] + move * perpendicular_direction[1]
        return new

    def reduceNumEdges(self, edges):
        # if there are more than 3 edges we need to reduce the edges
        # returns new edges, either intersection point (2 Ties) or distance to intersection point from center point, original edge numbering, case, new edge index
        # cases: '1' two ties, '2' else
        # 1. Check whether two tension ties are in one plane, replace by one
        for i in range(len(edges) - 1):
            if edges[i].type > 0:
                for j in range(i + 1, len(edges)):
                    if edges[j].type > 0 and math.isclose(abs(edges[i].line.angle - edges[j].line.angle), math.pi):
                        print('Two ties are on the same plane! One tie is reduced.')
                        # replace both edges by one, remove the larger one (as opposite is compressive force)
                        if edges[i].force <= edges[j].force:
                            edges[i].UpdateForce(edges[i].force - edges[j].force)
                            edges[i].UpdateAreaAuto()
                            end_point = edges[j].end_node.point
                            edges.remove(edges[j])
                            return edges, end_point, j, 1, i
                        else:
                            edges[j].UpdateForce(edges[j].force - edges[i].force)
                            edges[j].UpdateAreaAuto()
                            end_point = edges[i].end_node.point
                            edges.remove(edges[i])
                            return edges, end_point, i, 1, i
        # 2. determine two "smallest" edges by choosing the ones with the smallest angle between the lines
        min_sum = TS.BoundAngles(edges[0].line.angle - edges[-1].line.angle)
        min_pair = [-1, 0]
        for i in range(len(edges)-1):
            current_sum = edges[i+1].line.angle - edges[i].line.angle
            if current_sum < min_sum: # and current_sum > 0
                min_sum = current_sum
                min_pair = [i, i+1]
        # add the two edges
        UVL1 = TS.UnitVector(edges[min_pair[0]].line)
        UVL2 = TS.UnitVector(edges[min_pair[1]].line)
        Ftot = [edges[min_pair[0]].force * UVL1.ux + edges[min_pair[1]].force * UVL2.ux,
                edges[min_pair[0]].force * UVL1.uy + edges[min_pair[1]].force * UVL2.uy,
                edges[min_pair[0]].force * UVL1.uz + edges[min_pair[1]].force * UVL2.uz]
        angle = TS.getAngle(Ftot[0], Ftot[1]) + math.pi #only in x-y plane, compression is positive --> add pi
        # since all forces were converted into compressive forces
        sgn = -1
        # define new start_node
        L_new_start_node = TS.Node(TS.Point(edges[min_pair[0]].end_node.point.x +sgn* math.cos(angle)*100, edges[min_pair[0]].end_node.point.y +sgn* math.sin(angle)*100, edges[min_pair[0]].end_node.point.z))
        # define equivalent edge
        Fmag = math.sqrt(sum([f**2 for f in Ftot]))
        area_tot = abs(Fmag / edges[min_pair[0]].mat.fy)
        new_edge = TS.Edge(L_new_start_node, edges[0].end_node, edges[min_pair[0]].mat, A = area_tot, F = sgn*Fmag, b = area_tot/(edges[min_pair[0]].area/edges[min_pair[0]].width))
        # replace edge
        if min_pair[0] >= 0:
            new_edges = edges[:min_pair[0]]+[new_edge] + edges[min_pair[1]+1:]
        else:
            new_edges = edges[min_pair[1]+1:min_pair[0]] + [new_edge]
        # calculate relative discontinuous intersection point from center of new edge [-uy, ux]
        # width_new / 2 + width_ol
        UV_new_edge = TS.UnitVector(new_edge.line)
        dp1 = TS.Point(-1*(-UVL1.uy * edges[min_pair[0]].width) + (-UV_new_edge.uy * new_edge.width / 2),
                       -1*(UVL1.ux * edges[min_pair[0]].width) + (UV_new_edge.ux * new_edge.width / 2), 0)
        dp2 = TS.Point((-UVL2.uy * edges[min_pair[1]].width) - (-UV_new_edge.uy * new_edge.width / 2),
                       (UVL2.ux * edges[min_pair[1]].width) - (UV_new_edge.ux * new_edge.width / 2), 0)
        # dp1 and dp2 should be equal
        if dp1 != dp2:
            raise ValueError("The distances of the discontinuous points are not equal.")
        return new_edges, dp1, min_pair[0], 2, min_pair[0]


    def HydrostaticNZCheck(self, edges, fck):
        print('Hydrostatic nodal zones are assumed.')
        # checks whether hydrostatic stress is fulfilled in nodes
        # 1. all end nodes must be the same
        edges = copy.deepcopy(edges)
        for e in edges:
            if e.force > 0:
                e.line.mirror(False)
                e.force = -1*e.force
            if e.mat.fy != fck: # change force
                e.mat.fy = float(fck)
                e.UpdateAreaAuto()
        edges_indices = sorted(range(len(edges)), key = lambda index: edges[index].line.angle)
        edges = [edges[edges_indices[i]] for i in range(len(edges_indices))]
        # 2. if there are more than three intersecting edges they must be reduced to three
        # relative_disc_inter_p_red stores information about reduced edges
        relative_disc_inter_p_red = []
        #print('Current number of edges:', len(edges))
        while len(edges) > 3:
            # returns new edges, either intersection point (2 Ties) or distance to intersection point from center point, original edge numbering, case, new edge index
            edges, rel_disc_intersect_p, ind, case, new_edge_index = self.reduceNumEdges(edges)
            relative_disc_inter_p_red.append([rel_disc_intersect_p,ind, case])
            if edges[new_edge_index].force > 0:
                edges[new_edge_index].line.mirror(False)
                edges[new_edge_index].force = -1*edges[new_edge_index].force
        #print('Reduced number of edges:', len(edges))
        # 3. check for hydrostatic node
        # Nodal zone check for hydrostatic nodes
        if fck < 0:
            raise ValueError("fck should be a positive input parameter.")
        #since hydrostatic, we just need to check force equilibrium in the node
        fx, fy, fz = 0, 0, 0
        min_sig = edges[0].force/edges[0].area
        direc = []
        for e in edges:
            e.getprojectedForce()
            fx = fx + e.fx
            fy = fy + e.fy
            fz = fz + e.fz
            if abs(abs(e.force/e.area) - abs(min_sig)) > self.tol:
                print('Min stress: ', min_sig, 'Curr stress: ', abs(e.force/e.area))
                raise ValueError("This is not a hydrostatic node.")
            min_sig = min(e.force/e.area, min_sig)
            direc.append([e.end_node.point.x, e.end_node.point.y, e.line.ux, e.line.uy, e.width])
        # 4. determine discontinuity points (borders of nodal zone)
        disc_points = []
        for i in range(len(edges)):
            disc_points.append(self.generate_disc_points(direc, edges, i))
        if abs(fx)+abs(fy)+abs(fz) < self.tol: # check for force equilibrium with certain tolerance
            for i in range(len(relative_disc_inter_p_red)-1, -1, -1):
                # relative_disc_inter_p_red contains [either intersection point (2 Ties) or distance to intersection point from center point, original edge numbering, case]
                if relative_disc_inter_p_red[i][2] == 1: # two ties
                    if relative_disc_inter_p_red[i][1] >= len(disc_points):
                        disc_points = disc_points[:relative_disc_inter_p_red[i][1]-1] + [
                            disc_points[relative_disc_inter_p_red[i][1]-1]] + [disc_points[0]]
                    else:
                        disc_points = disc_points[:relative_disc_inter_p_red[i][1] + 1] +  [disc_points[relative_disc_inter_p_red[i][1]]] + disc_points[relative_disc_inter_p_red[i][1] + 1:]
                else: # other cases
                    # find discontinuity point between reduced edges
                    if 0 <= relative_disc_inter_p_red[i][1] < len(disc_points) - 1:
                        equiv_line = TS.Line(disc_points[relative_disc_inter_p_red[i][1]], disc_points[relative_disc_inter_p_red[i][1]+1])
                    else:
                        equiv_line = TS.Line(disc_points[relative_disc_inter_p_red[i][1]], disc_points[0])
                    center_equiv_line = TS.getCenterofLine(equiv_line)
                    disc_intersect = TS.movePoint(center_equiv_line, relative_disc_inter_p_red[i][0])
                    disc_points = disc_points[:relative_disc_inter_p_red[i][1] + 1] + [disc_intersect] + disc_points[relative_disc_inter_p_red[i][1] + 1:]
            # include if you want to print the locations of the discontinuity points
            if self.print_discontinuity_points:
                print('Discontinuity points:')
                for p in disc_points:
                    print(p.printPoint())
            if self.check_discontinuity_points:
                for poly in self.polygon_in:
                    for p in disc_points:
                        if not poly.containsPoint(p):
                            print('Discontinuity point p: ', p.printPoint(), ' is not inside the geometry.')
                if not self.polygon_out is None:
                    for poly in self.polygon_out:
                        for p in disc_points:
                            if not poly.containsPoint(p):
                                print('Discontinuity point p: ', p.printPoint(), ' is inside an opening.')
            return min_sig >= -1*fck-self.tol, disc_points, edges_indices
        else: # if force equilibrium is not fulfilled
            print('Here: fx: {:.2f}, fy: {:.2f}, fz: {:.2f}, given the minimum stress: {:.2f}'.format(fx, fy, fz, min_sig))
            print('If the stresses are only slightly off, consider changing the tolerance tol to a higher value.')
            raise ValueError("This node is not in equilibrium.")

    def generate_disc_points(self, direc, edges, i):
        # generates discontinuity points based on the direction of the edges, the edges and the index of the edge
        # move edges by edge width/2 of neighboring struts
        edge1 = self.applybtoEdge(direc[i - 1], positiveside=False)
        edge2 = self.applybtoEdge(direc[i], positiveside=True)
        # determine intersection point through linear algebra
        A = np.array([[edge1[2], -1 * edge2[2]], [edge1[3], -1 * edge2[3]]])
        if not np.linalg.det(A):
            lmd = np.array([[1]])
        else:
            b = np.array([[edge2[0] - edge1[0]], [edge2[1] - edge1[1]]])
            lmd = np.dot(np.linalg.inv(A), b)
        return TS.Point(x=edge1[0] + lmd[0, 0] * edge1[2], y=edge1[1] + lmd[0, 0] * edge1[3],
                        z=edges[i - 1].end_node.point.z)



