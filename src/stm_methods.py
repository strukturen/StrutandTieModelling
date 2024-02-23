"""
This module contains the methods to validate strut-and-tie models.
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
from abc import ABC, abstractmethod

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

