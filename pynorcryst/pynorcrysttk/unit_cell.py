import numpy as np
from enum import Enum
from collections import namedtuple

from crystal_exceptions import *


class UnitCell(object):

    def __init__(self, lattice):
        self.volume = None
        #Angles in radians for convenience
        self.al_r, self.be_r, self.ga_r = None, None, None
        #Lattices
        self.lattice, self.reciprocal_lattice = None, None
        self.crystal_family = None
        self.orthonormalisation_matrix = None
        #Tensors
        self.metric_tensor, self.reciprocal_metric_tensor = None, None
    
        #Set the lattices up values
        self.update_cell(lattice)

    def update_cell(self, lattice):
        self.lattice = self.__evaluate_lattice(lattice)
        self.al_r, self.be_r, self.ga_r = map(np.radians, [self.lattice.al, 
                                                           self.lattice.be,
                                                           self.lattice.ga])
        self.metric_tensor = self.__determine_metric_tensor()
        self.volume = np.sqrt(np.linalg.det(self.metric_tensor))
        self.crystal_family = self.__get_crystal_family_for_lattice_system()
        self.reciprocal_metric_tensor = np.linalg.inv(self.metric_tensor)
        self.reciprocal_volume = np.sqrt(np.linalg.det(self.reciprocal_metric_tensor))
        self.reciprocal_lattice = self.__determine_reciprocal_lattice()

        self.orthonormalisation_matrix = np.matrix([[self.lattice.a,0,0],
                                                   [self.lattice.b*np.cos(self.ga_r), self.lattice.b*np.sin(self.ga_r),0],
                                                   [self.lattice.c*np.cos(self.be_r), -self.lattice.c*np.sin(self.be_r)*np.cos(np.radians(self.reciprocal_lattice.al)), 1/self.reciprocal_lattice.c]])


    def __evaluate_lattice(self, lattice):
        a, b, c = lattice.a, lattice.b, lattice.c
        angles = [lattice.al, lattice.be, lattice.ga]
        p_axis = PrincipleAxis.NONE
        lat_sys = CrystalSystem.TRICLINIC

        #Implicit rhombohedral...
        if ([a,b,c].count(None) == 2 and angles.count(None) == 2):
            #Set all angles to be same
            for val in angles:
                if val is not None:
                    angles = [val,]*3
            for val in [a,b,c]:
                if val is not None:
                    a = val
                    b = val
                    c = val
                    break
            lat_sys = CrystalSystem.RHOMBOHEDRAL
        elif (a == b and a == c and al == be and al == ga):
            #Explicit rhombohedral
            lat_sys = CrystalSystem.RHOMBOHEDRAL
            pass
        else:
            #We've got some 90degree angles or all angles have been specified
            #Set any angles with value None to 90 & find out which angles are identical
            for i in range(3):
                if angles[i] is None:
                    angles[i] = 90
            angles_compared = [angles[0]==angles[1],
                               angles[0]==angles[2],
                               angles[1]==angles[2]]
            
            #Analyse the user input
            if False not in angles_compared:
                if (a == b and a == c) or (b is None and c is None):
                    b = c = a
                    if angles[0] != 90:
                        #Rhombohedral
                        lat_sys = CrystalSystem.RHOMBOHEDRAL
                        pass
                    else:
                        #Cubic
                        lat_sys = CrystalSystem.CUBIC
                        pass
                elif (a == b and a != c) or (b is None and c is not None):
                    #Tetragonal
                    b = a
                    p_axis = PrincipleAxis.C
                    lat_sys = CrystalSystem.TETRAGONAL
                elif (a != b and c is None):
                    #Tetragonal (b -> c)
                    c = b
                    b = a
                    p_axis = PrincipleAxis.C
                    lat_sys = CrystalSystem.TETRAGONAL
                else:
                    #Orthorhombic (a != b != c)
                    lat_sys = CrystalSystem.ORTHORHOMBIC
                    if (lattice.b is None or lattice.c is None):
                        raise LatticeException("Orthorhombic requires all three lengths to be given")
            elif True not in angles_compared:
                #Triclinic
                if (b is None or c is None) and None not in angles:
                    raise LatticeException("Triclinic requires three lengths three angles to be given")
            else:
                if (120 in angles and angles_compared.count(True) == 1):
                    #Hexagonal
                    if (b is None and c is not None):
                        b = a
                    elif (b is not None and c is None):
                        c = b
                        b = a
                    else:
                        raise LatticeException("Hexagonal requires at least two lengths to be given")
                    angles = [90,90,120]
                    p_axis = PrincipleAxis.C
                    lat_sys = CrystalSystem.HEXAGONAL
                else:
                    #Monoclinic
                    if lattice.principle_axis != PrincipleAxis.NONE:
                        p_axis = lattice.principle_axis
                        old_angles = angles
                        angles = [90]*3

                        if angles_compared[0]:
                            #ga different
                            angles[p_axis.value] = old_angles[2]
                        elif angles_compared[1]:
                            #be different
                            angles[p_axis.value] = old_angles[1]
                        else:
                            #al different
                            angles[p_axis.value] = old_angles[0]
                    else:
                        if (lattice.b is None or lattice.c is None):
                            raise LatticeException("Monoclinic requires all three lengths to be given")
                        if angles_compared[0]:
                            #ga different
                            p_axis = PrincipleAxis.C
                        elif angles_compared[1]:
                            #be different
                            p_axis = PrincipleAxis.B
                        else:
                            #al different
                            p_axis = PrincipleAxis.A
                    lat_sys = CrystalSystem.MONOCLINIC

        return Lattice(a, b, c, angles[0], angles[1], angles[2], p_axis, lat_sys)

    def __determine_metric_tensor(self):
        def __calc_offaxis(a, b, angle):
            result = a * b * np.cos(angle)
            if abs(result) < 1e-10:
                return 0.0
            return result

        #Generate metric tensor and store the.lattice space lattice
        p00 = self.lattice.a*self.lattice.a
        p11 = self.lattice.b*self.lattice.b
        p22 = self.lattice.c*self.lattice.c
        #Off axis
        p01 = __calc_offaxis(self.lattice.a, self.lattice.b, self.ga_r)
        p02 = __calc_offaxis(self.lattice.a, self.lattice.c, self.be_r)
        p12 = __calc_offaxis(self.lattice.b, self.lattice.c, self.al_r)

        return np.matrix([
            [p00, p01, p02],
            [p01, p11, p12],
            [p02, p12, p22]])

    def __determine_reciprocal_lattice(self):
        #Reciprocal_lattice lattice lengths
        r_a, r_b, r_c = map(lambda x: np.sqrt(self.reciprocal_metric_tensor[x,x]), [0, 1, 2])

        #Reciprocal_lattice lattice angles
        r_al = np.degrees(np.arccos(self.reciprocal_metric_tensor[1,2] / (r_b * r_c)))
        r_be = np.degrees(np.arccos(self.reciprocal_metric_tensor[0,2] / (r_a * r_c)))
        r_ga = np.degrees(np.arccos(self.reciprocal_metric_tensor[0,1] / (r_a * r_b)))
        
        return Lattice(r_a, r_b, r_c, r_al, r_be, r_ga)

    def __convert_to_vector(self, vector):
        #Check we've been given a lattice vector
        vector = np.matrix(vector)
        if vector.shape != (3L,1L,):
            if vector.size == 3:
                vector = vector.reshape(3,1)
            else:
                raise MatrixException("Given vector cannot be formed into a 3x1 vector.")
        return vector

    def lattice(self):
        return self.lattice

    def find_vector_magnitude(self, vector, tensor=None):
        #Default to use the metric_tensor
        if (tensor is None):
            tensor=self.metric_tensor

        #convert our input to a column vector
        vector = self.__convert_to_vector(vector)

        #Use a given metric tensor (tensor) to calculate vector length: |r| = sqrt(r^T . G . r)
        fix_zeros = np.vectorize(lambda x: 0 if (abs(x) <= 1e-10) else x)
        fixed_product = fix_zeros(vector.T.dot(tensor).dot(vector))
        return np.linalg.norm(np.sqrt(fixed_product))

    def find_plane_dspacing(self, plane_indices):
        one_on_d = self.find_vector_magnitude(plane_indices, self.reciprocal_metric_tensor)
        return 1/one_on_d

    def find_cartesian_coordinates(self, crystal_coordinates, crystal_origin=[0,0,0]):
        return self.__find_coordinates(crystal_coordinates, crystal_origin, self.orthonormalisation_matrix)

    def find_crystal_coordinates(self, cartesian_coordinates, cartesian_origin=[0,0,0]):
        return self.__find_coordinates(cartesian_coordinates, cartesian_origin, np.linalg.inv(self.orthonormalisation_matrix))

    def __find_coordinates(self, coordinates, origin, matrix):
        #Ensure we have two column vectors as input
        coordinates = self.__convert_to_vector(coordinates)
        origin = self.__convert_to_vector(origin)
        coordinates = coordinates - origin

        return matrix.dot(coordinates)

    def __get_crystal_family_for_lattice_system(self):
        crystal_family_to_lattice_system = {
            CrystalSystem.TRICLINIC : CrystalSystem.TRICLINIC,
            CrystalSystem.MONOCLINIC : CrystalSystem.MONOCLINIC,
            CrystalSystem.ORTHORHOMBIC : CrystalSystem.ORTHORHOMBIC,
            CrystalSystem.TETRAGONAL : CrystalSystem.TETRAGONAL,
            CrystalSystem.RHOMBOHEDRAL : CrystalSystem.HEXAGONAL,
            CrystalSystem.HEXAGONAL : CrystalSystem.HEXAGONAL,
            CrystalSystem.CUBIC : CrystalSystem.CUBIC
            }
        return crystal_family_to_lattice_system[self.lattice.lattice_system]


class PrincipleAxis(Enum):
    __order__ = "A B C" #Needed for python 2.7
    A = 0
    B = 1 #Standard setting (be != 90)
    C = 2
    NONE = -1

class CrystalSystem(Enum):
    __order__ = "TRICLINIC MONOCLINIC ORTHORHOMBIC TETRAGONAL TRIGONAL RHOMBOHEDRAL HEXAGONAL CUBIC"
    TRICLINIC = 0
    MONOCLINIC = 1
    ORTHORHOMBIC = 2
    TETRAGONAL = 3
    TRIGONAL = 4
    RHOMBOHEDRAL = 5
    HEXAGONAL = 6
    CUBIC = 7


Lattice = namedtuple('Lattice', 'a, b, c, al, be, ga, principle_axis, lattice_system')
Lattice.__new__.__defaults__ = (None, None, None, None, None, PrincipleAxis.NONE, CrystalSystem.TRICLINIC)