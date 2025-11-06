"""
    Algorithmic verification of topological MBQC circuits.
    Copyright (C) 2025  Marvin Schwiering

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

###############
### IMPORTS ###
###############

# Typing
from __future__ import annotations
from typing import Any, Type, Literal, Callable, Tuple, overload
from numbers import Real
from enum import Enum

# General
from warnings import warn
from pathlib import Path
import json

# Parallel Computing
from multiprocessing.pool import Pool
from multiprocessing import Pool as ProcessPool
from multiprocessing.dummy import Pool as ThreadPool

# Data Structures and Helful Subroutines
from itertools import product, chain
from functools import reduce
from bidict import bidict

# Math
import numpy as np
import galois as gl
from galois import GF2, FieldArray


##############
### TYPING ###
##############

# Global Toggle for Type Checking
safe_mode_global:bool = False

# Define Order and a Function to Check
class Order(Enum):
    SITE = 0
    LINK = 1
    FACE = 2
    CUBE = 3

def is_order(order: Any) -> bool:
    return isinstance(order, Order)

# Define Mode and a Function to Check
class Mode(Enum):
    PRIMAL = 0
    PRIM = 0
    DUAL = 1

def is_mode(mode: Any) -> bool:
    return isinstance(mode, Mode)

# Define Shape and a Function to Check
type Shape = Tuple[int, int, int]

def is_shape(shape: Any) -> bool:
    # Shape is a tuple ...
    if not type(shape) is tuple:
        return False
    # ... with three components ...
    if not len(shape) == 3:
        return False
    # ... that are all ints ...
    if not all([type(i) is int for i in shape]):
        return False
    # ... that are greater than one
    if not all([i > 1 for i in shape]):
        return False

    return True

# Define Coordinates and a Function to Check
type Coordinates = Tuple[Real, Real, Real]

def is_coordinates(coordinates: Any, shape:Shape = None, lattice:Lattice = None) -> bool:
    # Shape XOR Lattice must be non-trivial
    assert (shape is None and type(lattice) is Lattice) != (lattice is None and is_shape(shape))
    shape = lattice.shape if shape is None else shape

    # Coordinates is a tuple ...
    if not type(coordinates) is tuple:
        return False
    # ... with three components ...
    if not len(coordinates) == 3:
        return False
    # ... that are all ints or floats ...
    if not all([type(i) in {int, float} for i in coordinates]):
        return False
    # ... that are non-negative ...
    if not all([i >= 0 for i in coordinates]):
        return False
    # ... and divisible by 1/2 ...
    if not all(i % 0.5 == 0 for i in coordinates):
        return False
    # ... and contained within shape (if given).
    if not all([c <= s-1 for c, s in zip(coordinates, shape)]):
        return False

    return True

# Function to Determine Whether a Structure is At The Edge of the Lattice
def is_at_edge(coordinates:Coordinates, shape:Shape = None, lattice:Lattice = None) -> bool:
    # Shape XOR Lattice must be non-trivial
    assert (shape is None and type(lattice) is Lattice) != (lattice is None and is_shape(shape))
    shape = lattice.shape if shape is None else shape

    return sum([int(coordinates[i] in {0, shape[i] - 1}) for i in range(3)]) > 0

########################
### HELPER-FUNCTIONS ###
########################

def dimensions(shape: Shape, order:Order = None, safe_mode:bool = safe_mode_global) -> int | tuple[int, int, int, int]:
    # Input- and Type-Checking
    if safe_mode:
        assert is_shape(shape), f"{shape} is not a Shape."
        assert is_order(order) or order is None, f"{order} is not a valid chain complex order."

    # Calculate Dimensions for All Orders
    dimensions_list = [shape[0] * shape[1] * shape[2], # Number of Sites
                       (shape[0] - 1) * shape[1] * shape[2] + shape[0] * (shape[1] - 1) * shape[2] + shape[0] * shape[1] * (shape[2] - 1), # Number of Links
                       (shape[0] - 1) * (shape[1] - 1) * shape[2] + (shape[0] - 1) * shape[1] * (shape[2] - 1) + shape[0] * (shape[1] - 1) * (shape[2] - 1), # Number of Faces
                       (shape[0] - 1) * (shape[1] - 1) * (shape[2] - 1) # Number of Cubes
                       ]

    return dimensions_list if order is None else dimensions_list[order.value]

def determine_order(coordinates:Coordinates, shape:Shape = None, lattice:Lattice = None, safe_mode:bool = safe_mode_global) -> Order:
    if safe_mode:
        if shape is None and lattice is None:
            warn("determine_order could not run in safe mode due to insufficient parameters.")
        else:
            assert is_coordinates(coordinates, shape = shape, lattice = lattice)

    return Order( sum([int(i % 1 == 0.5) for i in coordinates]) )


###############
### CLASSES ###
###############

class Chain:

    def __init__(self, lattice:Lattice, safe_mode:bool = safe_mode_global) -> None:
        # Input- and Type-Checking
        if safe_mode:
            assert type(lattice) is Lattice

        # Elementary Properties
        self._lattice:Lattice = lattice
        self._order:Order = None
        self._tuples:set = None
        self._vector: FieldArray = None

    def __radd__(self, other:Chain, safe_mode:bool = safe_mode_global) -> Chain:
        # TODO: This Is An Ugly Work-Around to Make sum(chain_list) Work If chain_list:list[Chain]. Refactor?!
        if other == 0:
            if safe_mode:
                warn(f"Addition 0 + Chain returns the Chain Object.")
            return self
        else:
            return self.__add__(other = other, safe_mode = safe_mode)

    def __add__(self, other:Chain, safe_mode:bool = safe_mode_global) -> Chain:
        # Input- and Type-Checking
        if safe_mode:
            assert type(other) is Chain
            assert self._lattice == other._lattice
            assert self.order == other.order
            assert self._vector.shape == other._vector.shape

        # Calculate Summand
        result = Chain(lattice = self._lattice)
        result._order = self.order
        result._tuples = self._tuples ^ other._tuples # Set Symmetric Difference
        result._vector = self._vector + other._vector

        return result

    def __sub__(self, other:Chain, safe_mode:bool = safe_mode_global) -> Chain:
        # With Addition Modulo 2, Every Integer/Vector Is Its Own Additive Inverse.
        if safe_mode:
            warn("Subtraction of Chains is Equivalent to Their Addition; Any Reason You Are Using the Former?")
        return self.__add__(other = other, safe_mode = safe_mode)

    def __mul__(self, other:Chain, safe_mode:bool = safe_mode_global) -> Chain:
        # TODO: This Could be Implemented, But Will It Actually Be Used Anywhere? I Doubt It.
        raise NotImplementedError()

    def setminus(self, to_exclude:Chain, safe_mode:bool = safe_mode_global) -> Chain:
        # Returns the Chain Corresponding to the Set Difference of self._tuples - to_exclude._tuples.

        # Input- and Type-Checking
        if safe_mode:
            assert type(to_exclude) is Chain
            assert self._lattice == to_exclude._lattice
            assert self.order == to_exclude.order
            assert self._vector.shape == to_exclude._vector.shape

        # Calculate Set Difference
        result = Chain(lattice = self._lattice)
        result._order = self._order
        result._tuples = self._tuples - to_exclude._tuples
        result._vector = self._vector.copy()
        result._vector[ np.where(to_exclude == 1) ] = 0

        return result

    @overload
    def boundary(self, mode:Mode, relative_boundary:bool = False, smoothen_dual_lattice:bool = False, all_limits:Literal[False] = False, safe_mode:bool = safe_mode_global) -> Chain: ...

    @overload
    def boundary(self, mode:Mode, relative_boundary:bool = False, smoothen_dual_lattice:bool = False, all_limits:Literal[True] = True, safe_mode:bool = safe_mode_global) -> set[tuple[Real, Real, Real]]: ...

    # Real Function with Actual Implementation
    def boundary(self, mode:Mode, relative_boundary:bool = False, smoothen_dual_lattice:bool = False, all_limits:bool = False, safe_mode:bool = safe_mode_global) -> Chain | set[tuple[Real, Real, Real]]:
        # Input- and Type-Checking
        if safe_mode:
            assert is_mode(mode)
            assert type(relative_boundary) is bool
            assert type(smoothen_dual_lattice) is bool

        # Functions to Calculate the Primal Boundary of a Site/Link/Face/Cube From Its Coordinates
        def boundary_prim_site(coordinates:Coordinates) -> set[Coordinates]:
            raise RuntimeError("0-Chains Have No Primal Boundary.")

        def boundary_prim_link(coordinates:Coordinates) -> set[Coordinates]:
            if safe_mode:
                assert is_coordinates(coordinates, lattice = self._lattice)
                assert determine_order(coordinates, safe_mode = False) == Order.LINK

            differing_index = [c % 1 == 0.5 for c in coordinates].index(True)
            if differing_index == 0:
                tuples = { (coordinates[0] - 0.5, coordinates[1], coordinates[2]), (coordinates[0] + 0.5, coordinates[1], coordinates[2]) }
            elif differing_index == 1:
                tuples = { (coordinates[0], coordinates[1] - 0.5, coordinates[2]), (coordinates[0], coordinates[1] + 0.5, coordinates[2]) }
            elif differing_index == 2:
                tuples = { (coordinates[0], coordinates[1], coordinates[2] - 0.5), (coordinates[0], coordinates[1], coordinates[2] + 0.5) }
            else:
                raise RuntimeError(f"Coordinates {coordinates} yielded a differing index of {differing_index}; invalid.")

            return tuples

        def boundary_prim_face(coordinates:Coordinates) -> set[Coordinates]:
            if safe_mode:
                assert is_coordinates(coordinates, lattice = self._lattice)
                assert determine_order(coordinates, safe_mode = False) == Order.FACE

            differing_index = [c % 1 == 0 for c in coordinates].index(True)
            if differing_index == 0:
                tuples = { (coordinates[0], coordinates[1] - 0.5, coordinates[2]),
                           (coordinates[0], coordinates[1] + 0.5, coordinates[2]),
                           (coordinates[0], coordinates[1], coordinates[2] - 0.5),
                           (coordinates[0], coordinates[1], coordinates[2] + 0.5) }
            elif differing_index == 1:
                tuples = { (coordinates[0] - 0.5, coordinates[1], coordinates[2]),
                           (coordinates[0] + 0.5, coordinates[1], coordinates[2]),
                           (coordinates[0], coordinates[1], coordinates[2] - 0.5),
                           (coordinates[0], coordinates[1], coordinates[2] + 0.5) }
            elif differing_index == 2:
                tuples = { (coordinates[0] - 0.5, coordinates[1], coordinates[2]),
                           (coordinates[0] + 0.5, coordinates[1], coordinates[2]),
                           (coordinates[0], coordinates[1] - 0.5, coordinates[2]),
                           (coordinates[0], coordinates[1] + 0.5, coordinates[2]) }
            else:
                raise RuntimeError(f"Coordinates {coordinates} yielded a differing index of {differing_index}; invalid.")

            return tuples

        def boundary_prim_cube(coordinates:Coordinates) -> set[Coordinates]:
            if safe_mode:
                assert is_coordinates(coordinates, lattice = self._lattice)
                assert determine_order(coordinates, safe_mode = False) == Order.CUBE

            tuples = { (coordinates[0] - 0.5, coordinates[1], coordinates[2]),
                       (coordinates[0] + 0.5, coordinates[1], coordinates[2]),
                       (coordinates[0], coordinates[1] - 0.5, coordinates[2]),
                       (coordinates[0], coordinates[1] + 0.5, coordinates[2]),
                       (coordinates[0], coordinates[1], coordinates[2] - 0.5),
                       (coordinates[0], coordinates[1], coordinates[2] + 0.5) }

            return tuples


        # Functions to Calculate the Dual Boundary of a Site/Link/Face/Cube From Its Coordinates
        def boundary_dual_site(coordinates:Coordinates) -> set[Coordinates] | set[tuple[Real, Real, Real]]:
            if safe_mode:
                assert is_coordinates(coordinates, lattice = self._lattice)
                assert determine_order(coordinates, safe_mode = False) == Order.SITE

            tuples = { (coordinates[0] - 0.5, coordinates[1], coordinates[2]),
                       (coordinates[0] + 0.5, coordinates[1], coordinates[2]),
                       (coordinates[0], coordinates[1] - 0.5, coordinates[2]),
                       (coordinates[0], coordinates[1] + 0.5, coordinates[2]),
                       (coordinates[0], coordinates[1], coordinates[2] - 0.5),
                       (coordinates[0], coordinates[1], coordinates[2] + 0.5) }

            # Filter Out Coordinates Outside the Lattice (If Desired)
            if not all_limits:
                tuples = {t for t in tuples if is_coordinates(t, lattice = self._lattice)}

            return tuples

        def boundary_dual_link(coordinates:Coordinates) -> set[Coordinates] | set[tuple[Real, Real, Real]]:
            if safe_mode:
                assert is_coordinates(coordinates, lattice = self._lattice)
                assert determine_order(coordinates, safe_mode = False) == Order.LINK

            differing_index = [c % 1 == 0.5 for c in coordinates].index(True)
            if differing_index == 0:
                tuples = { (coordinates[0], coordinates[1] - 0.5, coordinates[2]),
                           (coordinates[0], coordinates[1] + 0.5, coordinates[2]),
                           (coordinates[0], coordinates[1], coordinates[2] - 0.5),
                           (coordinates[0], coordinates[1], coordinates[2] + 0.5) }
            elif differing_index == 1:
                tuples = { (coordinates[0] - 0.5, coordinates[1], coordinates[2]),
                           (coordinates[0] + 0.5, coordinates[1], coordinates[2]),
                           (coordinates[0], coordinates[1], coordinates[2] - 0.5),
                           (coordinates[0], coordinates[1], coordinates[2] + 0.5) }
            elif differing_index == 2:
                tuples = { (coordinates[0] - 0.5, coordinates[1], coordinates[2]),
                           (coordinates[0] + 0.5, coordinates[1], coordinates[2]),
                           (coordinates[0], coordinates[1] - 0.5, coordinates[2]),
                           (coordinates[0], coordinates[1] + 0.5, coordinates[2]) }
            else:
                raise RuntimeError(f"Coordinates {coordinates} yielded a differing index of {differing_index}; invalid.")

            # Filter Out Coordinates Outside the Lattice (If Desired)
            if not all_limits:
                tuples = {t for t in tuples if is_coordinates(t, lattice = self._lattice)}

            return tuples

        def boundary_dual_face(coordinates:Coordinates) -> set[Coordinates] | set[tuple[Real, Real, Real]]:
            if safe_mode:
                assert is_coordinates(coordinates, lattice = self._lattice)
                assert determine_order(coordinates, safe_mode = False) == Order.FACE

            differing_index = [c % 1 == 0 for c in coordinates].index(True)
            if differing_index == 0:
                tuples = { (coordinates[0] - 0.5, coordinates[1], coordinates[2]), (coordinates[0] + 0.5, coordinates[1], coordinates[2]) }
            elif differing_index == 1:
                tuples = { (coordinates[0], coordinates[1] - 0.5, coordinates[2]), (coordinates[0], coordinates[1] + 0.5, coordinates[2]) }
            elif differing_index == 2:
                tuples = { (coordinates[0], coordinates[1], coordinates[2] - 0.5), (coordinates[0], coordinates[1], coordinates[2] + 0.5) }
            else:
                raise RuntimeError(f"Coordinates {coordinates} yielded a differing index of {differing_index}; invalid.")

            # Filter Out Coordinates Outside the Lattice (If Desired)
            if not all_limits:
                tuples = {t for t in tuples if is_coordinates(t, lattice = self._lattice)}

            return tuples

        def boundary_dual_cube(coordinates:Coordinates) -> set[Coordinates] | set[tuple[Real, Real, Real]]:
            raise RuntimeError("3-Chains Have No Dual Boundary.")

        # Master Functions that Call the Relevant Prior Function,, then ...
        # 1. Remove Structure in the Relative Boundary (If Desired),
        # 1. Remove Structures at the Edge to Smoothen the Dual Lattice (If Desired),
        # and Return the Result as a Chain or Set of Limiting Points.
        order_in = self.order
        order_to = Order(self.order.value - 1) if mode == Mode.PRIM else Order(self.order.value + 1)

        methods = {(Mode.PRIM, Order.SITE): boundary_prim_site, (Mode.PRIM, Order.LINK): boundary_prim_link,
                   (Mode.PRIM, Order.FACE): boundary_prim_face, (Mode.PRIM, Order.CUBE): boundary_prim_cube,
                   (Mode.DUAL, Order.SITE): boundary_dual_site, (Mode.DUAL, Order.LINK): boundary_dual_link,
                   (Mode.DUAL, Order.FACE): boundary_dual_face, (Mode.DUAL, Order.CUBE): boundary_dual_cube}

        def boundary_as_Chain(coordinates:Coordinates) -> Chain:
            # Compute the Boundary (Excluding Limiting Points Outside the Lattice)
            tuples = methods[mode, self.order](coordinates)

            # Compute the Relative Boundary (If Desired)
            if relative_boundary and len({order_in, order_to} & {Order.LINK, Order.FACE}) > 0:
                # If the Order Mapping In Contains a Relative Boundary, Remove All Boundary-Coordinates If the In-Coordinates Are In It
                if order_in in {Order.LINK, Order.FACE}:
                    tuples = set() if coordinates in self._lattice.defects(order_in).as_tuples() else tuples
                # If the Order Mapping To Contains a Relative Boundary, Remove All Boundary-Coordinates That Are In It
                if order_to in {Order.LINK, Order.FACE}:
                    tuples = {c_boundary for c_boundary in tuples if not(c_boundary in self._lattice.defects(order_to).as_tuples())}

            # Smoothen The Dual Lattice (If Desired)
            if smoothen_dual_lattice and mode == Mode.DUAL:
                # If the Structure Represented by coordinates is at the Edge, Remove It (Yielding an Empty Boundary)
                # Else, Remove All Structures in Its Boundary that are at the Edge
                tuples = {c_boundary for c_boundary in tuples if not(is_at_edge(coordinates, lattice = self._lattice) or is_at_edge(c_boundary, lattice = self._lattice))}

            return Chain.fromCoordinates(tuples, self._lattice) if len(tuples) > 0 else Chain.createEmpty(order_to, self._lattice)

        def limits_as_set(coordinates:Coordinates) -> set[tuple[Real, Real, Real]]:
            # Compute the Boundary (Including Limiting Points Inside the Lattice)
            tuples = methods[mode, self.order](coordinates)

            return tuples

        # Call the Relevant Function for All Structures and Calculate Their Sum/Symmetric Set Difference
        if not all_limits:
            result = sum([ boundary_as_Chain(coordinates) for coordinates in self.as_tuples() ])
        else:
            result = reduce(set.symmetric_difference, [limits_as_set(coordinates) for coordinates in self.as_tuples()])

        return result

    @classmethod
    def fromCoordinates(cls, coordinates:Coordinates | set[Coordinates], lattice:Lattice, safe_mode:bool = safe_mode_global) -> Chain:
        # Input- and Type-Checking
        if safe_mode:
            assert is_coordinates(coordinates, lattice = lattice) or (type(coordinates) is set and all([is_coordinates(c, lattice = lattice) for c in coordinates]))
            assert type(lattice) is Lattice

        # Normalize Input to a Set so that Covering Sets Covers All Inputs
        if not type(coordinates) is set:
            coordinates = {coordinates}

        # Create the Chain and Add All Tuples
        chain = Chain(lattice = lattice)
        for c in coordinates:
            chain.append(c)
        return chain

    @classmethod
    def fromVector(cls, vector:FieldArray, lattice:Lattice, safe_mode:bool = safe_mode_global) -> Chain:
        # Input- and Type-Checking
        if safe_mode:
            assert isinstance(vector, FieldArray)
            assert type(lattice) is Lattice

        # Create the Chain and Add the Vector
        chain = Chain(lattice = lattice)
        chain.append(vector)
        return chain

    @classmethod
    def createEmpty(cls, order:Order, lattice:Lattice, safe_mode:bool = safe_mode_global) -> Chain:
        # Input- and Type-Checking
        if safe_mode:
            assert is_order(order)
            assert type(lattice) is Lattice

        chain = Chain(lattice = lattice)
        chain._order = order
        chain._tuples = {}
        chain._vector = GF2.Zeros(shape = dimensions(shape = lattice.shape, order = order))

        return chain

    def append(self, to_append:Coordinates|FieldArray, safe_mode:bool = safe_mode_global) -> None:
        if safe_mode:
            assert is_coordinates(to_append, lattice = self._lattice) or isinstance(to_append, FieldArray)

        def appendCoordinates(coordinates:Coordinates) -> None:
            # Check if Chain Was Empty Until Now
            if self._order is None:
                self._order = determine_order(coordinates)
                self._tuples = set()
                self._vector = GF2.Zeros(shape = dimensions(shape = self._lattice.shape, order = self.order))

            # Check If the Structure is Already in the Chain
            if coordinates in self._tuples:
                raise ValueError(f"{coordinates} already part of chain {self}.")

            # Check If the Structure has the Same Order as the Chain
            if not self.order == (coordinates_order := determine_order(coordinates)):
                raise ValueError(f"{coordinates} is of order {coordinates_order}, but this chain is of order {self.order}.")

            # Adjust Tuples and Vector
            self._tuples.add(coordinates)
            self._vector[self._lattice.bijection(self.order)[coordinates]] = 1

        def appendVector(vector:FieldArray) -> None:
            # Check if Chain Was Empty Until Now
            if self._order is None:
                self._order = Order( dimensions(shape = self._lattice.shape).index( vector.shape[0] ) )
                self._tuples = set()
                self._vector = GF2.Zeros(shape = vector.shape)

            # Check If the Structure is Already in the Chain
            if not set(self._vector.nonzero()[0]).isdisjoint( set(vector.nonzero()[0]) ):
                raise ValueError(f"{vector} is already (partially) part of chain {self}.")

            # Check If the Vector has the Same Dimension as the Chain Vector
            if not self._vector.shape == vector.shape:
                raise ValueError(f"{vector} has shape {vector.shape}, but this chain has dimension {self._vector.shape}.")

            # Adjust Tuples and Vector
            self._tuples.update([self._lattice.bijection(self._order).inverse[index] for index in vector.nonzero()[0]])
            self._vector += vector

        if isinstance(to_append, FieldArray):
            appendVector(to_append)
        else:
            appendCoordinates(to_append)

    @property
    def order(self) -> Order:
        if self._order is None:
            raise AttributeError("Empty Chain has no order.")

        return self._order

    @property
    def state_id(self) -> tuple[Real]:
        # Why Is This Always Recalculated and Not Stored? So That Not All Operations Changing the Chain Must Remember to Reset This
        return tuple(chain.from_iterable( [self._lattice.shape] + sorted(self.as_tuples())))

    def as_tuples(self) -> set[Coordinates]:
        return self._tuples

    def as_vector(self) -> FieldArray:
        return self._vector

    def as_column_vector(self) -> FieldArray:
        return np.reshape(self._vector, newshape = (self._vector.shape[0], 1))


class Lattice:

    def __init__(self, shape:Shape, safe_mode:bool = safe_mode_global) -> None:
        if safe_mode:
            assert is_shape(shape)

        # Elementary Properties
        self._shape:Shape = shape
        self._defects_prim:Chain = Chain.createEmpty(order = Order.LINK, lattice = self)
        self._defects_dual:Chain = Chain.createEmpty(order = Order.FACE, lattice = self)
        self._partial:dict[tuple[Mode, Order, bool], FieldArray] = {(mode, order, relative_boundary): None for relative_boundary in {False, True} for order in Order for mode in Mode}
        self._bijections:dict[Order, bidict] = {order: None for order in Order}

        # Target Vectors
        self._target_vectors_prim: dict[str, Chain] = dict()
        self._target_vectors_dual: dict[str, Chain] = dict()

        self._verify_prim_each:dict[str, bool] = dict()
        self._verify_dual_each:dict[str, bool] = dict()

        self._verify_prim_bulk:bool = None
        self._verify_dual_bulk:bool = None

        # Correlation Surfaces
        self._surfaces_prim:dict[str, Chain] = dict()
        self._surfaces_dual:dict[str, Chain] = dict()

        self._optimality_prim:dict[str, bool] = dict()
        self._optimality_dual:dict[str, bool] = dict()

    def populate(self, defects_prim:Chain = None, defects_dual:Chain = None, target_vectors:dict[str, Chain] = None, surfaces:dict[str, Chain] = None, optimality:dict[str, Chain] = None, safe_mode:bool = safe_mode_global) -> None:
        # PRIMAL DEFECTS
        if not(defects_prim is None):
            # Check Input
            if safe_mode:
                assert type(defects_prim) is Chain and defects_prim.order == Order.LINK
                assert len(self._defects_prim.as_tuples()) == 0

            self._defects_prim = defects_prim

        # DUAL DEFECTS
        if not(defects_dual is None):
            # Check Input
            if safe_mode:
                assert type(defects_dual) is Chain and defects_dual.order == Order.FACE
                assert len(self._defects_dual.as_tuples()) == 0

            self._defects_dual = defects_dual

        # TARGET VECTORS
        if not(target_vectors is None):
            # Check Input
            if safe_mode:
                assert type(target_vectors) is dict and all([type(tv) is Chain and tv.order in {Order.LINK, Order.FACE} for tv in target_vectors.values()])
                assert len(self._target_vectors_prim) == 0 and len(self._target_vectors_dual) == 0

            self._target_vectors_prim = {name:target_vector for name, target_vector in target_vectors.items() if target_vector.order == Order.LINK}
            self._target_vectors_dual = {name:target_vector for name, target_vector in target_vectors.items() if target_vector.order == Order.FACE}

            self._verify_prim_each = {name:None for name in self._target_vectors_prim.keys()}
            self._verify_dual_each = {name:None for name in self._target_vectors_dual.keys()}

        # SURFACES
        if not(surfaces is None):
            # Check Input
            if safe_mode:
                assert type(surfaces) is dict and all([type(s) is Chain and s.order in {Order.LINK, Order.FACE} for s in surfaces.values()])
                assert len(self._surfaces_prim) == 0 and len(self._surfaces_dual) == 0
                assert len(surfaces) == len(self._target_vectors_prim) + len(self._target_vectors_dual)
                assert set(surfaces.keys()) == set(self._target_vectors_prim.keys()) | set(self._target_vectors_dual.keys())
                # TODO: Check Here If They Are Actually Solutions to A x = b?!

            self._surfaces_prim = {name:surface for name, surface in surfaces.items() if surface.order == Order.FACE}
            self._surfaces_dual = {name:surface for name, surface in surfaces.items() if surface.order == Order.LINK}

        elif not(target_vectors) is None:
            # If Target Vectors Were Initialized, But Not the Surfaces, We Set the Surfaces to Empty
            self._surfaces_prim = {name:None for name in self._target_vectors_prim.keys()}
            self._surfaces_dual = {name:None for name in self._target_vectors_dual.keys()}

        # OPTIMALITY
        if not(optimality is None):
            # Check Input
            if safe_mode:
                assert type(optimality) is dict and all([type(o) is bool for o in optimality.values()])
                assert len(self._optimality_prim) == 0 and len(self._optimality_dual) == 0
                assert len(optimality) == len(self._surfaces_prim) + len(self._surfaces_dual)
                assert set(optimality.keys()) == set(self._surfaces_prim.keys()) | set(self._surfaces_dual.keys())

            self._optimality_prim = {name:optimal for name, optimal in optimality.items() if name in self._surfaces_prim.keys()}
            self._optimality_dual = {name:optimal for name, optimal in optimality.items() if name in self._surfaces_dual.keys()}

        elif not(target_vectors) is None:
            # If Target Vectors Were Initialized, But Not the Optimality, We Set the Optimality to None
            self._optimality_prim = {name:None for name in self._surfaces_prim.keys()}
            self._optimality_dual = {name:None for name in self._surfaces_dual.keys()}


    @classmethod
    def load(cls, save_name:str = None, path:Path = None, safe_mode:bool = True) -> Lattice:
        """Note: load Will Raise an Error If Not Executed in safe_mode."""

        # Check Input
        if safe_mode:
            assert (type(save_name) is str and path is None) != (save_name is None and isinstance(path, Path))

        # ACQUIRE RAW DATA
        try:
            # Load the JSON Into the Arguments Dict
            path = path if save_name is None else Path(f"saves/{save_name}.json")
            with open(path) as f:
                arguments = json.load(f)

            # Extract the Data From the Arguments Dict
            shape = tuple(arguments["shape"])
            defects_prim_tuples = {tuple(coordinates) for coordinates in arguments["defects_prim"]}
            defects_dual_tuples = {tuple(coordinates) for coordinates in arguments["defects_dual"]}
            target_vectors_tuples = {name:{tuple(coordinates) for coordinates in target_vector} for name, target_vector in arguments["target_vectors"].items()}
            surfaces_tuples = {name: ({tuple(coordinates) for coordinates in target_vector} if not(target_vector is None) else None) for name, target_vector in arguments["surfaces"].items() }
            optimality = {name:optimal for name, optimal in arguments["optimality"].items()}
        except IOError as e:
            warn(f"Could Not Load Lattice File; Got Exception {e}")
            return None

        except KeyError as e:
            warn(f"Could Not Load Lattice Due to Incomplete Data; Got Exception {e}")
            return None

        # VALIDATE THE DATA
        try:
            if safe_mode:
                # VALIDATE THE SHAPE
                assert is_shape(shape)

                # VALIDATE THE DEFECTS
                # For All Coordinates in defects_prim, Assert That They Are Indeed Coordinates of Order Faces
                for coordinates in defects_prim_tuples:
                    assert is_coordinates(coordinates = coordinates, shape = shape)
                    assert determine_order(coordinates = coordinates, shape = shape) == Order.LINK
                # For All Coordinates in defects_dual, Assert That They Are Indeed Coordinates of Order Faces
                for coordinates in defects_dual_tuples:
                    assert is_coordinates(coordinates = coordinates, shape = shape)
                    assert determine_order(coordinates = coordinates, shape = shape) == Order.FACE

                # VALIDATE THE TARGET VECTORS
                # For All Coordinates in target_vector, Assert That They Are Indeed Coordinates of the Same Order ∈ {Links, Faces}
                for name, target_vector in target_vectors_tuples.items():
                    reference_order = determine_order(coordinates = next(iter(target_vector)), shape = shape)
                    assert reference_order in {Order.LINK, Order.FACE}
                    for coordinates in target_vector:
                        assert is_coordinates(coordinates = coordinates, shape = shape)
                        assert determine_order(coordinates = coordinates, shape = shape) == reference_order

                # VALIDATE THE SURFACES
                # Surfaces Must Have the Same Names as the Target Vectors
                assert set(surfaces_tuples.keys()) == set(target_vectors_tuples.keys())
                # Surfaces Are Allowed to Be Unspecified
                # (It's Assumed Either All or None are Given; If Only Some are Specified, None Will be Loaded.)
                surfaces_unspecified = any([surface is None for surface in surfaces_tuples.values()])
                if not surfaces_unspecified:
                    # For All Coordinates in surface, Assert That They Are Indeed Coordinates of the Same Order ∈ {Links, Faces}
                    for name, surface in surfaces_tuples.items():
                        reference_order = determine_order(coordinates = next(iter(surface)), shape = shape)
                        assert reference_order in {Order.LINK, Order.FACE}
                        for coordinates in surface:
                            assert is_coordinates(coordinates = coordinates, shape = shape)
                            assert determine_order(coordinates = coordinates, shape = shape) == reference_order

                # VALIDATE THE OPTIMALITY
                # Optimals Must Have the Same Names as the Surfaces
                assert set(optimality.keys()) == set(surfaces_tuples.keys())
                # Optimals Are Allowed to Be Unspecified
                optimals_unspecified = all([o is None for o in optimality.values()])
                if not optimals_unspecified:
                    assert all([o in {True, False} for o in optimality.values()])

        except AssertionError as e:
            warn(f"Could Not Load Lattice Due to Failed Data Validation; Got Exception {e}")
            return None

        # CREATE AND RETURN THE LATTICE
        lattice = Lattice(shape = shape)
        lattice.populate(defects_prim = Chain.fromCoordinates(coordinates = defects_prim_tuples, lattice = lattice) if len(defects_prim_tuples) > 0 else None,
                         defects_dual = Chain.fromCoordinates(coordinates = defects_dual_tuples, lattice = lattice) if len(defects_dual_tuples) > 0 else None,
                         target_vectors = {name: Chain.fromCoordinates(coordinates = coordinates, lattice = lattice) for name, coordinates in target_vectors_tuples.items()},
                         surfaces = None if surfaces_unspecified else {name: Chain.fromCoordinates(coordinates = coordinates, lattice = lattice) for name, coordinates in surfaces_tuples.items()},
                         optimality = None if optimals_unspecified else optimality)
        return lattice

    def save(self, save_name:str = None, path:Path = None, description:str = None, safe_mode:bool = safe_mode_global) -> bool:
        # Check Input
        if safe_mode:
            assert (type(save_name) is str and path is None) != (save_name is None and isinstance(path, Path))

        # Create the JSON Object in Python
        json_object = json.dumps(indent = 4, obj = {
            "description":      description,
            "shape":            self.shape,
            "defects_prim":     sorted(self._defects_prim.as_tuples()) if type(self._defects_dual) is Chain else None,
            "defects_dual":     sorted(self._defects_dual.as_tuples()) if type(self._defects_dual) is Chain else None,
            "target_vectors":   {name: (sorted(target_vector.as_tuples()) if type(target_vector) is Chain else None) for name, target_vector in (self._target_vectors_prim | self._target_vectors_dual).items()},
            "surfaces":         {name: (sorted(surface.as_tuples()) if type(surface) is Chain else None) for name, surface in (self._surfaces_prim | self._surfaces_dual).items()},
            "optimality":       {name: optimal for name, optimal in (self._optimality_prim | self._optimality_dual).items()}
        })

        # Write to the JSON File
        try:
            with open(path, "w", newline="\n") as f:
                f.write(json_object)
            return True

        except IOError as e:
            warn(f"Could Not Save Lattice; Got Exception {e}")
            return False

    @property
    def shape(self) -> Shape:
        return self._shape

    def defects(self, order:Order, safe_mode:bool = safe_mode_global) -> Chain:
        # Check Input
        if safe_mode:
            assert is_order(order)

        if order == Order.LINK:
            return self._defects_prim
        elif order == Order.FACE:
            return self._defects_dual
        else:
            raise ValueError(f"Lattice Can Only Have Primal Defects (Order 1) or Dual Defects (Order 2), Not Defects of Order {order}.")

    def verification(self, return_each:bool = False, method:int = None, safe_mode:bool = safe_mode_global) -> bool | dict[str, bool]:
        """Note for Testing: Setting return_each to False Ensures that the Desired Method is the (Only) One Executed."""

        # Functions that Perform the Verification
        def verification_each_sequential() -> tuple[dict[str, bool], dict[str, bool]]:
            partial_prim = self.partial(mode = Mode.PRIM, order = Order.FACE, relative_boundary = True, smoothen_dual_lattice = True)
            partial_dual = self.partial(mode = Mode.DUAL, order = Order.LINK, relative_boundary = True, smoothen_dual_lattice = True)

            rank_prim = np.linalg.matrix_rank(partial_prim) # This Is ´Intercepted´ By Galois to Return the Rank over GF2
            rank_dual = np.linalg.matrix_rank(partial_dual) # This Is ´Intercepted´ By Galois to Return the Rank over GF2

            verify_prim_each = {name: rank_prim == np.linalg.matrix_rank(np.hstack([partial_prim, target_vector.as_column_vector()])) for name, target_vector in self._target_vectors_prim.items() }
            verify_dual_each = {name: rank_dual == np.linalg.matrix_rank(np.hstack([partial_dual, target_vector.as_column_vector()])) for name, target_vector in self._target_vectors_dual.items() }

            return verify_prim_each, verify_dual_each

        def verification_each_parallelized(PoolToUse:Type[Pool]) -> tuple[dict[str, bool], dict[str, bool]]:
            partial_prim = self.partial(mode = Mode.PRIM, order = Order.FACE, relative_boundary = True, smoothen_dual_lattice = True)
            partial_dual = self.partial(mode = Mode.DUAL, order = Order.LINK, relative_boundary = True, smoothen_dual_lattice = True)

            partial_prim_augmented = [np.hstack([partial_prim, target_vector.as_column_vector()]) for target_vector in self._target_vectors_prim.values()]
            partial_dual_augmented = [np.hstack([partial_dual, target_vector.as_column_vector()]) for target_vector in self._target_vectors_dual.values()]

            # Multithreaded Computation of Ranks
            with PoolToUse() as pool:
                results = pool.imap(np.linalg.matrix_rank, [partial_prim] + partial_prim_augmented + [partial_dual] + partial_dual_augmented)

                rank_list = [rank for rank in results]

            border_index = len(self._target_vectors_prim) + 1

            ranks_prim = rank_list[:border_index]
            ranks_dual = rank_list[border_index:]

            verify_prim_each = {name: ranks_prim[0] == ranks_prim[index + 1] for index, name in enumerate(self._target_vectors_prim.keys())}
            verify_dual_each = {name: ranks_dual[0] == ranks_dual[index + 1] for index, name in enumerate(self._target_vectors_dual.keys())}

            return verify_prim_each, verify_dual_each

        def verification_each_multithreaded() -> tuple[dict[str, bool], dict[str, bool]]:
            return verification_each_parallelized(PoolToUse = ThreadPool)

        def verification_each_multiprocessed() -> tuple[dict[str, bool], dict[str, bool]]:
            return verification_each_parallelized(PoolToUse = ProcessPool)

        def verification_each_adaptive() -> tuple[bool, bool]:
            shape_prim = self.partial(mode = Mode.PRIM, order = Order.FACE, relative_boundary = True, smoothen_dual_lattice = True).shape
            shape_dual = self.partial(mode = Mode.DUAL, order = Order.LINK, relative_boundary = True, smoothen_dual_lattice = True).shape

            # 2500 Axes Was Chosen As a Threshold Based On Anecdotal Evidence from Really Rough Test Runs. More Finetuned Testing Might Yield a Better Heuristic.
            if np.average([*shape_prim, *shape_dual]) > 2500:
                return verification_each_multiprocessed()
            else:
                return verification_each_sequential()

        def verification_bulk_sequential() -> tuple[bool, bool]:
            partial_prim = self.partial(mode = Mode.PRIM, order = Order.FACE, relative_boundary = True, smoothen_dual_lattice = True)
            partial_dual = self.partial(mode = Mode.DUAL, order = Order.LINK, relative_boundary = True, smoothen_dual_lattice = True)

            partial_prim_augmented = np.hstack([partial_prim] + [target_vector.as_column_vector() for target_vector in self._target_vectors_prim.values()])
            partial_dual_augmented = np.hstack([partial_dual] + [target_vector.as_column_vector() for target_vector in self._target_vectors_dual.values()])

            rank_prim = np.linalg.matrix_rank(partial_prim) # This Is ´Intercepted´ By Galois to Return the Rank over GF2
            rank_dual = np.linalg.matrix_rank(partial_dual) # This Is ´Intercepted´ By Galois to Return the Rank over GF2

            rank_prim_augmented = np.linalg.matrix_rank(partial_prim_augmented) # This Is ´Intercepted´ By Galois to Return the Rank over GF2
            rank_dual_augmented = np.linalg.matrix_rank(partial_dual_augmented) # This Is ´Intercepted´ By Galois to Return the Rank over GF2

            verify_prim_bulk = (rank_prim == rank_prim_augmented)
            verify_dual_bulk = (rank_dual == rank_dual_augmented)

            return verify_prim_bulk, verify_dual_bulk

        def verification_bulk_parallelized(PoolToUse:Type[Pool]) -> tuple[bool, bool]:
            partial_prim = self.partial(mode = Mode.PRIM, order = Order.FACE, relative_boundary = True, smoothen_dual_lattice = True)
            partial_dual = self.partial(mode = Mode.DUAL, order = Order.LINK, relative_boundary = True, smoothen_dual_lattice = True)

            partial_prim_augmented = np.hstack([partial_prim] + [target_vector.as_column_vector() for target_vector in self._target_vectors_prim.values()])
            partial_dual_augmented = np.hstack([partial_dual] + [target_vector.as_column_vector() for target_vector in self._target_vectors_dual.values()])

            # Parallelized Computation of Ranks
            with PoolToUse() as pool:
                results = pool.imap(np.linalg.matrix_rank, [partial_prim, partial_dual, partial_prim_augmented, partial_dual_augmented])

                rank_list = [rank for rank in results]

            verify_prim_bulk = (rank_list[0] == rank_list[2])
            verify_dual_bulk = (rank_list[1] == rank_list[3])

            return verify_prim_bulk, verify_dual_bulk

        def verification_bulk_multithreaded() -> tuple[bool, bool]:
            return verification_bulk_parallelized(ThreadPool)

        def verification_bulk_multiprocessed() -> tuple[bool, bool]:
            return verification_bulk_parallelized(ProcessPool)

        def verification_bulk_adaptive() -> tuple[bool, bool]:
            shape_prim = self.partial(mode = Mode.PRIM, order = Order.FACE, relative_boundary = True, smoothen_dual_lattice = True).shape
            shape_dual = self.partial(mode = Mode.DUAL, order = Order.LINK, relative_boundary = True, smoothen_dual_lattice = True).shape

            # 5000 Axes Was Chosen As a Threshold Based On Anecdotal Evidence from Really Rough Test Runs. More Finetuned Testing Might Yield a Better Heuristic.
            if np.average([*shape_prim, *shape_dual]) > 5000:
                return verification_bulk_multiprocessed()
            else:
                return verification_bulk_sequential()

        # Lists of All Available Methods for Verification
        methods_each = [verification_each_sequential, verification_each_multithreaded, verification_each_multiprocessed, verification_each_adaptive]
        methods_bulk = [verification_bulk_sequential, verification_bulk_multithreaded, verification_bulk_multiprocessed, verification_bulk_adaptive]
        methods = methods_each + methods_bulk

        # Check Input
        if safe_mode:
            assert method is None or (type(method) is int and 0 <= method and method < len(methods))

        each_method_selected = (type(method) is int and method < len(methods_each))
        bulk_method_selected = (type(method) is int and method >= len(methods_each))

        # INDIVIDUAL VERIFICATION VALUES ARE DESIRED
        if return_each:
            # 1st: Bulk-Verification Method Was Wished For
            not_done = any([v is None for v in (self._verify_prim_each | self._verify_dual_each).values()])
            if not_done and bulk_method_selected:
                # Check If Bulk Verification Was Not Yet Attempted Before Executing
                if self._verify_prim_bulk is None and self._verify_dual_bulk is None:
                    self._verify_prim_bulk, self._verify_dual_bulk = methods[method]()
                    self._verify_prim_each = {name:(True if self._verify_prim_bulk else None) for name in self._target_vectors_prim.keys()}
                    self._verify_dual_each = {name:(True if self._verify_dual_bulk else None) for name in self._target_vectors_dual.keys()}

            # 2nd: Each-Verification Method Was Wished For
            not_done = any([v is None for v in (self._verify_prim_each | self._verify_dual_each).values()])
            if not_done and each_method_selected:
                self._verify_prim_each, self._verify_dual_each = methods[method]()
                self._verify_prim_bulk = all(self._verify_prim_each.values())
                self._verify_dual_bulk = all(self._verify_dual_each.values())

            # 3rd: Default Option: Adaptive Each-Verification
            not_done = any([v is None for v in (self._verify_prim_each | self._verify_dual_each).values()])
            if not_done:
                self._verify_prim_each, self._verify_dual_each = verification_each_adaptive()
                self._verify_prim_bulk = all(self._verify_prim_each.values())
                self._verify_dual_bulk = all(self._verify_dual_each.values())

            assert not(not_done := any([v is None for v in (self._verify_prim_each | self._verify_dual_each).values()]))
            return self._verify_prim_each | self._verify_dual_each

        # TOTAL VERIFICATION VALUE IS DESIRED
        else:
            # 1st: Bulk-Verification Method Was Wished For
            not_done = any([self._verify_prim_bulk is None, self._verify_dual_bulk is None])
            if not_done and bulk_method_selected:
                self._verify_prim_bulk, self._verify_dual_bulk = methods[method]()
                self._verify_prim_each = {name:(True if self._verify_prim_bulk else None) for name in self._target_vectors_prim.keys()}
                self._verify_dual_each = {name:(True if self._verify_dual_bulk else None) for name in self._target_vectors_dual.keys()}

            # 2nd: Each-Verification Method Was Wished For
            not_done = any([self._verify_prim_bulk is None, self._verify_dual_bulk is None])
            if not_done and each_method_selected:
                if safe_mode:
                    warn("While Accessing a Total Verification Value, an Individual Verification Method was Wished for ... are you sure?")
                self._verify_prim_each, self._verify_dual_each = methods[method]()
                self._verify_prim_bulk = all(self._verify_prim_each.values())
                self._verify_dual_bulk = all(self._verify_dual_each.values())

            # 3rd: Default Option: Adaptive Bulk-Verification
            not_done = any([self._verify_prim_bulk is None, self._verify_dual_bulk is None])
            if not_done:
                self._verify_prim_bulk, self._verify_dual_bulk = verification_bulk_adaptive()
                self._verify_prim_each = {name:(True if self._verify_prim_bulk else None) for name in self._target_vectors_prim.keys()}
                self._verify_dual_each = {name:(True if self._verify_dual_bulk else None) for name in self._target_vectors_dual.keys()}

            assert not(not_done := any([self._verify_prim_bulk is None, self._verify_dual_bulk is None]))
            return self._verify_prim_bulk and self._verify_dual_bulk

    def partial(self, mode:Mode, order:Order, relative_boundary:bool = False, smoothen_dual_lattice:bool = True, safe_mode:bool = safe_mode_global) -> FieldArray:
        # Check Input
        if safe_mode:
            assert is_mode(mode)
            assert is_order(order)
            assert type(relative_boundary) is bool

        # Check If the Relevant Boundary Operator Exists
        if mode == Mode.PRIM and order == Order.SITE:
            raise RuntimeError("0-Chains Have No Primal Boundary.")
        if mode == Mode.DUAL and order == Order.CUBE:
            raise RuntimeError("3-Chains Have No Dual Boundary.")

        # Check if Relevant Boundary Operator Was Calculated
        if self._partial[mode, order, relative_boundary] is None:

            # Initialize the Boundary Operator
            order_in = order
            order_to = Order(order.value - 1) if mode == Mode.PRIM else Order(order.value + 1)

            dimension_in = dimensions(shape = self.shape, order = order_in)
            dimension_to = dimensions(shape = self.shape, order = order_to)

            partial = GF2.Zeros(shape = (dimension_to, dimension_in))

            # Calculate the Boundary Operator
            for index_in in range(dimension_in):
                tuple_in = self.bijection(order = order_in).inverse[index_in]
                boundary = Chain.fromCoordinates(tuple_in, lattice = self).boundary(mode = mode, relative_boundary = relative_boundary, smoothen_dual_lattice = smoothen_dual_lattice)
                partial[:, index_in] = boundary.as_vector()

            self._partial[mode, order, relative_boundary] = partial

        # Return the Existing Relevant Boundary Operator
        return self._partial[mode, order, relative_boundary]

    def all_coordinates(self, order: Order, safe_mode:bool = safe_mode_global) -> list[Coordinates]:
        # Check Input
        if safe_mode:
            assert is_order(order)

        # Create List of All Structures/Coordinates of the Given Order
        """
        Note: Why sometimes use np.arange followed by mapping to floats? We want tuples rather than np.ndarrays because they are hashable. However,
        we cannot generate them all as products of range(...) because range doesn't support half-integers. However, when using numpy.arange instead,
        numpy will default to np.float64, which -- when displayed within a standard tuple -- will become unreadable, e.g. in debugging.
        """
        coordinates_list = list()
        if order == Order.SITE:
            coordinates_list.extend( product(range(self.shape[0]), range(self.shape[1]), range(self.shape[2])) )
        elif order == Order.LINK:
            coordinates_list.extend( product( map(float, np.arange(0.5, self.shape[0] - 1)), range(self.shape[1]), range(self.shape[2]) ))
            coordinates_list.extend( product( range(self.shape[0]), map(float, np.arange(0.5, self.shape[1] - 1)), range(self.shape[2]) ))
            coordinates_list.extend( product( range(self.shape[0]), range(self.shape[1]), map(float, np.arange(0.5, self.shape[2] - 1)) ))
        elif order == Order.FACE:
            coordinates_list.extend( product( range(self.shape[0]), map(float, np.arange(0.5, self.shape[1] - 1)), map(float, np.arange(0.5, self.shape[2] - 1)) ))
            coordinates_list.extend( product( map(float, np.arange(0.5, self.shape[0] - 1)), range(self.shape[1]), map(float, np.arange(0.5, self.shape[2] - 1)) ))
            coordinates_list.extend( product( map(float, np.arange(0.5, self.shape[0] - 1)), map(float, np.arange(0.5, self.shape[1] - 1)), range(self.shape[2]) ))
        elif order == Order.CUBE:
            coordinates_list.extend( product(map(float, np.arange(0.5, self.shape[0] - 1)), map(float, np.arange(0.5, self.shape[1] - 1)), map(float, np.arange(0.5, self.shape[2] - 1))))

        # Sanity-Check Output Length
        if safe_mode:
            length, dimension = len(coordinates_list), dimensions(shape = self.shape, order = order)
            assert length == dimension, f"The number of structures {length} doesn't match the dimension of the space {dimension}"

        return coordinates_list

    def bijection(self, order:Order, safe_mode:bool = safe_mode_global) -> bidict:
        # Check Input
        if safe_mode:
            assert is_order(order)

        # Check if Relevant Bijection Was Calculated
        if self._bijections[order] is None:
            # Enumerate All Coordinates and Save Them
            new_bijection = bidict()
            for index, coordinates in enumerate(self.all_coordinates(order)):
                new_bijection[coordinates] = index

            self._bijections[order] = new_bijection

        return self._bijections[order]


############
### MAIN ###
############

def main():
    return

if __name__ == "__main__":
    main()
