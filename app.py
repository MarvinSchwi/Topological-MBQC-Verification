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
from typing import Any, Callable
from enum import Enum

# Logic
from main import *

# Tkinter
from tkinter import Tk
from tkinter.filedialog import askopenfilename

# Panda3D
from panda3d.core import PandaNode, WindowProperties, TransparencyAttrib, AntialiasAttrib
from panda3d.core import Geom, GeomNode, GeomVertexFormat, GeomVertexData, GeomVertexWriter, GeomLines, GeomTriangles
from direct.showbase.ShowBase import ShowBase
from direct.filter.CommonFilters import CommonFilters

# Parsing
from argparse import ArgumentParser


##############
### TYPING ###
##############

# Global Toggle for Type Checking
safe_mode_global:bool = False

# Define DrawType and a Function to Check
class DrawType(Enum):
    LATTICE = 0
    L = 0
    DEFECT = 1
    D = 1
    TARGET = 2
    T = 2
    SURFACE = 3
    S = 3

def is_drawtype(drawtype:Any) -> bool:
    return isinstance(drawtype, DrawType)

# Global Dictionaries for Default Draw Parameters
draw_args_global = {
    # PRIMAL OBJECTS
    (DrawType.L, Order.LINK): {"color_4f": (0.8, 0.8, 0.8, 0.3)},
    (DrawType.D, Order.LINK): {"color_4f": (0.0, 0.0, 0.0, 1.0)},
    (DrawType.T, Order.LINK): {"color_4f": (0.85, 0.13, 0.13, 0.80)},
    (DrawType.S, Order.FACE): {"color_4f": (0.3, 0.3, 0.3, 0.1)},
    # DUAL OBJECTS
    (DrawType.L, Order.FACE): {"color_4f": (0.2, 0.2, 0.2, 0.3)},
    (DrawType.D, Order.FACE): {"color_4f": (0.3, 0.3, 1.0, 1.0)},
    (DrawType.T, Order.FACE): {"color_4f": (1.00, 0.65, 0.00, 0.80)},
    (DrawType.S, Order.LINK): {"color_4f": (0.0, 0.0, 0.6, 0.1)},
}

thickness_args_global = {
    DrawType.L: 1,
    DrawType.D: 3,
    DrawType.T: 5,
}

background_color_global = (0.6, 0.6, 0.6)

# Define a Main Function That Checks Whether a PandaNode Has a Specified Set of Python Tags of Correct Type
# The Reason This is Done Instead of Simply Subclassing PandaNode is Due to Subclassing Not Working Properly
# in Panda3d Due to Its C++ Backend; see https://docs.panda3d.org/1.10/python/programming/object-management/subclassing.
def is_nodeType(node:Any, tag_types:dict[str, Callable]) -> bool:
    # Check That The node is Of Type PandaNode
    if not isinstance(node, PandaNode):
        return False

    for tag, check_function in tag_types.items():
        # Check That The node has the Specified Tag ...
        if not(node.hasPythonTag(tag)):
            return False
        # ... and That The Item is Either None or of the Specified Type.
        if not(node.getPythonTag(tag) is None or check_function(node.getPythonTag(tag))):
            return False

    return True

# Define ToggleNode (as a PandaNode with Specified Python Tags)
type ToggleNode = PandaNode

def is_toggleNode(toggleNode:Any) -> bool:
    toggle_tag_types = {"toggle_methods": lambda x: isinstance(x, dict) and all([isinstance(y, int) for y in x.keys()]) and all([isinstance(y, Callable) for y in x.values()]),
                        "toggle_options": lambda x: isinstance(x, dict) and all([isinstance(y, int) for y in x.keys()]) and all([isinstance(y, list) for y in x.values()]),
                        "toggle_current": lambda x: isinstance(x, dict) and all([isinstance(y, int) for y in x.keys()]) and all([isinstance(y, str) for y in x.values()])}

    return is_nodeType(node = toggleNode, tag_types = toggle_tag_types)

# Define ChainNode (as a PandaNode with Specified Python Tags)
type ChainNode = PandaNode

def is_chainNode(chainNode:Any) -> bool:
    chain_tag_types = {"chain_ref": lambda x: isinstance(x, Chain),
                       "state_id_current": lambda x: isinstance(x, tuple) and all([isinstance(y, Real) for y in x])}

    return is_nodeType(node = chainNode, tag_types = chain_tag_types)

# Define StructureNode (as a PandaNode with Specified Python Tags)
type StructureNode = PandaNode

def is_structureNode(structureNode:Any) -> bool:
    structure_tag_types = {"coordinates": lambda x: isinstance(x, tuple) and len(x) == 3 } # is_coordinates Function Would Be Better, But Requires Shape or Lattice

    return is_nodeType(node = structureNode, tag_types = structure_tag_types)


########################
### HELPER-FUNCTIONS ###
########################

def generateLinkNode(points:set[tuple[Real, Real, Real]], **kwargs:Any) -> GeomNode:
    # Check Input
    assert len(points) == 2

    # Obtain and Organize the Raw Data
    p1, p2 = points
    color_4f = kwargs["color_4f"]

    # Define the Vertex Format Including Color and Obtain the Necessary Writers
    format = GeomVertexFormat.getV3c4()
    vdata = GeomVertexData("vertexData", format, Geom.UHStatic)
    vertex = GeomVertexWriter(vdata, "vertex")
    color = GeomVertexWriter(vdata, "color")

    # Write the Positions and Colors
    vertex.addData3f(p1)
    vertex.addData3f(p2)

    color.addData4f(*color_4f)
    color.addData4f(*color_4f)

    # Create Geometry and Line Primitive
    geom_lines = GeomLines(Geom.UHStatic)
    geom_lines.addVertices(0, 1)

    # Create Geometry and Node
    geom = Geom(vdata)
    geom.addPrimitive(geom_lines)
    geom_node = GeomNode(f"link {points} {color_4f}")
    geom_node.addGeom(geom)

    return geom_node

def generateFaceNode(points:set[tuple[Real, Real, Real]], **kwargs:Any) -> GeomNode:
    # Check Input
    assert len(points) == 4

    # Obtain and Organize the Raw Data
    points = list(points)

    ind_e = [all([points[p_index][c_index] == points[p_index + 1][c_index] for p_index in range(3)]) for c_index in range(3)].index(True)
    ind_1, ind_2 = {0, 1, 2} - {ind_e}
    min_1, min_2 = min([point[ind_1] for point in points]), min([point[ind_2] for point in points])

    p1 = [point for point in points if point[ind_1] == min_1 and point[ind_2] == min_2]
    p2 = [point for point in points if point[ind_1] != min_1 and point[ind_2] == min_2]
    p3 = [point for point in points if point[ind_1] != min_1 and point[ind_2] != min_2]
    p4 = [point for point in points if point[ind_1] == min_1 and point[ind_2] != min_2]

    points = p1 + p2 + p3 + p4

    color_4f = kwargs["color_4f"]

    # Define the Vertex Format Including Color and Obtain the Necessary Writers
    format = GeomVertexFormat.getV3c4()
    vdata = GeomVertexData("vertexData", format, Geom.UHStatic)
    vertex = GeomVertexWriter(vdata, "vertex")
    color = GeomVertexWriter(vdata, "color")

    # Write the Positions and Colors
    for point in points:
        vertex.addData3f(*point)
        color.addData4f(*color_4f)

    # Create Geometry and Line Primitive
    geom_triangles = GeomTriangles(Geom.UHStatic)
    geom_triangles.addVertices(0, 1, 2)
    geom_triangles.addVertices(2, 3, 0)

    # Create Geometry and Node
    geom = Geom(vdata)
    geom.addPrimitive(geom_triangles)
    geom_node = GeomNode(f"face {points} {color_4f}")
    geom_node.addGeom(geom)

    return geom_node

def generateCoordinatesNode(coordinates:Coordinates, lattice:Lattice, drawtype:DrawType, draw_args:dict[tuple[DrawType, Order], dict], safe_mode:bool = safe_mode_global) -> StructureNode:
    # Determine Order of Coordinates
    order = determine_order(coordinates)

    # Check Input
    if safe_mode:
        assert isinstance(lattice, Lattice)
        assert is_coordinates(coordinates, lattice = lattice)
        assert is_drawtype(drawtype)
        assert order in {Order.LINK, Order.FACE}

    selected_drawargs = draw_args[drawtype, order]

    if drawtype == DrawType.LATTICE:
        """
        One Doesn't Always Want to Draw All Possible Lattice Links.
        - In the Primal Lattice, All Links With Both End-Points Within the Shape Should Be Drawn. (This Happens to Coincide With All Links Because All Primal Links Have Both Endpoints in the Lattice.)
        - In the Dual Lattice, All Links Within the "Smoothened" Lattice Should Be Drawn. This Constitutes All Dual Links (i.e., Primal Faces) Not at the Edge.
        """
        # Determine the Points Defining the Link
        mode = Mode.PRIM if order == Order.LINK else Mode.DUAL
        points = Chain.fromCoordinates(coordinates, lattice).boundary(mode = mode, all_limits = True)

        # Smoothen The Dual Lattice
        if mode == Mode.DUAL and is_at_edge(coordinates = coordinates, lattice = lattice):
            return None

        # Construct the GeomNode for the Link
        coordinatesNode = generateLinkNode(points, **selected_drawargs)

    elif drawtype == DrawType.DEFECT:
        """
        One Always Wants to Draw All Defects, i.e., We Use all_limits = True.
        """
        # Determine the Points Defining the Link
        mode = Mode.PRIM if order == Order.LINK else Mode.DUAL
        points = Chain.fromCoordinates(coordinates, lattice).boundary(mode = mode, all_limits = True)

        # Construct the GeomNode for the Link
        coordinatesNode = generateLinkNode(points, **selected_drawargs)

    elif drawtype == DrawType.TARGET:
        """
        One Always Wants to Draw All Targets, i.e., We Use all_limits = True.
        TODO (Low Priority): That Being Said, Due to How Target Vectors Are Normally Are Constructed, This Should Be Equal to the More Restrictive Boundary Point Pairs. (Right?) One Could/Should Assert That.
        """
        # Determine the Points Defining the Link
        mode = Mode.PRIM if order == Order.LINK else Mode.DUAL
        points = Chain.fromCoordinates(coordinates, lattice).boundary(mode = mode, all_limits = True)

        # Construct the GeomNode for the Link
        coordinatesNode = generateLinkNode(points, **selected_drawargs)

    elif drawtype in {DrawType.SURFACE}:
        """
        One Always Wants to Draw All Surfaces, i.e., We Use all_limits = True.
        TODO (Medium to High Priority): Rather Said, We Want to Use all_limits = True. However, Here We Need the Union of the Boundary of the Individual Coordinates of the Boundary of the Face.
             Since We Can't Determine the Second Boundary For Arbitrary Points (Yet?), We Simply Use all_limits = False. In General, These Two Should Probably Be Equivalent. (Right?)
             But to Confirm This I'll Need to Think About It Again.
        TODO (Low Priority): Write Union Function for Chains.
        """
        # Determine the Points Defining the Face
        mode = Mode.PRIM if order == Order.FACE else Mode.DUAL
        boundary = Chain.fromCoordinates(coordinates, lattice).boundary(mode = mode, all_limits = True)
        points = {coordinates for structure in boundary if is_coordinates(coordinates = structure, lattice = lattice) for coordinates in Chain.fromCoordinates(structure, lattice).boundary(mode = mode, all_limits = True)}

        # Construct the GeomNode for the Face
        coordinatesNode = generateFaceNode(points, **selected_drawargs)

    else:
        raise RuntimeError(f"Invalid DrawType \'{drawtype}\'")

    # Augment Relevant Python Tags
    coordinatesNode.setPythonTag("coordinates", coordinates)

    # Check Output
    if safe_mode:
        assert is_structureNode(coordinatesNode)

    return coordinatesNode

def generateChainNode(chain:Chain, drawtype:DrawType, draw_args:dict[tuple[DrawType, Order], dict], safe_mode:bool = safe_mode_global) -> ChainNode:
    # Return a Trivial ChainNode if chain is None
    if chain is None:
        if safe_mode:
            warn("Generating dead ChainNode for a None chain argument.")
        NoneNode = PandaNode("Dead ChainNode.")
        NoneNode.setPythonTag("chain_ref", None)
        NoneNode.setPythonTag("state_id_current", None)
        return NoneNode

    # Check Input
    if safe_mode:
        assert isinstance(chain, Chain)
        assert is_drawtype(drawtype)

    chainNode = PandaNode(f"chain; original state_id: {chain.state_id}")

    # Attach a coordinatesNode For All Coordinates (Whose Defining Points Exist in the Lattice)
    for coordinates in chain.as_tuples():
        if (coordinatesNode := generateCoordinatesNode(coordinates, chain._lattice, drawtype, draw_args)) is None:
            continue

        chainNode.addChild( coordinatesNode )

    # Augment Relevant Python Tags
    chainNode.setPythonTag("chain_ref", chain)
    chainNode.setPythonTag("state_id_current", chain.state_id)

    # TODO: Add a Task for Updating

    # Check Output
    if safe_mode:
        assert is_chainNode(chainNode)

    return chainNode

def generateLatticesNode(lattice:Lattice, draw_args:dict[tuple[DrawType, Order], dict], safe_mode:bool = safe_mode_global) -> ToggleNode:
    # Check Input
    if safe_mode:
        assert isinstance(lattice, Lattice)
        assert isinstance(draw_args, dict)

    # Create the ToggleNode for the Lattice
    main_lattice_node = PandaNode("main lattice")

    # Generate and Attach the ChainNodes for the Primal and Dual Lattice
    prim_lattice = Chain.fromVector(vector = GF2.Ones(shape = dimensions(shape = lattice.shape, order = Order.LINK)), lattice = lattice)
    dual_lattice = Chain.fromVector(vector = GF2.Ones(shape = dimensions(shape = lattice.shape, order = Order.FACE)), lattice = lattice)

    prim_lattice_node = generateChainNode(chain = prim_lattice, drawtype = DrawType.LATTICE, draw_args = draw_args)
    dual_lattice_node = generateChainNode(chain = dual_lattice, drawtype = DrawType.LATTICE, draw_args = draw_args)

    prim_lattice_node.overall_hidden = False
    dual_lattice_node.overall_hidden = True

    main_lattice_node.addChild(prim_lattice_node)
    main_lattice_node.addChild(dual_lattice_node)

    # Augment Python Tags for the Back-End
    main_lattice_node.setPythonTag("child_nodes", {"prim": prim_lattice_node, "dual": dual_lattice_node})

    # Augment Python Tags for Toggling
    def toggle_lattice(toggle_option_:str = None) -> None:
        # Determine the Current Option and the Next Option
        toggle_options = main_lattice_node.getPythonTag("toggle_options")[0]
        toggle_current = main_lattice_node.getPythonTag("toggle_current")[0]

        assert toggle_option_ is None or toggle_option_ in toggle_options

        if toggle_option_ is None:
            toggle_option_ = toggle_options[(toggle_options.index(toggle_current) + 1) % len(toggle_options)]

        # Hide Current Option and Show Next Option
        child_nodes = main_lattice_node.getPythonTag("child_nodes")

        if not (toggle_current == "none"):
            child_nodes[toggle_current].overall_hidden = True

        if not (toggle_option_ == "none"):
            child_nodes[toggle_option_].overall_hidden = False

        # Update the State
        main_lattice_node.setPythonTag("toggle_current", {0: toggle_option_})

    main_lattice_node.setPythonTag("toggle_methods", {0: toggle_lattice})
    main_lattice_node.setPythonTag("toggle_options", {0: ["none", "prim", "dual"]})
    main_lattice_node.setPythonTag("toggle_current", {0: "prim"})

    # Check Output
    if safe_mode:
        assert is_toggleNode(main_lattice_node)

    return main_lattice_node

def generateDefectsNode(lattice:Lattice, draw_args:dict[tuple[DrawType, Order], dict], safe_mode:bool = safe_mode_global) -> ToggleNode:
    # Check Input
    if safe_mode:
        assert isinstance(lattice, Lattice)
        assert isinstance(draw_args, dict)

    # Create the ToggleNode for the Lattice
    main_defects_node = PandaNode("main defects")

    # Generate and Attach the ChainNodes for the Primal and Dual Lattice
    prim_defects_node = generateChainNode(chain = lattice.defects(order = Order.LINK), drawtype = DrawType.DEFECT, draw_args = draw_args)
    dual_defects_node = generateChainNode(chain = lattice.defects(order = Order.FACE), drawtype = DrawType.DEFECT, draw_args = draw_args)

    prim_defects_node.overall_hidden = False
    dual_defects_node.overall_hidden = False

    main_defects_node.addChild(prim_defects_node)
    main_defects_node.addChild(dual_defects_node)

    # Augment Python Tags for the Back-End
    main_defects_node.setPythonTag("child_nodes", {"prim": prim_defects_node, "dual": dual_defects_node})

    # Augment Python Tags for Toggling
    def toggle_defects(toggle_option_:str = None) -> None:
        # Determine the Current Option and the Next Option
        toggle_options = main_defects_node.getPythonTag("toggle_options")[0]
        toggle_current = main_defects_node.getPythonTag("toggle_current")[0]

        assert toggle_option_ is None or toggle_option_ in toggle_options

        if toggle_option_ is None:
            toggle_option_ = toggle_options[(toggle_options.index(toggle_current) + 1) % len(toggle_options)]

        # Hide Current Option and Show Next Option
        child_nodes = main_defects_node.getPythonTag("child_nodes")

        if not (toggle_current in {"none", "all"}):
            child_nodes[toggle_current].overall_hidden = True
        elif toggle_current == "all":
            for child_node in child_nodes.values():
                child_node.overall_hidden = True

        if not (toggle_option_ in {"none", "all"}):
            child_nodes[toggle_option_].overall_hidden = False
        elif toggle_option_ == "all":
            for child_node in child_nodes.values():
                child_node.overall_hidden = False

        # Update the State
        main_defects_node.setPythonTag("toggle_current", {0: toggle_option_})

    main_defects_node.setPythonTag("toggle_methods", {0: toggle_defects})
    main_defects_node.setPythonTag("toggle_options", {0: ["none", "prim", "dual", "all"]})
    main_defects_node.setPythonTag("toggle_current", {0: "all"})

    # Check Output
    if safe_mode:
        assert is_toggleNode(main_defects_node)

    return main_defects_node

def generateObjectiveNode(lattice:Lattice, draw_args:dict[tuple[DrawType, Order], dict], safe_mode:bool = safe_mode_global) -> ToggleNode:
    # Check Input
    if safe_mode:
        assert isinstance(lattice, Lattice)
        assert isinstance(draw_args, dict)

    # Create the ToggleNode for the Lattice
    main_objective_node = PandaNode("main objective")

    # Generate and Attach Main Nodes for the Target Vectors and Surfaces

    main_target_v_node = PandaNode("main target vectors")
    main_surfaces_node = PandaNode("main surfaces")

    main_target_v_node.overall_hidden = False # This Always Stays False
    main_surfaces_node.overall_hidden = True

    main_objective_node.addChild(main_target_v_node)
    main_objective_node.addChild(main_surfaces_node)

    # Generate and Attach the ChainNodes for the Target Vectors and Surfaces
    target_v_nodes = {name: generateChainNode(chain = target_v, drawtype = DrawType.T, draw_args = draw_args)
                      for name, target_v in (lattice._target_vectors_prim | lattice._target_vectors_dual).items()}
    surfaces_nodes = {name: generateChainNode(chain = surface_, drawtype = DrawType.S, draw_args = draw_args)
                      for name, surface_ in (lattice._surfaces_prim | lattice._surfaces_dual).items()}

    for target_v_node in target_v_nodes.values():
        main_target_v_node.addChild(target_v_node)
        target_v_node.overall_hidden = True

    for surface__node in surfaces_nodes.values():
        main_surfaces_node.addChild(surface__node)
        surface__node.overall_hidden = True

    # Augment Python Tags for the Back-End
    main_objective_node.setPythonTag("t_child_nodes", target_v_nodes)
    main_objective_node.setPythonTag("s_child_nodes", surfaces_nodes)

    # Augment Python Tags for Toggling
    def toggle_targets(toggle_option_:str = None) -> None:
        # Determine the Current Option and the Next Option
        toggle_options = main_objective_node.getPythonTag("toggle_options")[0]
        toggle_current = main_objective_node.getPythonTag("toggle_current")[0]

        assert toggle_option_ is None or toggle_option_ in toggle_options

        if toggle_option_ is None:
            toggle_option_ = toggle_options[(toggle_options.index(toggle_current) + 1) % len(toggle_options)]

        # Hide Current Option and Show Next Option
        t_child_nodes = main_objective_node.getPythonTag("t_child_nodes")
        s_child_nodes = main_objective_node.getPythonTag("s_child_nodes")

        if not (toggle_current == "none"):
            t_child_nodes[toggle_current].overall_hidden = True
            s_child_nodes[toggle_current].overall_hidden = True

        if not (toggle_option_ == "none"):
            t_child_nodes[toggle_option_].overall_hidden = False
            s_child_nodes[toggle_option_].overall_hidden = False

        # Update the State
        main_objective_node.setPythonTag("toggle_current", {0: toggle_option_, 1:main_objective_node.getPythonTag("toggle_current")[1]})

    def toggle_surfaces(toggle_option_:str = None) -> None:
        # Determine the Current Option and the Next Option
        toggle_options = main_objective_node.getPythonTag("toggle_options")[1]
        toggle_current = main_objective_node.getPythonTag("toggle_current")[1]

        assert toggle_option_ is None or toggle_option_ in toggle_options

        if toggle_option_ is None:
            toggle_option_ = toggle_options[(toggle_options.index(toggle_current) + 1) % len(toggle_options)]

        # Hide Current Option and Show Next Option
        if toggle_option_ == "on":
            main_surfaces_node.overall_hidden = False
        elif toggle_option_ == "off":
            main_surfaces_node.overall_hidden = True

        # Update the State
        main_objective_node.setPythonTag("toggle_current", {0:main_objective_node.getPythonTag("toggle_current")[0], 1: toggle_option_})

    main_objective_node.setPythonTag("toggle_methods", {0: toggle_targets, 1: toggle_surfaces})
    main_objective_node.setPythonTag("toggle_options", {0: ["none"] + list((lattice._target_vectors_prim | lattice._target_vectors_dual).keys()), 1: ["off", "on"]})
    main_objective_node.setPythonTag("toggle_current", {0: "none", 1: "off"})

    # Check Output
    if safe_mode:
        assert is_toggleNode(main_objective_node)

    return main_objective_node


##################
### VISUAL APP ###
##################

class LatticeApp(ShowBase):

    def __init__(self, lattice:Lattice = None, fullscreen:bool = False, background_color:tuple[Real] = background_color_global) -> None:
        ShowBase.__init__(self)

        # Set Up the Screen
        width, height = self.pipe.getDisplayWidth(), self.pipe.getDisplayHeight()
        properties = WindowProperties(fullscreen = fullscreen, size = (width, height) if fullscreen else (width//2, height//2))
        self.win.requestProperties(properties)

        # Disable HDR
        filters = CommonFilters(self.win, self.cam)
        filters.delHighDynamicRange()

        # Set the Background Color
        self.setBackgroundColor(background_color)

        # Keep Track of Whether a Lattice is Loaded
        self.lattice_initialized = False

        # Initialize Lattice
        self.initialize_lattice(lattice)

        # Draw Both Sides of All Objects
        # See https://discourse.panda3d.org/t/texture-front-and-back-solved/12734.
        self.render.setTwoSided(True)

        # Allow Transparent Objects
        self.render.setTransparency(TransparencyAttrib.MAlpha)

        # Enable Antialiasing
        self.render.setAntialias(AntialiasAttrib.MAuto)


    def initialize_lattice(self, lattice:Lattice, draw_args:dict[tuple[DrawType, Order], dict] = draw_args_global, thickness_args:dict[DrawType, int] = thickness_args_global) -> None:
        assert lattice is None or isinstance(lattice, Lattice)

        if lattice is None:
            raise NotImplementedError()

        else:
            # If the Lattice Was Already Initialized, Remove All Nodes Except the Camera
            if self.lattice_initialized:
                for child in self.render.getChildren():
                    if not (child == self.camera):
                        child.removeNode()

            # Add the Lattice, Defects, and Objective Nodes
            self.main_lattice_node = generateLatticesNode(lattice = lattice, draw_args = draw_args)
            main_lattice_node_path = self.render.attachNewNode(self.main_lattice_node)
            main_lattice_node_path.setRenderModeThickness(thickness_args[DrawType.L])
            self.accept("l", self.main_lattice_node.getPythonTag("toggle_methods")[0])

            self.main_defects_node = generateDefectsNode(lattice = lattice, draw_args = draw_args)
            main_defects_node_path = self.render.attachNewNode(self.main_defects_node)
            main_defects_node_path.setRenderModeThickness(thickness_args[DrawType.D])
            self.accept("d", self.main_defects_node.getPythonTag("toggle_methods")[0])

            self.main_objective_node = generateObjectiveNode(lattice = lattice, draw_args = draw_args)
            main_objective_node_path = self.render.attachNewNode(self.main_objective_node)
            main_objective_node_path.setRenderModeThickness(thickness_args[DrawType.T])
            self.accept("t", self.main_objective_node.getPythonTag("toggle_methods")[0])
            self.accept("s", self.main_objective_node.getPythonTag("toggle_methods")[1])

            # Add the Function to Open/Load a Lattice
            self.accept("o", self.load_lattice)

            self.lattice_initialized = True

    def load_lattice(self) -> None:
        path_file = Path( askopenfilename(filetypes=[("Lattice JSON files", ".json")]) )

        try:
            lattice = Lattice.load(path = path_file)
            self.initialize_lattice(lattice)
        except Exception as e:
            warn(f"Couldn't load lattice with path {path_file}; got error {e}.")


##############
### PARSER ###
##############

def parse() -> dict[str, Any]:
    parser = ArgumentParser(prog = "MBG app.py", description = "A visualizer for topological circuits.")

    parser.add_argument("save", type = str)
    parser.add_argument("--fullscreen", action = "store_true", default = False, help = "Whether to start the visualizer in fullscreen.")

    return vars(parser.parse_args())


############
### MAIN ###
############

def main():
    # Parse the Arguments and Load the Lattice
    args = parse()
    lattice = Lattice.load(path = Path(args["save"]))

    # Start the App and Run It
    app = LatticeApp(lattice = lattice, fullscreen = args["fullscreen"])
    app.run()

if __name__ == "__main__":
    main()
