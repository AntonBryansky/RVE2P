'''
RVE model:
    N7  _________________  N8
      /|                /|
 N5  /_|_______________/ |
    |  |           N6 |  |
    |  |              |  |
    |  |              |  |
    |  |              |  |
    |  _______________|__|
    | / N3            | / N4
    |/________________|/
    N1                N2
    

  RP 3 (RP_Z)
         RP 2 (RP_Y)
    z   y
    |  / 
    | /
    |/_____x    RP 1 (RP_X)    RP 4 (RP_X_1)    RP 5 (RP_X_2)    RP 6 (RP_X_3)

Preliminary
1) create RUC as a single part
2) set material properties and define sections
3) mesh part with the symmetrical mesh on opposite surfaces

Code process:
1) copy the part to assembly as instance
2) find coordinates of vertices and center
3) calculate strains for each direction of deformation
4) create sets at vertices
5) create sets at surfaces
6) create sets at edges
7) create sets at reference points
8) create jobs and calculate shear moduli:
    9.1) set periodic boundary conditions
    9.2) Set boundary conditions
    9.3) Create field outputs for displacements and reaction forces
    9.4) Create and run job
    9.5) Open ODB
    9.6) Calculate RUC sizes and get force at the reference point
    9.7) Calculate cross-sectional area, axial strains and stress
    9.8) Calculate Shear modulus
    9.9) Enable optimal point of view to the RUC in th Viewport
'''


from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import odbAccess
from abaqusConstants import *


# Field-output version of Calculate_E_nu_G_PBC.py.
# U and RF use the last increment; dimensions use the undeformed mesh.
# All history output requests are removed before each job.

# Input parameters
model_name = 'Model-1'
part_name = 'RUC-HC'
instance_name = 'RUC'
# Parameter that defines accuracy of node coordinates
delta = 1e-6
# Parameters of job
num_CPUs = 16
num_Domains = num_CPUs

# Select options for calculations
calculate_E11 = True
calculate_E22 = True
calculate_E33 = True
calculate_G12 = True
calculate_G13 = True
calculate_G23 = True
calculate_p = True


# Functions
def set_part_to_assembly(model_name, part_name, instance_name):
    RUC_assembly = mdb.models[model_name].rootAssembly
    RUC_assembly.DatumCsysByDefault(CARTESIAN)
    RUC_part = mdb.models[model_name].parts[part_name]
    RUC_assembly.Instance(name=instance_name, part=RUC_part, dependent=ON)

def get_RUC_vertices(model_name, instance_name):
    nodes = mdb.models[model_name].rootAssembly.instances[instance_name].nodes
    X_min = 1e+10
    X_max = -1e+10
    Y_min = 1e+10
    Y_max = -1e+10
    Z_min = 1e+10
    Z_max = -1e+10

    for i in range(len(nodes)):
        x = nodes[i].coordinates[0]
        y = nodes[i].coordinates[1]
        z = nodes[i].coordinates[2]
        if X_min > x:
            X_min = x
        elif X_max < x:
            X_max = x
        
        if Y_min > y:
            Y_min = y
        elif Y_max < y:
            Y_max = y
        
        if Z_min > z:
            Z_min = z
            continue
        if Z_max < z:
            Z_max = z
            continue

    N1 = (X_min, Y_min, Z_min)
    N2 = (X_max, Y_min, Z_min)
    N3 = (X_min, Y_max, Z_min)
    N4 = (X_max, Y_max, Z_min)
    N5 = (X_min, Y_min, Z_max)
    N6 = (X_max, Y_min, Z_max)
    N7 = (X_min, Y_max, Z_max)
    N8 = (X_max, Y_max, Z_max)
    Nc = ((X_max+X_min)/2, (Y_max+Y_min)/2, (Z_max+Z_min)/2)
    
    return (N1, N2, N3, N4, N5, N6, N7, N8, Nc)

def calculate_strains(N1, N8, strain):
    X_min = N1[0]
    Y_min = N1[1]
    Z_min = N1[2]
    X_max = N8[0]
    Y_max = N8[1]
    Z_max = N8[2]
    dist_x = X_max - X_min
    dist_y = Y_max - Y_min
    dist_z = Z_max - Z_min
    return (dist_x*strain, dist_y*strain, dist_z*strain)

def set_vertices_to_sets(model_name, instance_name, delta, N1, N2, N3, N4, N5, N6, N7, N8):
    nodes = mdb.models[model_name].rootAssembly.instances[instance_name].nodes

    N1_nodes = nodes.getByBoundingBox(N1[0] - delta, N1[1] - delta, N1[2] - delta, N1[0] + delta, N1[1] + delta, N1[2] + delta)
    N2_nodes = nodes.getByBoundingBox(N2[0] - delta, N2[1] - delta, N2[2] - delta, N2[0] + delta, N2[1] + delta, N2[2] + delta)
    N3_nodes = nodes.getByBoundingBox(N3[0] - delta, N3[1] - delta, N3[2] - delta, N3[0] + delta, N3[1] + delta, N3[2] + delta)
    N4_nodes = nodes.getByBoundingBox(N4[0] - delta, N4[1] - delta, N4[2] - delta, N4[0] + delta, N4[1] + delta, N4[2] + delta)
    N5_nodes = nodes.getByBoundingBox(N5[0] - delta, N5[1] - delta, N5[2] - delta, N5[0] + delta, N5[1] + delta, N5[2] + delta)
    N6_nodes = nodes.getByBoundingBox(N6[0] - delta, N6[1] - delta, N6[2] - delta, N6[0] + delta, N6[1] + delta, N6[2] + delta)
    N7_nodes = nodes.getByBoundingBox(N7[0] - delta, N7[1] - delta, N7[2] - delta, N7[0] + delta, N7[1] + delta, N7[2] + delta)
    N8_nodes = nodes.getByBoundingBox(N8[0] - delta, N8[1] - delta, N8[2] - delta, N8[0] + delta, N8[1] + delta, N8[2] + delta)
    
    set_N1 = mdb.models[model_name].rootAssembly.Set(name="N1", nodes=N1_nodes)
    set_N2 = mdb.models[model_name].rootAssembly.Set(name="N2", nodes=N2_nodes)
    set_N3 = mdb.models[model_name].rootAssembly.Set(name="N3", nodes=N3_nodes)
    set_N4 = mdb.models[model_name].rootAssembly.Set(name="N4", nodes=N4_nodes)
    set_N5 = mdb.models[model_name].rootAssembly.Set(name="N5", nodes=N5_nodes)
    set_N6 = mdb.models[model_name].rootAssembly.Set(name="N6", nodes=N6_nodes)
    set_N7 = mdb.models[model_name].rootAssembly.Set(name="N7", nodes=N7_nodes)
    set_N8 = mdb.models[model_name].rootAssembly.Set(name="N8", nodes=N8_nodes)
    set_vertices = ('N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8')

    return set_N1, set_N2, set_N3, set_N4, set_N5, set_N6, set_N7, set_N8, set_vertices

def set_RP_point_to_set(model_name, instance_name, coordinates, RP_name):
    RP_center_id = mdb.models[model_name].rootAssembly.ReferencePoint(point=coordinates).id
    RPs = mdb.models[model_name].rootAssembly.referencePoints
    RP_center = RPs[RP_center_id]
    set_RP = mdb.models[model_name].rootAssembly.Set(referencePoints=(RP_center, ), name=RP_name)
    return set_RP

def get_node_sets_edges(model_name, instance_name, delta, N1, N2, N3, N4, N5, N6, N7, N8):
    nodes = mdb.models[model_name].rootAssembly.instances[instance_name].nodes
    X_min = N1[0]
    X_max = N8[0]
    Y_min = N1[1]
    Y_max = N8[1]
    Z_min = N1[2]
    Z_max = N8[2]

    # Edge nodes
    nodes_N1_N2 = nodes.getByBoundingBox(X_min + delta, Y_min - delta, Z_min - delta, X_max - delta, Y_min + delta, Z_min + delta)
    nodes_N1_N3 = nodes.getByBoundingBox(X_min - delta, Y_min + delta, Z_min - delta, X_min + delta, Y_max - delta, Z_min + delta)
    nodes_N1_N5 = nodes.getByBoundingBox(X_min - delta, Y_min - delta, Z_min + delta, X_min + delta, Y_min + delta, Z_max - delta)
    nodes_N2_N4 = nodes.getByBoundingBox(X_max - delta, Y_min + delta, Z_min - delta, X_max + delta, Y_max - delta, Z_min + delta)
    nodes_N2_N6 = nodes.getByBoundingBox(X_max - delta, Y_min - delta, Z_min + delta, X_max + delta, Y_min + delta, Z_max - delta)
    nodes_N3_N4 = nodes.getByBoundingBox(X_min + delta, Y_max - delta, Z_min - delta, X_max - delta, Y_max + delta, Z_min + delta)
    nodes_N3_N7 = nodes.getByBoundingBox(X_min - delta, Y_max - delta, Z_min + delta, X_min + delta, Y_max + delta, Z_max - delta)
    nodes_N4_N8 = nodes.getByBoundingBox(X_max - delta, Y_max - delta, Z_min + delta, X_max + delta, Y_max + delta, Z_max - delta)
    nodes_N5_N6 = nodes.getByBoundingBox(X_min + delta, Y_min - delta, Z_max - delta, X_max - delta, Y_min + delta, Z_max + delta)
    nodes_N5_N7 = nodes.getByBoundingBox(X_min - delta, Y_min + delta, Z_max - delta, X_min + delta, Y_max - delta, Z_max + delta)
    nodes_N6_N8 = nodes.getByBoundingBox(X_max - delta, Y_min + delta, Z_max - delta, X_max + delta, Y_max - delta, Z_max + delta)
    nodes_N7_N8 = nodes.getByBoundingBox(X_min + delta, Y_max - delta, Z_max - delta, X_max - delta, Y_max + delta, Z_max + delta)
    
    set_edges_names = ('edge_N1_N2', 'edge_N1_N3', 'edge_N1_N5', 'edge_N2_N4', 'edge_N2_N6', 'edge_N3_N4', 'edge_N3_N7', 'edge_N4_N8', 'edge_N5_N6',
                       'edge_N5_N7', 'edge_N6_N8', 'edge_N7_N8')
    set_edges_nodes = (nodes_N1_N2, nodes_N1_N3, nodes_N1_N5, nodes_N2_N4, nodes_N2_N6, nodes_N3_N4, nodes_N3_N7, nodes_N4_N8, nodes_N5_N6,
                       nodes_N5_N7, nodes_N6_N8, nodes_N7_N8)
    set_edges = []

    for i in range(len(set_edges_nodes)):
        buff = [];
        direction = 0
        if i == 0 or i == 5 or i == 8 or i == 11:
            direction = 0
        elif i == 1 or i == 3 or i == 9 or i == 10:
            direction = 1
        elif i == 2 or i == 4 or i == 6 or i == 7:
            direction = 2

        for j in range(len(set_edges_nodes[i])):
            x = set_edges_nodes[i][j].coordinates[0]
            y = set_edges_nodes[i][j].coordinates[1]
            z = set_edges_nodes[i][j].coordinates[2]
            buff_name = set_edges_names[i] + '_' + str(j)
            buff_node = nodes.getByBoundingBox(x - delta, y - delta, z - delta, x + delta, y + delta, z + delta)
            mdb.models[model_name].rootAssembly.Set(name=buff_name, nodes=buff_node)
            buff.append((set_edges_nodes[i][j].coordinates[direction], buff_name))

        buff.sort()
        buff_sorted = [x for y, x in buff]
        set_edges.append(buff_sorted)
    return set_edges

def get_node_pair_sets_surfaces(model_name, instance_name, delta, N1, N2, N3, N4, N5, N6, N7, N8):
    nodes = mdb.models[model_name].rootAssembly.instances[instance_name].nodes
    X_min = N1[0]
    X_max = N8[0]
    Y_min = N1[1]
    Y_max = N8[1]
    Z_min = N1[2]
    Z_max = N8[2]

    X_pairs = []
    Y_pairs = []
    Z_pairs = []

    # Surfaces nodes
    for i in range(len(nodes)):
        x = nodes[i].coordinates[0]
        y = nodes[i].coordinates[1]
        z = nodes[i].coordinates[2]

        # Front X plane
        if x >= (X_max - delta) and x <= (X_max + delta) and y < (Y_max-delta) and y > (Y_min+delta) and z < (Z_max-delta) and z > (Z_min+delta):
            name_x = 'pair_' + str(i) + '_x_min'
            name_X = 'pair_' + str(i) + '_x_max'
            node_min = nodes.getByBoundingBox(X_min - delta, y - delta, z - delta, X_min + delta, y + delta, z + delta)
            set_X_min = mdb.models[model_name].rootAssembly.Set(name=name_x, nodes=node_min)
            node_max = nodes.getByBoundingBox(x - delta, y - delta, z - delta, x + delta, y + delta, z + delta)
            set_X_max = mdb.models[model_name].rootAssembly.Set(name=name_X, nodes=node_max)
            pair_names = (name_x, name_X)
            X_pairs.append(pair_names)                   
        # Front Y plane
        if y >= (Y_max - delta) and y <= (Y_max + delta) and x > (X_min + delta) and x < (X_max - delta) and z < (Z_max-delta) and z > (Z_min+delta):
            name_y = 'pair_' + str(i) + '_y_min'
            name_Y = 'pair_' + str(i) + '_y_max'
            node_min = nodes.getByBoundingBox(x - delta, Y_min - delta, z - delta, x + delta, Y_min + delta, z + delta)
            set_Y_min = mdb.models[model_name].rootAssembly.Set(name=name_y, nodes=node_min)
            node_max = nodes.getByBoundingBox(x - delta, y - delta, z - delta, x + delta, y + delta, z + delta)
            set_Y_max = mdb.models[model_name].rootAssembly.Set(name=name_Y, nodes=node_max)
            pair_names = (name_y, name_Y)
            Y_pairs.append(pair_names)
        # Front Z plane
        if z >= (Z_max - delta) and z <= (Z_max + delta) and y < (Y_max-delta) and y > (Y_min+delta) and x > (X_min + delta) and x < (X_max - delta):
            name_z = 'pair_' + str(i) + '_z_min'
            name_Z = 'pair_' + str(i) + '_z_max'
            node_min = nodes.getByBoundingBox(x - delta, y - delta, Z_min - delta, x + delta, y + delta, Z_min + delta)
            set_Z_min = mdb.models[model_name].rootAssembly.Set(name=name_z, nodes=node_min)
            node_max = nodes.getByBoundingBox(x - delta, y - delta, z - delta, x + delta, y + delta, z + delta)
            set_Z_max = mdb.models[model_name].rootAssembly.Set(name=name_Z, nodes=node_max)
            #pair_nodes = (set_Z_min, set_Z_max)
            pair_names = (name_z, name_Z)
            Z_pairs.append(pair_names)
        
    return X_pairs, Y_pairs, Z_pairs

def get_surface_pairs(model_name, instance_name, delta, N1, N2, N3, N4, N5, N6, N7, N8):
    nodes = mdb.models[model_name].rootAssembly.instances[instance_name].nodes
    X_min = N1[0]
    X_max = N8[0]
    Y_min = N1[1]
    Y_max = N8[1]
    Z_min = N1[2]
    Z_max = N8[2]

    X_surf_pairs = []
    Y_surf_pairs = []
    Z_surf_pairs = []

    #X
    nodes_N1_N3_N7_N5 = nodes.getByBoundingBox(X_min - delta, Y_min + delta, Z_min + delta, X_min + delta, Y_max - delta, Z_max - delta)
    nodes_N2_N4_N8_N6 = nodes.getByBoundingBox(X_max - delta, Y_min + delta, Z_min + delta, X_max + delta, Y_max - delta, Z_max - delta)
    #Y
    nodes_N1_N2_N6_N5 = nodes.getByBoundingBox(X_min + delta, Y_min - delta, Z_min + delta, X_max - delta, Y_min + delta, Z_max - delta)
    nodes_N3_N4_N8_N7 = nodes.getByBoundingBox(X_min + delta, Y_max - delta, Z_min + delta, X_max - delta, Y_max + delta, Z_max - delta)
    #Z
    nodes_N1_N2_N4_N3 = nodes.getByBoundingBox(X_min + delta, Y_min + delta, Z_min - delta, X_max - delta, Y_max - delta, Z_min + delta)
    nodes_N5_N6_N8_N7 = nodes.getByBoundingBox(X_min + delta, Y_min + delta, Z_max - delta, X_max - delta, Y_max - delta, Z_max + delta)

    set_surf_X_min = mdb.models[model_name].rootAssembly.Set(name='surf_X_min', nodes=nodes_N1_N3_N7_N5)
    set_surf_X_max = mdb.models[model_name].rootAssembly.Set(name='surf_X_max', nodes=nodes_N2_N4_N8_N6)
    X_surf_pairs = ('surf_X_min', 'surf_X_max')

    set_surf_Y_min = mdb.models[model_name].rootAssembly.Set(name='surf_Y_min', nodes=nodes_N1_N2_N6_N5)
    set_surf_Y_max = mdb.models[model_name].rootAssembly.Set(name='surf_Y_max', nodes=nodes_N3_N4_N8_N7)
    Y_surf_pairs = ('surf_Y_min', 'surf_Y_max')

    set_surf_Z_min = mdb.models[model_name].rootAssembly.Set(name='surf_Z_min', nodes=nodes_N1_N2_N4_N3)
    set_surf_Z_max = mdb.models[model_name].rootAssembly.Set(name='surf_Z_max', nodes=nodes_N5_N6_N8_N7)
    Z_surf_pairs = ('surf_Z_min', 'surf_Z_max')

    return X_surf_pairs, Y_surf_pairs, Z_surf_pairs

def set_PBC(model_name, set_RP, X_pairs, Y_pairs, Z_pairs, set_edges, set_vertices):
    direction_X = 1;
    direction_Y = 2;
    direction_Z = 3;

    # Set surface node sets
    for i in range(len(X_pairs)):
        mdb.models[model_name].Equation(name='Eq_surfaces_X_in_dir_X' + str(i),
                                        terms = ((1, X_pairs[i][1], direction_X),
                                                 (-1, X_pairs[i][0], direction_X),
                                                 (-1, set_RP[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_surfaces_X_in_dir_Y' + str(i),
                                        terms = ((1, X_pairs[i][1], direction_Y),
                                                 (-1, X_pairs[i][0], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_surfaces_X_in_dir_Z' + str(i),
                                        terms = ((1, X_pairs[i][1], direction_Z),
                                                 (-1, X_pairs[i][0], direction_Z)))
        
    for i in range(len(Y_pairs)):
        mdb.models[model_name].Equation(name='Eq_surfaces_Y_in_dir_X' + str(i),
                                        terms = ((1, Y_pairs[i][1], direction_X),
                                                 (-1, Y_pairs[i][0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_surfaces_Y_in_dir_Y' + str(i),
                                        terms = ((1, Y_pairs[i][1], direction_Y),
                                                 (-1, Y_pairs[i][0], direction_Y),
                                                 (-1, set_RP[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_surfaces_Y_in_dir_Z' + str(i),
                                        terms = ((1, Y_pairs[i][1], direction_Z),
                                                 (-1, Y_pairs[i][0], direction_Z)))
    
    for i in range(len(Z_pairs)):
        mdb.models[model_name].Equation(name='Eq_surfaces_Z_in_dir_X' + str(i),
                                        terms = ((1, Z_pairs[i][1], direction_X),
                                                 (-1, Z_pairs[i][0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_surfaces_Z_in_dir_Y' + str(i),
                                        terms = ((1, Z_pairs[i][1], direction_Y),
                                                 (-1, Z_pairs[i][0], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_surfaces_Z_in_dir_Z' + str(i),
                                        terms = ((1, Z_pairs[i][1], direction_Z),
                                                 (-1, Z_pairs[i][0], direction_Z),
                                                 (-1, set_RP[2], direction_Z)))
        
    # Set edges sets
    #   0        1        2        3        4        5        6        7        8        9        10       11
    #'N1_N2', 'N1_N3', 'N1_N5', 'N2_N4', 'N2_N6', 'N3_N4', 'N3_N7', 'N4_N8', 'N5_N6', 'N5_N7', 'N6_N8', 'N7_N8'
    for i in range(len(set_edges[10])):
        mdb.models[model_name].Equation(name='Eq_edges_N6_8_to_N5_7_in_dir_X' + str(i),
                                        terms = ((1, set_edges[10][i], direction_X),
                                                (-1, set_edges[9][i], direction_X),
                                                (-1, set_RP[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N6_8_to_N5_7_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[10][i], direction_Y),
                                                (-1, set_edges[9][i], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N6_8_to_N5_7_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[10][i], direction_Z),
                                                (-1, set_edges[9][i], direction_Z)))
        
    for i in range(len(set_edges[9])):
        mdb.models[model_name].Equation(name='Eq_edges_N5_7_to_N1_3_in_dir_X' + str(i),
                                        terms = ((1, set_edges[9][i], direction_X),
                                                 (-1, set_edges[1][i], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N5_7_to_N1_3_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[9][i], direction_Y),
                                                 (-1, set_edges[1][i], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N5_7_to_N1_3_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[9][i], direction_Z),
                                                 (-1, set_edges[1][i], direction_Z),
                                                 (-1, set_RP[2], direction_Z)))

    for i in range(len(set_edges[1])):
        mdb.models[model_name].Equation(name='Eq_edges_N1_3_to_N2_4_in_dir_X' + str(i),
                                        terms = ((1, set_edges[1][i], direction_X),
                                                 (-1, set_edges[3][i], direction_X),
                                                 (1, set_RP[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N1_3_to_N2_4_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[1][i], direction_Y),
                                                 (-1, set_edges[3][i], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N1_3_to_N2_4_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[1][i], direction_Z),
                                                 (-1, set_edges[3][i], direction_Z)))

    for i in range(len(set_edges[4])):
        mdb.models[model_name].Equation(name='Eq_edges_N2_6_to_N1_5_in_dir_X' + str(i),
                                        terms = ((1, set_edges[4][i], direction_X),
                                                 (-1, set_edges[2][i], direction_X),
                                                 (-1, set_RP[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N2_6_to_N1_5_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[4][i], direction_Y),
                                                 (-1, set_edges[2][i], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N2_6_to_N1_5_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[4][i], direction_Z),
                                                 (-1, set_edges[2][i], direction_Z)))

    for i in range(len(set_edges[2])):
        mdb.models[model_name].Equation(name='Eq_edges_N1_5_to_N3_7_in_dir_X' + str(i),
                                        terms = ((1, set_edges[2][i], direction_X),
                                                 (-1, set_edges[6][i], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N1_5_to_N3_7_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[2][i], direction_Y),
                                                 (-1, set_edges[6][i], direction_Y),
                                                 (1, set_RP[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N1_5_to_N3_7_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[2][i], direction_Z),
                                                 (-1, set_edges[6][i], direction_Z)))

    for i in range(len(set_edges[6])):
        mdb.models[model_name].Equation(name='Eq_edges_N3_7_to_N4_8_in_dir_X' + str(i),
                                        terms = ((1, set_edges[6][i], direction_X),
                                                 (-1, set_edges[7][i], direction_X),
                                                 (1, set_RP[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N3_7_to_N4_8_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[6][i], direction_Y),
                                                 (-1, set_edges[7][i], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N3_7_to_N4_8_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[6][i], direction_Z),
                                                 (-1, set_edges[7][i], direction_Z)))

    for i in range(len(set_edges[8])):
        mdb.models[model_name].Equation(name='Eq_edges_N5_6_to_N1_2_in_dir_X' + str(i),
                                        terms = ((1, set_edges[8][i], direction_X),
                                                 (-1, set_edges[0][i], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N5_6_to_N1_2_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[8][i], direction_Y),
                                                 (-1, set_edges[0][i], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N5_6_to_N1_2_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[8][i], direction_Z),
                                                 (-1, set_edges[0][i], direction_Z),
                                                 (-1, set_RP[2], direction_Z)))

    for i in range(len(set_edges[0])):
        mdb.models[model_name].Equation(name='Eq_edges_N1_2_to_N3_4_in_dir_X' + str(i),
                                        terms = ((1, set_edges[0][i], direction_X),
                                                 (-1, set_edges[5][i], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N1_2_to_N3_4_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[0][i], direction_Y),
                                                 (-1, set_edges[5][i], direction_Y),
                                                 (1, set_RP[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N1_2_to_N3_4_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[0][i], direction_Z),
                                                 (-1, set_edges[5][i], direction_Z)))

    for i in range(len(set_edges[5])):
        mdb.models[model_name].Equation(name='Eq_edges_N3_4_to_N7_8_in_dir_X' + str(i),
                                        terms = ((1, set_edges[5][i], direction_X),
                                                 (-1, set_edges[11][i], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N3_4_to_N7_8_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[5][i], direction_Y),
                                                 (-1, set_edges[11][i], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N3_4_to_N7_8_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[5][i], direction_Z),
                                                 (-1, set_edges[11][i], direction_Z),
                                                 (1, set_RP[2], direction_Z)))

    # Set vertices PBC
    mdb.models[model_name].Equation(name='Eq_vertices_N1_to_N5_in_dir_X',
                                        terms = ((1, set_vertices[0], direction_X),
                                                 (-1, set_vertices[4], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N1_to_N5_in_dir_Y',
                                        terms = ((1, set_vertices[0], direction_Y),
                                                 (-1, set_vertices[4], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N1_to_N5_in_dir_Z',
                                        terms = ((1, set_vertices[0], direction_Z),
                                                 (-1, set_vertices[4], direction_Z),
                                                 (1, set_RP[2], direction_Z)))

    mdb.models[model_name].Equation(name='Eq_vertices_N5_to_N7_in_dir_X',
                                        terms = ((1, set_vertices[4], direction_X),
                                                 (-1, set_vertices[6], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N5_to_N7_in_dir_Y',
                                        terms = ((1, set_vertices[4], direction_Y),
                                                 (-1, set_vertices[6], direction_Y),
                                                 (1, set_RP[1], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N5_to_N7_in_dir_Z',
                                        terms = ((1, set_vertices[4], direction_Z),
                                                 (-1, set_vertices[6], direction_Z)))
    
    mdb.models[model_name].Equation(name='Eq_vertices_N7_to_N8_in_dir_X',
                                        terms = ((1, set_vertices[6], direction_X),
                                                 (-1, set_vertices[7], direction_X),
                                                 (1, set_RP[0], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N7_to_N8_in_dir_Y',
                                        terms = ((1, set_vertices[6], direction_Y),
                                                 (-1, set_vertices[7], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N7_to_N8_in_dir_Z',
                                        terms = ((1, set_vertices[6], direction_Z),
                                                 (-1, set_vertices[7], direction_Z)))
    
    
    mdb.models[model_name].Equation(name='Eq_vertices_N8_to_N4_in_dir_X',
                                        terms = ((1, set_vertices[7], direction_X),
                                                 (-1, set_vertices[3], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N8_to_N4_in_dir_Y',
                                        terms = ((1, set_vertices[7], direction_Y),
                                                 (-1, set_vertices[3], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N8_to_N4_in_dir_Z',
                                        terms = ((1, set_vertices[7], direction_Z),
                                                 (-1, set_vertices[3], direction_Z),
                                                 (-1, set_RP[2], direction_Z)))
    
    mdb.models[model_name].Equation(name='Eq_vertices_N4_to_N2_in_dir_X',
                                        terms = ((1, set_vertices[3], direction_X),
                                                 (-1, set_vertices[1], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N4_to_N2_in_dir_Y',
                                        terms = ((1, set_vertices[3], direction_Y),
                                                 (-1, set_vertices[1], direction_Y),
                                                 (-1, set_RP[1], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N4_to_N2_in_dir_Z',
                                        terms = ((1, set_vertices[3], direction_Z),
                                                 (-1, set_vertices[1], direction_Z)))
    
    
    mdb.models[model_name].Equation(name='Eq_vertices_N2_to_N6_in_dit_X',
                                        terms = ((1, set_vertices[1], direction_X),
                                                 (-1, set_vertices[5], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N2_to_N6_in_dit_Y',
                                        terms = ((1, set_vertices[1], direction_Y),
                                                 (-1, set_vertices[5], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N2_to_N6_in_dit_Z',
                                        terms = ((1, set_vertices[1], direction_Z),
                                                 (-1, set_vertices[5], direction_Z),
                                                 (1, set_RP[2], direction_Z)))

    mdb.models[model_name].Equation(name='Eq_vertices_N6_to_N3_in_dir_X',
                                        terms = ((1, set_vertices[5], direction_X),
                                                 (-1, set_vertices[2], direction_X),
                                                 (-1, set_RP[0], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N6_to_N3_in_dir_Y',
                                        terms = ((1, set_vertices[5], direction_Y),
                                                 (-1, set_vertices[2], direction_Y),
                                                 (1, set_RP[1], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N6_to_N3_in_dir_Z',
                                        terms = ((1, set_vertices[5], direction_Z),
                                                 (-1, set_vertices[2], direction_Z),
                                                 (-1, set_RP[2], direction_Z)))
 
def set_PBC_G(model_name, set_RP, X_pairs, Y_pairs, Z_pairs, set_edges, set_vertices, X_pattern, Y_pattern, Z_pattern):
    # Definition of the PBC patterns for equations for X, Y and Z -directions
    '''
    X_pattern = (set_RP[0], set_RP[3], set_RP[4])
    Y_pattern = (set_RP[3], set_RP[1], set_RP[5])
    Z_pattern = (set_RP[4], set_RP[5], set_RP[2])
    '''
    direction_X = 1;
    direction_Y = 2;
    direction_Z = 3;

    # Set surface node sets
    for i in range(len(X_pairs)):
        mdb.models[model_name].Equation(name='Eq_surfaces_X_in_dir_X' + str(i),
                                        terms = ((1, X_pairs[i][1], direction_X),
                                                 (-1, X_pairs[i][0], direction_X),
                                                 (-1, X_pattern[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_surfaces_X_in_dir_Y' + str(i),
                                        terms = ((1, X_pairs[i][1], direction_Y),
                                                 (-1, X_pairs[i][0], direction_Y),
                                                 (-1, X_pattern[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_surfaces_X_in_dir_Z' + str(i),
                                        terms = ((1, X_pairs[i][1], direction_Z),
                                                 (-1, X_pairs[i][0], direction_Z),
                                                 (-1, X_pattern[2], direction_Z)))
      
    for i in range(len(Y_pairs)):
        mdb.models[model_name].Equation(name='Eq_surfaces_Y_in_dir_X' + str(i),
                                        terms = ((1, Y_pairs[i][1], direction_X),
                                                 (-1, Y_pairs[i][0], direction_X),
                                                 (-1, Y_pattern[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_surfaces_Y_in_dir_Y' + str(i),
                                        terms = ((1, Y_pairs[i][1], direction_Y),
                                                 (-1, Y_pairs[i][0], direction_Y),
                                                 (-1, Y_pattern[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_surfaces_Y_in_dir_Z' + str(i),
                                        terms = ((1, Y_pairs[i][1], direction_Z),
                                                 (-1, Y_pairs[i][0], direction_Z),
                                                 (-1, Y_pattern[2], direction_Z)))
    
    for i in range(len(Z_pairs)):
        mdb.models[model_name].Equation(name='Eq_surfaces_Z_in_dir_X' + str(i),
                                        terms = ((1, Z_pairs[i][1], direction_X),
                                                 (-1, Z_pairs[i][0], direction_X),
                                                 (-1, Z_pattern[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_surfaces_Z_in_dir_Y' + str(i),
                                        terms = ((1, Z_pairs[i][1], direction_Y),
                                                 (-1, Z_pairs[i][0], direction_Y),
                                                 (-1, Z_pattern[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_surfaces_Z_in_dir_Z' + str(i),
                                        terms = ((1, Z_pairs[i][1], direction_Z),
                                                 (-1, Z_pairs[i][0], direction_Z),
                                                 (-1, Z_pattern[2], direction_Z)))
    
    # Set edges sets
    #   0        1        2        3        4        5        6        7        8        9        10       11
    #'N1_N2', 'N1_N3', 'N1_N5', 'N2_N4', 'N2_N6', 'N3_N4', 'N3_N7', 'N4_N8', 'N5_N6', 'N5_N7', 'N6_N8', 'N7_N8'
    for i in range(len(set_edges[7])):
        mdb.models[model_name].Equation(name='Eq_edges_N4_8_to_N3_7_in_dir_X' + str(i),
                                        terms = ((1, set_edges[7][i], direction_X),
                                                (-1, set_edges[6][i], direction_X),
                                                (-1, X_pattern[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N4_8_to_N3_7_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[7][i], direction_Y),
                                                (-1, set_edges[6][i], direction_Y),
                                                (-1, X_pattern[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N4_8_to_N3_7_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[7][i], direction_Z),
                                                (-1, set_edges[6][i], direction_Z),
                                                (-1, X_pattern[2], direction_Z)))
    
    for i in range(len(set_edges[6])):
        mdb.models[model_name].Equation(name='Eq_edges_N3_7_to_N1_5_in_dir_X' + str(i),
                                        terms = ((1, set_edges[6][i], direction_X),
                                                 (-1, set_edges[2][i], direction_X),
                                                 (-1, Y_pattern[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N3_7_to_N1_5_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[6][i], direction_Y),
                                                 (-1, set_edges[2][i], direction_Y),
                                                 (-1, Y_pattern[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N3_7_to_N1_5_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[6][i], direction_Z),
                                                 (-1, set_edges[2][i], direction_Z),
                                                 (-1, Y_pattern[2], direction_Z)))
    
    for i in range(len(set_edges[2])):
        mdb.models[model_name].Equation(name='Eq_edges_N1_5_to_N2_6_in_dir_X' + str(i),
                                        terms = ((1, set_edges[2][i], direction_X),
                                                 (-1, set_edges[4][i], direction_X),
                                                 (1, X_pattern[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N1_5_to_N2_6_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[2][i], direction_Y),
                                                 (-1, set_edges[4][i], direction_Y),
                                                 (1, X_pattern[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N1_5_to_N2_6_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[2][i], direction_Z),
                                                 (-1, set_edges[4][i], direction_Z),
                                                 (1, X_pattern[2], direction_Z)))
    
    for i in range(len(set_edges[10])):
        mdb.models[model_name].Equation(name='Eq_edges_N6_8_to_N5_7_in_dir_X' + str(i),
                                        terms = ((1, set_edges[10][i], direction_X),
                                                 (-1, set_edges[9][i], direction_X),
                                                 (-1, X_pattern[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N6_8_to_N5_7_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[10][i], direction_Y),
                                                 (-1, set_edges[9][i], direction_Y),
                                                 (-1, X_pattern[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N6_8_to_N5_7_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[10][i], direction_Z),
                                                 (-1, set_edges[9][i], direction_Z),
                                                 (-1, X_pattern[2], direction_Z)))
    
    for i in range(len(set_edges[9])):
        mdb.models[model_name].Equation(name='Eq_edges_N5_7_to_N1_3_in_dir_X' + str(i),
                                        terms = ((1, set_edges[9][i], direction_X),
                                                 (-1, set_edges[1][i], direction_X),
                                                 (-1, Z_pattern[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N5_7_to_N1_3_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[9][i], direction_Y),
                                                 (-1, set_edges[1][i], direction_Y),
                                                 (-1, Z_pattern[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N5_7_to_N1_3_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[9][i], direction_Z),
                                                 (-1, set_edges[1][i], direction_Z),
                                                 (-1, Z_pattern[2], direction_Z)))
    
    for i in range(len(set_edges[1])):
        mdb.models[model_name].Equation(name='Eq_edges_N1_3_to_N2_4_in_dir_X' + str(i),
                                        terms = ((1, set_edges[1][i], direction_X),
                                                 (-1, set_edges[3][i], direction_X),
                                                 (1, X_pattern[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N1_3_to_N2_4_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[1][i], direction_Y),
                                                 (-1, set_edges[3][i], direction_Y),
                                                 (1, X_pattern[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N1_3_to_N2_4_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[1][i], direction_Z),
                                                 (-1, set_edges[3][i], direction_Z),
                                                 (1, X_pattern[2], direction_Z)))
    
    for i in range(len(set_edges[11])):
        mdb.models[model_name].Equation(name='Eq_edges_N7_8_to_N5_6_in_dir_X' + str(i),
                                        terms = ((1, set_edges[11][i], direction_X),
                                                 (-1, set_edges[8][i], direction_X),
                                                 (-1, Y_pattern[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N7_8_to_N5_6_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[11][i], direction_Y),
                                                 (-1, set_edges[8][i], direction_Y),
                                                 (-1, Y_pattern[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N7_8_to_N5_6_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[11][i], direction_Z),
                                                 (-1, set_edges[8][i], direction_Z),
                                                 (-1, Y_pattern[2], direction_Z)))
    
    for i in range(len(set_edges[8])):
        mdb.models[model_name].Equation(name='Eq_edges_N5_6_to_N1_2_in_dir_X' + str(i),
                                        terms = ((1, set_edges[8][i], direction_X),
                                                 (-1, set_edges[0][i], direction_X),
                                                 (-1, Z_pattern[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N5_6_to_N1_2_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[8][i], direction_Y),
                                                 (-1, set_edges[0][i], direction_Y),
                                                 (-1, Z_pattern[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N5_6_to_N1_2_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[8][i], direction_Z),
                                                 (-1, set_edges[0][i], direction_Z),
                                                 (-1, Z_pattern[2], direction_Z)))
    
    for i in range(len(set_edges[0])):
        mdb.models[model_name].Equation(name='Eq_edges_N1_2_to_N3_4_in_dir_X' + str(i),
                                        terms = ((1, set_edges[0][i], direction_X),
                                                 (-1, set_edges[5][i], direction_X),
                                                 (1, Y_pattern[0], direction_X)))
        mdb.models[model_name].Equation(name='Eq_edges_N1_2_to_N3_4_in_dir_Y' + str(i),
                                        terms = ((1, set_edges[0][i], direction_Y),
                                                 (-1, set_edges[5][i], direction_Y),
                                                 (1, Y_pattern[1], direction_Y)))
        mdb.models[model_name].Equation(name='Eq_edges_N1_2_to_N3_4_in_dir_Z' + str(i),
                                        terms = ((1, set_edges[0][i], direction_Z),
                                                 (-1, set_edges[5][i], direction_Z),
                                                 (1, Y_pattern[2], direction_Z)))
    
    # Set vertices PBC
    mdb.models[model_name].Equation(name='Eq_vertices_N5_to_N7_in_dir_X',
                                        terms = ((1, set_vertices[4], direction_X),
                                                 (-1, set_vertices[6], direction_X),
                                                 (1, Y_pattern[0], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N5_to_N7_in_dir_Y',
                                        terms = ((1, set_vertices[4], direction_Y),
                                                 (-1, set_vertices[6], direction_Y),
                                                 (1, Y_pattern[1], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N5_to_N7_in_dir_Z',
                                        terms = ((1, set_vertices[4], direction_Z),
                                                 (-1, set_vertices[6], direction_Z),
                                                 (1, Y_pattern[2], direction_Z)))
    
    mdb.models[model_name].Equation(name='Eq_vertices_N7_to_N3_in_dir_X',
                                        terms = ((1, set_vertices[6], direction_X),
                                                 (-1, set_vertices[2], direction_X),
                                                 (-1, Z_pattern[0], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N7_to_N3_in_dir_Y',
                                        terms = ((1, set_vertices[6], direction_Y),
                                                 (-1, set_vertices[2], direction_Y),
                                                 (-1, Z_pattern[1], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N7_to_N3_in_dir_Z',
                                        terms = ((1, set_vertices[6], direction_Z),
                                                 (-1, set_vertices[2], direction_Z),
                                                 (-1, Z_pattern[2], direction_Z)))
    
    mdb.models[model_name].Equation(name='Eq_vertices_N3_to_N4_in_dir_X',
                                        terms = ((1, set_vertices[2], direction_X),
                                                 (-1, set_vertices[3], direction_X),
                                                 (1, X_pattern[0], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N3_to_N4_in_dir_Y',
                                        terms = ((1, set_vertices[2], direction_Y),
                                                 (-1, set_vertices[3], direction_Y),
                                                 (1, X_pattern[1], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N3_to_N4_in_dir_Z',
                                        terms = ((1, set_vertices[2], direction_Z),
                                                 (-1, set_vertices[3], direction_Z),
                                                 (1, X_pattern[2], direction_Z)))
    
    mdb.models[model_name].Equation(name='Eq_vertices_N4_to_N2_in_dir_X',
                                        terms = ((1, set_vertices[3], direction_X),
                                                 (-1, set_vertices[1], direction_X),
                                                 (-1, Y_pattern[0], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N4_to_N2_in_dir_Y',
                                        terms = ((1, set_vertices[3], direction_Y),
                                                 (-1, set_vertices[1], direction_Y),
                                                 (-1, Y_pattern[1], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N4_to_N2_in_dir_Z',
                                        terms = ((1, set_vertices[3], direction_Z),
                                                 (-1, set_vertices[1], direction_Z),
                                                 (-1, Y_pattern[2], direction_Z)))
    
    mdb.models[model_name].Equation(name='Eq_vertices_N2_to_N6_in_dir_X',
                                        terms = ((1, set_vertices[1], direction_X),
                                                 (-1, set_vertices[5], direction_X),
                                                 (1, Z_pattern[0], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N2_to_N6_in_dir_Y',
                                        terms = ((1, set_vertices[1], direction_Y),
                                                 (-1, set_vertices[5], direction_Y),
                                                 (1, Z_pattern[1], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N2_to_N6_in_dir_Z',
                                        terms = ((1, set_vertices[1], direction_Z),
                                                 (-1, set_vertices[5], direction_Z),
                                                 (1, Z_pattern[2], direction_Z)))
    
    mdb.models[model_name].Equation(name='Eq_vertices_N6_to_N8_in_dir_X',
                                        terms = ((1, set_vertices[5], direction_X),
                                                 (-1, set_vertices[7], direction_X),
                                                 (1, Y_pattern[0], direction_X)))
    mdb.models[model_name].Equation(name='Eq_vertices_N6_to_N8_in_dir_Y',
                                        terms = ((1, set_vertices[5], direction_Y),
                                                 (-1, set_vertices[7], direction_Y),
                                                 (1, Y_pattern[1], direction_Y)))
    mdb.models[model_name].Equation(name='Eq_vertices_N6_to_N8_in_dir_Z',
                                        terms = ((1, set_vertices[5], direction_Z),
                                                 (-1, set_vertices[7], direction_Z),
                                                 (1, Y_pattern[2], direction_Z)))
    
    mdb.models[model_name].Equation(name='Eq_vertices_N8_to_N1_in_dir_X',
                                        terms = ((1, set_vertices[7], direction_X),
                                                 (-1, set_vertices[0], direction_X),
                                                 (-1, X_pattern[0], direction_X),
                                                 (-1, X_pattern[1], direction_X),
                                                 (-1, X_pattern[2], direction_X)))
    
    mdb.models[model_name].Equation(name='Eq_vertices_N8_to_N1_in_dir_Y',
                                        terms = ((1, set_vertices[7], direction_Y),
                                                 (-1, set_vertices[0], direction_Y),
                                                 (-1, Y_pattern[1], direction_Y),
                                                 (-1, Y_pattern[0], direction_Y),
                                                 (-1, Y_pattern[2], direction_Y)))
    
    mdb.models[model_name].Equation(name='Eq_vertices_N8_to_N1_in_dir_Z',
                                        terms = ((1, set_vertices[7], direction_Z),
                                                 (-1, set_vertices[0], direction_Z),
                                                 (-1, Z_pattern[2], direction_Z),
                                                 (-1, Z_pattern[0], direction_Z),
                                                 (-1, Z_pattern[1], direction_Z)))

def set_BC(model_name, set_RP, step_name, strain):
    region_RP_X = mdb.models[model_name].rootAssembly.sets[set_RP[0]]
    region_RP_Y = mdb.models[model_name].rootAssembly.sets[set_RP[1]]
    region_RP_Z = mdb.models[model_name].rootAssembly.sets[set_RP[2]]
    region_RP_4 = mdb.models[model_name].rootAssembly.sets[set_RP[3]]
    region_RP_5 = mdb.models[model_name].rootAssembly.sets[set_RP[4]]
    region_RP_6 = mdb.models[model_name].rootAssembly.sets[set_RP[5]]
    mdb.models[model_name].StaticStep(name=step_name, previous='Initial')
    step_name = step_name
    if step_name == 'Step_E11':
        mdb.models[model_name].DisplacementBC(name='BC_X', createStepName=step_name, region=region_RP_X,
                                              u1=strains[0], u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='X_disp', localCsys=None)
        mdb.models[model_name].DisplacementBC(name='BC_Y', createStepName=step_name, region=region_RP_Y,
                                              u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='Y_disp', localCsys=None)
        mdb.models[model_name].DisplacementBC(name='BC_Z', createStepName=step_name, region=region_RP_Z,
                                              u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='Z_disp', localCsys=None)
    elif step_name == 'Step_E22':
        mdb.models[model_name].DisplacementBC(name='BC_X', createStepName=step_name, region=region_RP_X,
                                              u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='X_disp', localCsys=None)
        mdb.models[model_name].DisplacementBC(name='BC_Y', createStepName=step_name, region=region_RP_Y,
                                              u1=UNSET, u2=strains[1], u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='Y_disp', localCsys=None)
        mdb.models[model_name].DisplacementBC(name='BC_Z', createStepName=step_name, region=region_RP_Z,
                                              u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='Z_disp', localCsys=None)
    elif step_name == 'Step_E33':
        mdb.models[model_name].DisplacementBC(name='BC_X', createStepName=step_name, region=region_RP_X,
                                              u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='X_disp', localCsys=None)
        mdb.models[model_name].DisplacementBC(name='BC_Y', createStepName=step_name, region=region_RP_Y,
                                              u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='Y_disp', localCsys=None)
        mdb.models[model_name].DisplacementBC(name='BC_Z', createStepName=step_name, region=region_RP_Z,
                                              u1=UNSET, u2=UNSET, u3=strains[2], ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='Z_disp', localCsys=None)
    elif step_name == 'Step_G12':
        mdb.models[model_name].DisplacementBC(name='BC_X', createStepName=step_name, region=region_RP_4,
                                              u1=strains[0], u2=strains[1], u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='X_disp', localCsys=None)
        mdb.models[model_name].DisplacementBC(name='BC_Y', createStepName=step_name, region=region_RP_5,
                                              u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='Y_disp', localCsys=None)
        mdb.models[model_name].DisplacementBC(name='BC_Z', createStepName=step_name, region=region_RP_6,
                                              u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='Z_disp', localCsys=None)
        
    elif step_name == 'Step_G23':
        mdb.models[model_name].DisplacementBC(name='BC_X', createStepName=step_name, region=region_RP_4,
                                              u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='X_disp', localCsys=None)
        mdb.models[model_name].DisplacementBC(name='BC_Y', createStepName=step_name, region=region_RP_5,
                                              u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='Y_disp', localCsys=None)
        mdb.models[model_name].DisplacementBC(name='BC_Z', createStepName=step_name, region=region_RP_6,
                                              u1=UNSET, u2=strains[1], u3=strains[2], ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='Z_disp', localCsys=None)
    elif step_name == 'Step_G13':
        mdb.models[model_name].DisplacementBC(name='BC_X', createStepName=step_name, region=region_RP_4,
                                              u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='X_disp', localCsys=None)
        mdb.models[model_name].DisplacementBC(name='BC_Y', createStepName=step_name, region=region_RP_5,
                                              u1=strains[0], u2=UNSET, u3=strains[2], ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='Y_disp', localCsys=None)
        mdb.models[model_name].DisplacementBC(name='BC_Z', createStepName=step_name, region=region_RP_6,
                                              u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                                              amplitude=UNSET, distributionType=UNIFORM, fieldName='Z_disp', localCsys=None)

def create_Job(model_name, job_name):
    current_job = mdb.Job(name=job_name, model=model_name, description='RUC', type=ANALYSIS, 
                          atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                          memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                          explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                          modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                          scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
                          multiprocessingMode=DEFAULT, numCpus=num_CPUs, numDomains=num_Domains, numGPUs=0)
    return current_job

def get_RUC_sizes_and_strains_from_FO(current_odb, step, instance_name, X_pairs, Y_pairs, Z_pairs):
    # Read U at the last increment, the same increment used for RF.
    instance = current_odb.rootAssembly.instances[instance_name.upper()]
    # Store the original coordinates by node label.
    coordinates = {}
    for i in range(len(instance.nodes)):
        node = instance.nodes[i]
        coordinates[node.label] = node.coordinates
    # Match values by node label, not by their order in the ODB.
    current_frame = step.frames[-1]
    U_values = current_frame.fieldOutputs['U'].getSubset(region=instance).values
    displacements = {}
    for i in range(len(U_values)):
        displacements[U_values[i].nodeLabel] = U_values[i].data

    Xs = []
    X_deltas = []
    Ys = []
    Y_deltas = []
    Zs = []
    Z_deltas = []
    
    for i in range(len(X_pairs)):
        node_label_min = current_odb.rootAssembly.nodeSets[X_pairs[i][0].upper()].nodes[0][0].label
        node_label_max = current_odb.rootAssembly.nodeSets[X_pairs[i][1].upper()].nodes[0][0].label

        U_delta = displacements[node_label_max][0] - displacements[node_label_min][0]
        C_delta = coordinates[node_label_max][0] - coordinates[node_label_min][0]
        Xs.append(C_delta)
        X_deltas.append(U_delta)
    X = sum(Xs) / len(Xs)
    X_delta = sum(X_deltas) / len(X_deltas)

    for i in range(len(Y_pairs)):
        node_label_min = current_odb.rootAssembly.nodeSets[Y_pairs[i][0].upper()].nodes[0][0].label
        node_label_max = current_odb.rootAssembly.nodeSets[Y_pairs[i][1].upper()].nodes[0][0].label

        U_delta = displacements[node_label_max][1] - displacements[node_label_min][1]
        C_delta = coordinates[node_label_max][1] - coordinates[node_label_min][1]
        Ys.append(C_delta)
        Y_deltas.append(U_delta)
    Y = sum(Ys) / len(Ys)
    Y_delta = sum(Y_deltas) / len(Y_deltas)

    for i in range(len(Z_pairs)):
        node_label_min = current_odb.rootAssembly.nodeSets[Z_pairs[i][0].upper()].nodes[0][0].label
        node_label_max = current_odb.rootAssembly.nodeSets[Z_pairs[i][1].upper()].nodes[0][0].label

        U_delta = displacements[node_label_max][2] - displacements[node_label_min][2]
        C_delta = coordinates[node_label_max][2] - coordinates[node_label_min][2]
        Zs.append(C_delta)
        Z_deltas.append(U_delta)
    Z = sum(Zs) / len(Zs)
    Z_delta = sum(Z_deltas) / len(Z_deltas)
    return X, X_delta, Y, Y_delta, Z, Z_delta

def get_RUC_sizes_from_mesh(current_odb, step, instance_name, X_pairs, Y_pairs, Z_pairs):
    # ODB node coordinates are the original, undeformed mesh coordinates.
    instance = current_odb.rootAssembly.instances[instance_name.upper()]
    # Store the original coordinates by node label.
    coordinates = {}
    for i in range(len(instance.nodes)):
        node = instance.nodes[i]
        coordinates[node.label] = node.coordinates
    Xs = []
    Ys = []
    Zs = []
    
    for i in range(len(X_pairs)):
        node_label_min = current_odb.rootAssembly.nodeSets[X_pairs[i][0].upper()].nodes[0][0].label
        node_label_max = current_odb.rootAssembly.nodeSets[X_pairs[i][1].upper()].nodes[0][0].label

        C_delta = coordinates[node_label_max][0] - coordinates[node_label_min][0]
        Xs.append(C_delta)
    X = sum(Xs) / len(Xs)
    
    for i in range(len(Y_pairs)):
        node_label_min = current_odb.rootAssembly.nodeSets[Y_pairs[i][0].upper()].nodes[0][0].label
        node_label_max = current_odb.rootAssembly.nodeSets[Y_pairs[i][1].upper()].nodes[0][0].label

        C_delta = coordinates[node_label_max][1] - coordinates[node_label_min][1]
        Ys.append(C_delta)
    Y = sum(Ys) / len(Ys)
    
    for i in range(len(Z_pairs)):
        node_label_min = current_odb.rootAssembly.nodeSets[Z_pairs[i][0].upper()].nodes[0][0].label
        node_label_max = current_odb.rootAssembly.nodeSets[Z_pairs[i][1].upper()].nodes[0][0].label

        C_delta = coordinates[node_label_max][2] - coordinates[node_label_min][2]
        Zs.append(C_delta)
    Z = sum(Zs) / len(Zs)
    
    return X, Y, Z

def get_RP_forces_from_FO(current_odb, current_step, RP_name):
    # Select by the existing RP set, not by history-region ordering.
    region = current_odb.rootAssembly.nodeSets[RP_name.upper()]
    values = current_step.frames[-1].fieldOutputs['RF'].getSubset(region=region).values
    if len(values) != 1:
        raise RuntimeError('Expected one reaction-force value in set ' + RP_name)
    RP_X = values[0].data[0]
    RP_Y = values[0].data[1]
    RP_Z = values[0].data[2]

    return RP_X, RP_Y, RP_Z


def get_RP_forces_from_FO_E(current_odb, current_step):
    X_forces = get_RP_forces_from_FO(current_odb, current_step, 'RP_X')
    Y_forces = get_RP_forces_from_FO(current_odb, current_step, 'RP_Y')
    Z_forces = get_RP_forces_from_FO(current_odb, current_step, 'RP_Z')
    return X_forces, Y_forces, Z_forces


def create_field_output(model_name, step_name):
    model = mdb.models[model_name]
    # Remove default/inherited histories too, so the node count cannot exceed
    # the history-request limit. Existing field requests are retained.
    for name in list(model.historyOutputRequests.keys()):
        del model.historyOutputRequests[name]
    if 'FO_PBC' in model.fieldOutputRequests:
        del model.fieldOutputRequests['FO_PBC']
    model.FieldOutputRequest(name='FO_PBC', createStepName=step_name,
                             variables=('U', 'RF'), frequency=LAST_INCREMENT)

# 1. Copy a part to assembly
set_part_to_assembly(model_name, part_name, instance_name)

# 2. Find coordinates of RUC vertexes and center
(N1, N2, N3, N4, N5, N6, N7, N8, Nc) = get_RUC_vertices(model_name, instance_name)

if calculate_p:
    props = mdb.models[model_name].parts[part_name].getMassProperties()
    density = props['mass']/((N8[0] - N1[0])*(N8[1] - N1[1])*(N8[2] - N1[2]))
    print('--------------------------------------------------------')
    print('p     : ' + str(density) + ' kg/m^3')


if calculate_E11 or calculate_E22 or calculate_E33 or calculate_G12 or calculate_G13 or calculate_G23:

    # 3. Define strains for x, y and z -directions
    strains = calculate_strains(N1, N8, 0.005)

    # 4. Create sets at vertices
    (set_N1, set_N2, set_N3, set_N4, set_N5, set_N6, set_N7, set_N8, set_vertices) = set_vertices_to_sets(model_name, instance_name, delta, N1, N2, N3, N4, N5, N6, N7, N8)

    # 5. Create node sets for nodes at surfaces
    (X_pairs, Y_pairs, Z_pairs) = get_node_pair_sets_surfaces(model_name, instance_name, delta, N1, N2, N3, N4, N5, N6, N7, N8)
    (X_surf_pairs, Y_surf_pairs, Z_surf_pairs) = get_surface_pairs(model_name, instance_name, delta, N1, N2, N3, N4, N5, N6, N7, N8)

    # 6. Create node sets for nodes at edges
    set_edges = get_node_sets_edges(model_name, instance_name, delta, N1, N2, N3, N4, N5, N6, N7, N8)

    # 7. Set reference points
    set_RP_X = set_RP_point_to_set(model_name, instance_name, (N2[0]*2, Nc[1], Nc[2]), 'RP_X')
    set_RP_Y = set_RP_point_to_set(model_name, instance_name, (Nc[0], N3[1]*2, Nc[2]), 'RP_Y')
    set_RP_Z = set_RP_point_to_set(model_name, instance_name, (Nc[0], Nc[1], N5[2]*2), 'RP_Z')
    set_RP_X_1 = set_RP_point_to_set(model_name, instance_name, (N2[0]*3, Nc[1], Nc[2]), 'RP_X_1')
    set_RP_X_2 = set_RP_point_to_set(model_name, instance_name, (N2[0]*4, Nc[1], Nc[2]), 'RP_X_2')
    set_RP_X_3 = set_RP_point_to_set(model_name, instance_name, (N2[0]*5, Nc[1], Nc[2]), 'RP_X_3')
    set_RPs = ('RP_X', 'RP_Y', 'RP_Z', 'RP_X_1', 'RP_X_2', 'RP_X_3')
    X_pattern = (set_RPs[0], set_RPs[3], set_RPs[4])
    Y_pattern = (set_RPs[3], set_RPs[1], set_RPs[5])
    Z_pattern = (set_RPs[4], set_RPs[5], set_RPs[2])


    if calculate_E11 or calculate_E22 or calculate_E33:
        # 8. Set PBC
        set_PBC(model_name, set_RPs, X_pairs, Y_pairs, Z_pairs, set_edges, set_vertices)
        # 9. Create jobs and calculate Young moduli
        # E11
        if calculate_E11:
            step_name = 'Step_E11'
            # 9.1 Set BC
            step_E = set_BC(model_name, set_RPs, step_name, strains[0])
            # 9.2 Create field outputs for displacements and reaction forces
            create_field_output(model_name, step_name)
            # 9.3 Create and run job
            Job_E = create_Job(model_name, 'Job_E11')
            Job_E.submit()
            Job_E.waitForCompletion()
            if Job_E.status != COMPLETED:
                raise RuntimeError('Job ' + Job_E.name + ' did not complete; check its .dat and .msg files.')
            print('Job E11 is done!')
            # 9.4 Open ODB
            current_odb = session.openOdb(name='Job_E11.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=current_odb)
            session.viewports['Viewport: 1'].makeCurrent()
            current_odb = session.odbs['Job_E11.odb']
            current_step = current_odb.steps[step_name]
            step_0 = current_odb.steps[step_name].frames[0]
            step_1 = current_odb.steps[step_name].frames[-1]
            # 9.5 Calculate RUC sizes and strains, get force at the reference point
            (X_size, X_delta, Y_size, Y_delta, Z_size, Z_delta) = get_RUC_sizes_and_strains_from_FO(current_odb, current_step, instance_name, X_pairs, Y_pairs, Z_pairs)
            (X_forces, Y_forces, Z_forces) = get_RP_forces_from_FO_E(current_odb, current_step)
            # 9.6 Calculate cross-sectional area, axial strains and stress
            cross_section_area = Z_size * Y_size
            axial_strain_x = X_delta / X_size
            force_x = X_forces[0]
            stress_x = force_x / cross_section_area
            # 9.7 Calculate Young modulus andPoisson's ratio
            E11 = stress_x / axial_strain_x
            nu12 = -(Y_delta/Y_size) / (X_delta/X_size)
            nu13 = -(Z_delta/Z_size) / (X_delta/X_size)
            del mdb.models[model_name].steps['Step_E11']

        # E22
        if calculate_E22:
            step_name = 'Step_E22'
            # 9.1 Set BC
            step_E = set_BC(model_name, set_RPs, step_name, strains[0])
            # 9.2 Create field outputs for displacements and reaction forces
            create_field_output(model_name, step_name)
            # 9.3 Create and run job
            Job_E = create_Job(model_name, 'Job_E22')
            Job_E.submit()
            Job_E.waitForCompletion()
            if Job_E.status != COMPLETED:
                raise RuntimeError('Job ' + Job_E.name + ' did not complete; check its .dat and .msg files.')
            print('Job E22 is done!')
            # 9.4 Open ODB
            current_odb = session.openOdb(name='Job_E22.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=current_odb)
            session.viewports['Viewport: 1'].makeCurrent()
            current_odb = session.odbs['Job_E22.odb']
            current_step = current_odb.steps[step_name]
            step_0 = current_odb.steps[step_name].frames[0]
            step_1 = current_odb.steps[step_name].frames[-1]
            # 9.5 Calculate RUC sizes and strains, get force at the reference point
            (X_size, X_delta, Y_size, Y_delta, Z_size, Z_delta) = get_RUC_sizes_and_strains_from_FO(current_odb, current_step, instance_name, X_pairs, Y_pairs, Z_pairs)
            (X_forces, Y_forces, Z_forces) = get_RP_forces_from_FO_E(current_odb, current_step)
            # 9.6 Calculate cross-sectional area, axial strains and stress
            cross_section_area = Z_size * X_size
            axial_strain_y = Y_delta / Y_size
            force_y = Y_forces[1]
            stress_y = force_y / cross_section_area
            # 9.7 Calculate Young modulus andPoisson's ratio
            E22 = stress_y / axial_strain_y
            nu21 = -(X_delta/X_size) / (Y_delta/Y_size)
            nu23 = -(Z_delta/Z_size) / (Y_delta/Y_size)
            del mdb.models[model_name].steps['Step_E22']

        # E33
        if calculate_E33:
            step_name = 'Step_E33'
            # 9.1 Set BC
            step_E = set_BC(model_name, set_RPs, step_name, strains[0])
            # 9.2 Create field outputs for displacements and reaction forces
            create_field_output(model_name, step_name)
            # 9.3 Create and run job
            Job_E = create_Job(model_name, 'Job_E33')
            Job_E.submit()
            Job_E.waitForCompletion()
            if Job_E.status != COMPLETED:
                raise RuntimeError('Job ' + Job_E.name + ' did not complete; check its .dat and .msg files.')
            print('Job E33 is done!')
            # 9.4 Open ODB
            current_odb = session.openOdb(name='Job_E33.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=current_odb)
            session.viewports['Viewport: 1'].makeCurrent()
            current_odb = session.odbs['Job_E33.odb']
            current_step = current_odb.steps[step_name]
            step_0 = current_odb.steps[step_name].frames[0]
            step_1 = current_odb.steps[step_name].frames[-1]
            # 9.5 Calculate RUC sizes and strains, get force at the reference point
            (X_size, X_delta, Y_size, Y_delta, Z_size, Z_delta) = get_RUC_sizes_and_strains_from_FO(current_odb, current_step, instance_name, X_pairs, Y_pairs, Z_pairs)
            (X_forces, Y_forces, Z_forces) = get_RP_forces_from_FO_E(current_odb, current_step)
            # 9.6 Calculate cross-sectional area, axial strains and stress
            cross_section_area = X_size * Y_size
            axial_strain_z = Z_delta / Z_size
            force_z = Z_forces[2]
            stress_z = force_z / cross_section_area
            # 9.7 Calculate Young modulus andPoisson's ratio
            E33 = stress_z / axial_strain_z
            nu31 = -(X_delta/X_size) / (Z_delta/Z_size)
            nu32 = -(Y_delta/Y_size) / (Z_delta/Z_size)
            del mdb.models[model_name].steps['Step_E33']

    if calculate_G12 or calculate_G13 or calculate_G23:
        for i in mdb.models[model_name].constraints.keys():
            del mdb.models[model_name].constraints[i]
        # 8. Set PBC
        set_PBC_G(model_name, set_RPs, X_pairs, Y_pairs, Z_pairs, set_edges, set_vertices, X_pattern, Y_pattern, Z_pattern)

        # 9. Create jobs and calculate shear moduli
        # G12
        if calculate_G12:
            step_name = 'Step_G12'
            # 9.1 Set BC
            step_G12 = set_BC(model_name, set_RPs, step_name, strains)
            # 9.2 Create field outputs for displacements and reaction forces
            create_field_output(model_name, step_name)
            # 9.3 Create and run job
            Job_G12 = create_Job(model_name, 'Job_G12')
            Job_G12.submit()
            Job_G12.waitForCompletion()
            if Job_G12.status != COMPLETED:
                raise RuntimeError('Job ' + Job_G12.name + ' did not complete; check its .dat and .msg files.')
            print('Job G12 is done!')
            # 9.4 Open ODB
            current_odb = session.openOdb(name='Job_G12.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=current_odb)
            session.viewports['Viewport: 1'].makeCurrent()
            current_odb = session.odbs['Job_G12.odb']
            current_step = current_odb.steps[step_name]
            step_0 = current_odb.steps[step_name].frames[0]
            step_1 = current_odb.steps[step_name].frames[-1]
            # 9.5 Calculate RUC sizes and get force at the reference point
            (X_size, Y_size, Z_size) = get_RUC_sizes_from_mesh(current_odb, current_step, instance_name, X_pairs, Y_pairs, Z_pairs)
            (X_force, Y_force, Z_force) = get_RP_forces_from_FO(current_odb, current_step, 'RP_X_1')
            # 9.6 Calculate cross-sectional area, axial strains and stress
            cross_section_area = X_size * Z_size
            axial_strain_x = strains[1] / X_size
            axial_strain_y = strains[0] / Y_size
            stress = X_force / cross_section_area
            # 9.7 Calculate shear modulus
            G_12 = stress / (axial_strain_x + axial_strain_y)
            '''
            # 9.8 Enable optimal point of view to the RUC in the viewport
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
            session.viewports['Viewport: 1'].view.setValues(
            nearPlane=Z_size*10, 
            farPlane=Z_size*11, width=X_size*2, height=Y_size*2,
            cameraPosition=(Nc[0], Nc[1], Z_size*10),
            cameraUpVector=(0, Nc[1], 1),
            cameraTarget=(Nc[0], Nc[1], Nc[2]))
            '''
            del mdb.models[model_name].steps[step_name]

        # G13
        if calculate_G13:
            step_name = 'Step_G13'
            # 9.1 Set BC
            step_G13 = set_BC(model_name, set_RPs, step_name, strains)
            # 9.2 Create field outputs for displacements and reaction forces
            create_field_output(model_name, step_name)
            # 9.3 Create and run job
            Job_G = create_Job(model_name, 'Job_G13')
            Job_G.submit()
            Job_G.waitForCompletion()
            if Job_G.status != COMPLETED:
                raise RuntimeError('Job ' + Job_G.name + ' did not complete; check its .dat and .msg files.')
            print('Job G13 is done!')
            # 9.4 Open ODB
            current_odb = session.openOdb(name='Job_G13.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=current_odb)
            session.viewports['Viewport: 1'].makeCurrent()
            current_odb = session.odbs['Job_G13.odb']
            current_step = current_odb.steps[step_name]
            step_0 = current_odb.steps[step_name].frames[0]
            step_1 = current_odb.steps[step_name].frames[-1]
            # 9.5 Calculate RUC sizes and get force at the reference point
            (X_size, Y_size, Z_size) = get_RUC_sizes_from_mesh(current_odb, current_step, instance_name, X_pairs, Y_pairs, Z_pairs)
            (X_force, Y_force, Z_force) = get_RP_forces_from_FO(current_odb, current_step, 'RP_X_2')
            # 9.6 Calculate cross-sectional area, axial strains and stress
            cross_section_area = X_size * Y_size
            axial_strain_x = strains[2] / X_size
            axial_strain_z = strains[0] / Z_size
            stress = X_force / cross_section_area
            # 9.7 Calculate shear modulus
            G_13 = stress / (axial_strain_x + axial_strain_z)
            '''
            # 9.8 Enable optimal point of view to the RUC in the viewport
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
            session.viewports['Viewport: 1'].view.setValues(
                nearPlane=Y_size*10, 
                farPlane=Y_size*11, width=X_size*2, height=Z_size*2,
                cameraPosition=(Nc[0], -Y_size*10, Nc[2]),
                cameraUpVector=(0, Nc[1], Nc[2]),
                cameraTarget=(Nc[0], Nc[1], Nc[2]))
            '''
            del mdb.models[model_name].steps['Step_G13']

        # G23
        if calculate_G23:
            step_name = 'Step_G23'
            # 9.1 Set BC
            step_G23 = set_BC(model_name, set_RPs, step_name, strains)
            # 9.2 Create field outputs for displacements and reaction forces
            create_field_output(model_name, step_name)
            # 9.3 Create and run job
            Job_G = create_Job(model_name, 'Job_G23')
            Job_G.submit()
            Job_G.waitForCompletion()
            if Job_G.status != COMPLETED:
                raise RuntimeError('Job ' + Job_G.name + ' did not complete; check its .dat and .msg files.')
            print('Job G23 is done!')
            # 9.4 Open ODB
            current_odb = session.openOdb(name='Job_G23.odb')
            session.viewports['Viewport: 1'].setValues(displayedObject=current_odb)
            session.viewports['Viewport: 1'].makeCurrent()
            current_odb = session.odbs['Job_G23.odb']
            current_step = current_odb.steps[step_name]
            step_0 = current_odb.steps[step_name].frames[0]
            step_1 = current_odb.steps[step_name].frames[-1]
            # 9.5 Calculate RUC sizes and get force at the reference point
            (X_size, Y_size, Z_size) = get_RUC_sizes_from_mesh(current_odb, current_step, instance_name, X_pairs, Y_pairs, Z_pairs)
            (X_force, Y_force, Z_force) = get_RP_forces_from_FO(current_odb, current_step, 'RP_X_3')
            # 9.6 Calculate cross-sectional area, axial strains and stress
            cross_section_area = X_size * Y_size
            axial_strain_y = strains[2] / Y_size
            axial_strain_z = strains[1] / Z_size
            stress = Y_force / cross_section_area
            # 9.7 Calculate shear modulus
            G_23 = stress / (axial_strain_y + axial_strain_z)
            print('G23 is: ', str(G_23))
            '''
            # 9.8 Enable optimal point of view to the RUC in the viewport
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
            session.viewports['Viewport: 1'].view.setValues(
                nearPlane=X_size*10, 
                farPlane=X_size*11, width=Y_size*2, height=Z_size*2,
                cameraPosition=(X_size*10, Nc[1], Nc[2]),
                cameraUpVector=(-1, 0, Nc[2]),
                cameraTarget=(Nc[0], Nc[1], Nc[2]))
            '''
            del mdb.models[model_name].steps['Step_G23']
       

    if calculate_E11 or calculate_E22 or calculate_E33:
        if calculate_E11:
            print('--------------------------------------------------------')
            print('Elastic properties in X-direction:')
            print('Exx   : ' + str(E11))
            print('nu_xy : ' + str(nu12))
            print('nu_xz : ' + str(nu13))
        if calculate_E22:
            print('--------------------------------------------------------')
            print('Elastic properties in Y-direction:')
            print('Eyy   : ' + str(E22))
            print('nu_yx : ' + str(nu21))
            print('nu_yz : ' + str(nu23))
        if calculate_E33:
            print('--------------------------------------------------------')
            print('Elastic properties in Z-direction:')
            print('Ezz   : ' + str(E33))
            print('nu_zx : ' + str(nu31))
            print('nu_zy : ' + str(nu32))
    if calculate_G12 or calculate_G13 or calculate_G23:
        if calculate_G12:
            print('--------------------------------------------------------')
            print('Shear modulus in XY-direction:')
            print('Gxy   : ' + str(G_12))
        if calculate_G13:
            print('--------------------------------------------------------')
            print('Shear modulus in XZ-direction:')
            print('Gxz   : ' + str(G_13))
        if calculate_G23:
            print('--------------------------------------------------------')
            print('Shear modulus in YZ-direction:')
            print('Gyz   : ' + str(G_23))
        print('--------------------------------------------------------')
    else:
        print('--------------------------------------------------------')
