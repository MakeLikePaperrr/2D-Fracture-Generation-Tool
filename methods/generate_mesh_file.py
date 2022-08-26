import numpy as np


def create_geo_file(act_frac_sys, filename, decimals,
                    height_res, char_len, box_data, char_len_boundary, export_frac=True):
    """
    Creates geo file which serves as input to gmsh
    :param act_frac_sys: list of fractures in domain in format [[x1, y1, x2, y2], [...], ..., [...]]
    :param filename: name of the resulting geo-file
    :param decimals: data is rounded off to this number of decimals
    :param height_res: height of the resulting 1-layer 3D reservoir
    :param char_len: characteristic length of the resulting mesh
    :param box_data: coordinates of the box-data around the fracture network
    :param char_len_boundary: characteristic length of mesh elements at the boundary
    :param export_frac: boolean which exports fractures into the meshed file
    :return:
    """
    act_frac_sys = np.round(act_frac_sys * 10 ** decimals) * 10 ** (-decimals)
    num_segm_tot = act_frac_sys.shape[0]
    unique_nodes = np.unique(np.vstack((act_frac_sys[:, :2], act_frac_sys[:, 2:])), axis=0)
    num_nodes_tot = unique_nodes.shape[0]
    f = open(filename, "w+")

    # Note: always end statement in GMSH with ";"
    # Note: comments in GMSH are as in C(++) "//"
    f.write('// Geo file which meshes the input mesh from act_frac_sys.\n')
    f.write('// Change mesh-elements size by varying "lc" below.\n\n')

    # Can specify the type of meshing algorithm for 2D meshing here:
    # f.write('// Specify meshing algorithm:\n')
    # f.write('-algo meshadapt;\n\n')

    # Set some parameters in the model:
    f.write('lc = {:1.3f};\n'.format(char_len))
    f.write('lc_box = {:1.3f};\n'.format(char_len_boundary))
    f.write('height_res = {:4.3f};\n\n'.format(height_res))

    # Allocate memory for points_created array and counters:
    points_created = np.zeros((num_nodes_tot,), dtype=bool)
    line_count = 0
    point_count = 0

    for ii in act_frac_sys:
        # Take two points per segment and write to .geo file:
        # e.g. point: Point(1) = {.1, 0, 0, lc};
        nodes = np.zeros((2,), dtype=int)
        nodes[0] = np.where(np.logical_and(ii[0] == unique_nodes[:, 0], ii[1] == unique_nodes[:, 1]))[0]
        nodes[1] = np.where(np.logical_and(ii[2] == unique_nodes[:, 0], ii[3] == unique_nodes[:, 1]))[0]

        # Check if first point is already created, if not, add it:
        if not points_created[nodes[0]]:
            points_created[nodes[0]] = True
            point_count += 1
            f.write('Point({:d}) = {{{:8.5f}, {:8.5f}, 0, lc}};\n'.format(nodes[0] + 1,
                                                                           unique_nodes[nodes[0], 0],
                                                                           unique_nodes[nodes[0], 1]))

        if not points_created[nodes[1]]:
            points_created[nodes[1]] = True
            point_count += 1
            f.write('Point({:d}) = {{{:8.5f}, {:8.5f}, 0, lc}};\n'.format(nodes[1] + 1,
                                                                           unique_nodes[nodes[1], 0],
                                                                           unique_nodes[nodes[1], 1]))

        line_count += 1
        f.write('Line({:d}) = {{{:d}, {:d}}};\n\n'.format(line_count, nodes[0] + 1, nodes[1] + 1))

    # Store some internal variables for gmsh (used later after extrude):
    f.write('num_points_frac = newp - 1;\n')
    f.write('num_lines_frac = newl - 1;\n\n')

    # Write the box_data (box around fracture network in which we embed the fractures)
    f.write('// Extra points for boundary of domain:\n')
    for ii in range(4):
        # For every corner of the box:
        point_count += 1
        f.write('Point({:d}) = {{{:8.5f}, {:8.5f}, 0, lc_box}};\n'.format(point_count,
                                                                          box_data[ii, 0], box_data[ii, 1]))

    # Add four lines for each side of the box:
    f.write('\n// Extra lines for boundary of domain:\n')
    line_count += 1
    f.write('Line({:d}) = {{{:d}, {:d}}};\n'.format(line_count, point_count - 3, point_count - 2))

    line_count += 1
    f.write('Line({:d}) = {{{:d}, {:d}}};\n'.format(line_count, point_count - 2, point_count - 1))

    line_count += 1
    f.write('Line({:d}) = {{{:d}, {:d}}};\n'.format(line_count, point_count - 1, point_count - 0))

    line_count += 1
    f.write('Line({:d}) = {{{:d}, {:d}}};\n'.format(line_count, point_count - 0, point_count - 3))

    # Make Curve loop for the boundary:
    f.write('\n// Create line loop for boundary surface:\n')
    f.write('Curve Loop(1) = {{{:d}, {:d}, {:d}, {:d}}};\n'.format(line_count - 3, line_count - 2,
                                                                   line_count - 1, line_count))
    f.write('Plane Surface(1) = {1};\n\n')
    f.write('Curve{1:num_lines_frac} In Surface{1};\n')

    # Extrude model to pseuo-3D:
    f.write('\n// Extrude surface with embedded features\n')
    f.write('Extrude {0, 0, height_res}{ Surface {1}; Layers{1}; Recombine;}\n')
    f.write('Physical Volume("matrix", 9991) = {1};\n')

    f.write('num_surfaces_before = news;\n')
    # f.write('Extrude {0, 0, height_res}{ Line {1:num_lines_frac}; Layers{1}; Recombine;}\n')
    f.write('num_surfaces_after = news - 1;\n')
    f.write('num_surfaces_fracs = num_surfaces_after - num_surfaces_before;\n\n')
    for ii in range(act_frac_sys.shape[0]):
        # f.write('Physical Surface({:d}) = {{num_surfaces_before + {:d}}};\n'.format(90000 + ii, ii))
        f.write('Extrude {{0, 0, height_res}}{{ Line {{{:d}}}; Layers{{1}}; Recombine;}}\n'.format(ii + 1))
        if export_frac:
            f.write('Physical Surface({:d}) = {{news - 1}};\n'.format(90000 + ii))

    # Create mesh and perform coherency check:
    f.write('Mesh 3;  // Generalte 3D mesh\n')
    f.write('Coherence Mesh;  // Remove duplicate entities\n')
    f.write('Mesh.MshFileVersion = 2.1;\n')
    f.close()
    return 0