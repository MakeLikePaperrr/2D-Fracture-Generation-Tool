import os
import numpy as np
from methods.generate_mesh_file import create_geo_file


NR_REAL = 10
DIR_INPUT = 'ensemble_2'
BASE_FILENAME_INPUT = lambda ith_real: f'real_{ith_real}.txt'
BASE_FILENAME_OUTPUT = lambda ith_real: f'{DIR_INPUT}_mesh_{ith_real}'
DIR_OUTPUT = f'{DIR_INPUT}_meshes'

char_len = 20
char_len_bound_mult = 1

margin = 25
x_min = 0 - margin
y_min = 0 - margin
x_max = 1000 + margin
y_max = 1000 + margin
box_data = np.array([[x_min, y_min], [x_max, y_min], [x_max, y_max], [x_min, y_max]])

if not os.path.exists(DIR_OUTPUT):
    os.makedirs(DIR_OUTPUT)

for i in range(1, NR_REAL + 1):
    act_frac_sys = np.genfromtxt(os.path.join(DIR_INPUT, BASE_FILENAME_INPUT(i)))
    filename = os.path.join(DIR_OUTPUT, BASE_FILENAME_OUTPUT(i))
    create_geo_file(act_frac_sys=act_frac_sys, filename=filename+'.geo', decimals=7,
                    height_res=25, char_len=char_len, box_data=box_data,
                    char_len_boundary=char_len*char_len_bound_mult)
    """
        NOTE: In gmsh you need to have under Options -> Geometry -> General -> uncheck "Remove duplicate ..." 
        otherwise meshing will crash/take too long!
    """
    os.system(f"gmsh {filename+'.geo'} -o {filename+'.msh'} -save")
