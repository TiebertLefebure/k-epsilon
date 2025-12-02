
# File paths for mesh and boundary data
mesh_files = {
    'MESH_DIRECTORY': 'Meshes/BackStep/Fine/mesh.xdmf',
    'FACET_DIRECTORY': 'Meshes/BackStep/Fine/facet.xdmf'
}

# Specify type of boundaries
boundary_markers = {
    'INFLOW': [4],
    'OUTFLOW': [2],
    'WALLS': [1, 3],
    'SYMMETRY': [5]
}

# Initial conditions
initial_conditions = {
    'U': (0.0, 0.0),
    'P': 0.0,
    'K': 1.73,
    'E': 1.46,
    'NU_TILDE': 1e-6
}

# Boundary conditions
boundary_conditions = {
    'INFLOW':{
        'U': (25.0, 0.0),
        'P': None,
        'K': 1.73,
        'E': 1.46,
        'NU_TILDE': None
    },
    'OUTFLOW':{
        'U': None,
        'P': 0.0,
        'K': None,
        'E': None,
        'NU_TILDE': None
    },
    'WALLS':{
        'U': (0.0, 0.0),
        'P': None,
        'K': 0.0,
        'E': None,
        'NU_TILDE': 0.0
    },
    'SYMMETRY':{
        'U': 0.0,
        'P': None,
        'K': 1.73,
        'E': 1.46,
        'NU_TILDE': None
    }
}

# Physical quantities
physical_prm = {
    'VISCOSITY': 0.000181818, # kinematic viscosity [m^2/s]
    'FORCE': (0.0, 0.0)
}

# Simulation parameters
simulation_prm = {
    'QUADRATURE_DEGREE': 2,
    'MAX_ITERATIONS': 3000, # initially: 3000
    'TOLERANCE': 1e-3, # initially: 1e-6
    'PICARD_RELAXATION': 0.1 # initially: 0.1
}

# Specify where results are saved for k-epsilon runs
saving_directory = {
    'PVD_FILES': 'Results/BackStep/PVD files/',
    'H5_FILES':  'Results/BackStep/H5 files/',
    'RESIDUALS': 'Results/BackStep/Residual files/'
}

# Dedicated saving directories for SA runs
saving_directory_SA = {
    'PVD_FILES': 'Results/BackStep_SA/PVD files/',
    'H5_FILES':  'Results/BackStep_SA/H5 files/',
    'RESIDUALS': 'Results/BackStep_SA/Residual files/'
}

# Summary result file for SA runs
results_file_SA = 'Results/BackStep_SA/results.txt'

# Specify what to do after simulation
post_processing = {
    'PLOT': True,
    'SAVE': True,
}
