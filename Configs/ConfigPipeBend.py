
# File paths for mesh and boundary data
mesh_files = {
    'MESH_DIRECTORY': 'Meshes/PipeBend/mesh.xdmf',
    'FACET_DIRECTORY': 'Meshes/PipeBend/facet.xdmf'
}

# Specify type of boundaries
boundary_markers = {
    'INFLOW': [2],
    'OUTFLOW': [3],
    'WALLS': [4],
    'SYMMETRY': None
}

# 2D U-bend (Ansys Manual VMFL048)
# radius = 14 mm
# radius of curvature = 125 mm
# straight section lenghts = 10*D = 280 mm


# inlet bulk velocity U = 4.0 m/s

# choose turbulence intensity I = 5%

# set turbulent length scale l ≈ 0.07*D = 2.0 mm

# initial & inflow boundary conditions for K & E:
# K = 1.5 * (U * I)^2 = 0.06  [m^2/s^2]
# E = Cµ^(3/4) * K^(3/2) / l ≈ 1.2  [m^2/s^3], using Cµ = 0.09



# Initial conditions
initial_conditions = {
    'U': (0.5, 0.0), #(20.0, 0.0)
    'P': 2.0,
    'K': 0.06,  #1.5
    'E': 1.2    #2.23
}

# Boundary conditions
boundary_conditions = {
    'INFLOW':{
        'U': (4.0, 0.0), #0.0
        'P': None, #2.0
        'K': 0.06, #None
        'E': 1.2   #None
    },
    'OUTFLOW':{
        'U': None,
        'P': 0.0,
        'K': None,
        'E': None
    },
    'WALLS':{
        'U': (0.0, 0.0),
        'P': None,
        'K': 0.0,
        'E': None
    },
    'SYMMETRY':{
        'U': None,
        'P': None,
        'K': None,
        'E': None
    }
}

# Physical quantities
physical_prm = {
    'VISCOSITY': 1.5e-5,  # air at room temperature   #0.00181818
    'FORCE': (0.0, 0.0),
    'STEP_SIZE': 5e-4 #0.005
}

# Reynolds number:
# Re = U * D / ν = 4 * 0.028 / 1.5e-5 ≈ 7.5 × 10^3

# Simulation parameters
simulation_prm = {
    'QUADRATURE_DEGREE': 2,
    'MAX_ITERATIONS': 3000,
    'TOLERANCE': 1e-5,
    'CFL_RELAXATION': 0.1 #0.25
}

# Specify where results are saved
saving_directory = {
    'PVD_FILES': 'Results/PipeBend/PVD files/',
    'H5_FILES':  'Results/PipeBend/H5 files/',
    'RESIDUALS': 'Results/PipeBend/Residual files/'
}

# Specify what to do after simulation
post_processing = {
    'PLOT': True,
    'SAVE': True,
}


#docker run -it --rm \
#    -v /Users/tiebertlefebure/Documents/FEniCSx/Turbulence_models:/home/fenics/shared \
#    --platform linux/amd64 \
#    quay.io/fenicsproject/stable:current


#cd ~/shared

#ls 

#python3 PipeBendSimulation.py