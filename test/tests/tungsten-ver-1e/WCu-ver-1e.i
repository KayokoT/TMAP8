# Verification Problem #1e from TMAP4/TMAP7 V&V document
# Permeation problem in a composite layer (Cu & W)
# No Soret effect, trapping, or solubility included.

# Numerical parameters
nx_num = 2000 # -
simulation_time = ${units 0.2 s}

T_Cu = ${units 50 mum -> m} #Cu layer thickness
T_W = ${units 25 mum -> m} #W layer thickness
D_ver = ${units 24 mum -> m} #probe in W layer
D_ver_Cu = ${units 49 mum -> m} #probe in Cu layer
Diffusivity_Cu = ${units 5.9968958217e-8 m^2/s} #at 1500K
Diffusivity_W = ${units 2.72e-8 m^2/s} #at 1500K
length_Cu = ${units 50 mum -> m}
initial_concentration = ${units 50.7079 mol/m^3} #default value 50.7079

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = ${nx_num}
  xmax = ${fparse ${T_Cu} + ${T_W} }
  allow_renumbering = false
[]

[Variables]
  [u]
  []
[]

[Functions]
  # Diffusivity assign based on different material domains
  [diffusivity_value]
    type = ParsedFunction
    expression = 'if(x < ${length_Cu}, ${Diffusivity_Cu}, ${Diffusivity_W} )'
  []
[]

[Kernels]
  [diff]
    type = FunctionDiffusion
    variable = u
    function = diffusivity_value
  []
  [time]
    type = TimeDerivative
    variable = u
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = ${initial_concentration}
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 0
  []
[]

# Used while obtaining steady-state solution
[VectorPostprocessors]
  [line]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '${Mesh/xmax} 0 0'
    num_points = ${Mesh/nx}
    sort_by = 'x'
    variable = u
    outputs = vector_postproc
  []
[]

[Postprocessors]
  # Used to obtain varying concentration with time at a
  # point in W layer 'x' um from Cu/W boundary
  [concentration_at_x_W]
    type = PointValue
    variable = u
    point = '${fparse ${T_Cu} + ${D_ver}} 0 0'
    outputs = 'csv'
  []
  [concentration_at_x_Cu]
    type = PointValue
    variable = u
    point = '${D_ver_Cu} 0 0'
    outputs = 'csv'
  []
[]

[Executioner]
  type = Transient
  end_time = ${simulation_time}
  dtmax = 0.1
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  scheme = 'bdf2'
  nl_rel_tol = 1e-50 # Make this really tight so that our absolute tolerance criterion is the one
  # we must meet
  nl_abs_tol = 1e-12
  abort_on_solve_fail = true
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-5
    optimal_iterations = 4
    growth_factor = 1.05
    cutback_factor = 0.9
  []
[]

[Outputs]
  [exodus]
    type = Exodus
  []
  [csv]
    type = CSV
  []
  [vector_postproc]
    type = CSV
    sync_times = ${simulation_time}
    sync_only = true
  []
[]
