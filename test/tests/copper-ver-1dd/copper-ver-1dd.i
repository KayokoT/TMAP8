# Verification Problem #1dd from TMAP7 V&V document
# Permeation Problem without Trapping sites
# No Soret effect or solubility included.

# modeling parameters
nx_num = 200 # (-)

simulation_time = ${units 1 s}
!include copper.params

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]
 
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = ${nx_num}
  xmax = '${Cu_thickness}'
[]

[Variables]
  [mobile]
  []
[]

[Kernels]
  [diff]
    type = MatDiffusion
    variable = mobile
    diffusivity = ${Cu_D}
  []
  [time]
    type = TimeDerivative
    variable = mobile
  []
[]

[BCs]
  [left]
    type = ADFunctionDirichletBC
    variable = mobile
    function = 'BC_func'
    boundary = left
    extra_vector_tags = ref
  []
  [right]
    type = ADDirichletBC
    variable = mobile
    value = 0
    boundary = right
    extra_vector_tags = ref
  []
[]

[Functions]
  [BC_func]
    type = ParsedFunction
    expression = '${fparse ${cl} * tanh(1e10 * t)}' #ramps up to 1
  []
[]

[Postprocessors]
  [outflux]
    type = SideDiffusiveFluxAverage
    boundary = 'right'
    diffusivity = ${Cu_D}
    variable = mobile
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  end_time = ${simulation_time}
  solve_type = NEWTON
  scheme = BDF2
  nl_abs_tol = 1e-13
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  verbose = true
  compute_scaling_once = false
    [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-5
    optimal_iterations = 6
    growth_factor = 1.05
    cutback_factor = 0.90
  []
[]

[Outputs]
  exodus = true
  csv = true
  [dof]
    type = DOFMap
    execute_on = initial
  []
  perf_graph = true
[]
