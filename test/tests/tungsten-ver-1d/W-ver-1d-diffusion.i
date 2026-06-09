# Verification Problem #1d from TMAP4/TMAP7 V&V document
# Permeation Problem with Trapping in Diffusion-limited Case
# No Soret effect, or solubility included.

# Modeling parameters
node_num = 2000
interval_time_min = 1e-15 #min
interval_time = 1e-4 #min
thickness = ${units 2.5e-2 mm}
end_time = 0.014 #min
temperature = '${units 1500 K}'
diffusivity = 1.632 #mm^2/min 2.72e-8 m^2/s

# Trapping parameters
density = '${units 6.11497e28 at/m^3 -> at/mm^3}' #lattice density
cl = ${units 6.11497e24 at/m^3 -> at/mm^3} #lattice density * 1e-4
trapping_prefactor = 6e16 #1/min
release_prefactor = 6e14 #1/min
release_energy = '${units 5802.259 K}' #trap energy / k
trapping_fraction = 0.01 # -

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = ${node_num}
  xmax = ${thickness}
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Variables]
  # mobile tritium variable
  [mobile]
  []
  # trapped tritium variable
  [trapped]
  []
[]

[Kernels]
  # kernel for mobile tritium
  [diff]
    type = MatDiffusion
    variable = mobile
    diffusivity = ${diffusivity}
    extra_vector_tags = ref
  []
  [time]
    type = TimeDerivative
    variable = mobile
    extra_vector_tags = ref
  []
  # kernel for trapped tritium
  [coupled_time]
    type = CoupledTimeDerivative
    variable = mobile
    v = trapped
    extra_vector_tags = ref
  []
[]

[NodalKernels]
  # kernel for trapped tritium in nodal sites
  [time]
    type = TimeDerivativeNodalKernel
    variable = trapped
  []
  [trapping]
    type = TrappingNodalKernel
    variable = trapped
    alpha_t = ${trapping_prefactor}
    N = '${fparse density / cl}'
    Ct0 = ${trapping_fraction}
    mobile_concentration = 'mobile'
    temperature = ${temperature}
    extra_vector_tags = ref
  []
  [release]
    type = ReleasingNodalKernel
    alpha_r = ${release_prefactor}
    temperature = ${temperature}
    detrapping_energy = ${release_energy}
    variable = trapped
  []
[]

[Functions]
  [BC_func]
    type = ParsedFunction
    expression = '${fparse cl / cl}*tanh(1000 * t )'
  []
[]

[BCs]
  [left]
    type = FunctionDirichletBC
    variable = mobile
    function = 'BC_func'
    boundary = left
  []
  # [left]
  #   type = DirichletBC
  #   variable = mobile
  #   value = 1
  #   boundary = left
  # []
  [right]
    type = DirichletBC
    variable = mobile
    value = 0
    boundary = right
  []
[]

[Postprocessors]
  [outflux]
    type = SideDiffusiveFluxAverage
    boundary = 'right'
    diffusivity = ${diffusivity}
    variable = mobile
  []
  [scaled_outflux]
    type = ScalePostprocessor
    value = outflux
    scaling_factor = ${cl}
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
  end_time = ${end_time}
  dt = ${interval_time}
  dtmin = ${interval_time_min}
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  verbose = true
  compute_scaling_once = false
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = ${interval_time}
    growth_factor = 1.5
    cutback_factor = 0.90
    optimal_iterations = 9
  []
[]
  

[Outputs]
  csv = true
  exodus = false
  [dof]
    type = DOFMap
    execute_on = initial
  []
  perf_graph = false
[]