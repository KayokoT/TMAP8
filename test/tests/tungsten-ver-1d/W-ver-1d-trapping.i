# Tungsten adaption of verification Problem #1d from TMAP4/TMAP7 V&V document
# Permeation Problem with Trapping in Trapping-limited Case
# No Soret effect, or solubility included.

# Modeling parameters
node_num = 2000
thickness = '${units 2.5e-5 m -> mm}'
end_time = '${units 3 s}' #min
temperature = '${units 1500 K}'
diffusivity = '${units 1.632 mm^2/s}' #mm^2/min

# Trapping parameters
cl = ${units 6.11497e24 at/m^3 -> at/mm^3} #mobile hydrogen concentration
N = ${units 6.11497e28 at/m^3 -> at/mm^3} #trap site density
trapping_prefactor = '${units 6e16 1/s}' #1/min
release_prefactor = '${units 6e14 1/s}' #1/min
epsilon=${units 13925.421746094702 K} #trapping energy, 0.5 eV
#0.85 eV = 9863.840403483748 K
#1.2 eV = 13925.421746094702 K
#1.5 eV = 17406.777182618378 K
trapping_fraction = 0.01 # -
trap_per_free = 1e6

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
  [mobile]
  []
  [trapped]
  []
[]

[AuxVariables]
  [empty_sites]
  []
  [scaled_empty_sites]
  []
  [trapped_sites]
  []
  [total_sites]
  []
[]

[AuxKernels]
  [empty_sites]
    variable = empty_sites
    type = EmptySitesAux
    N = '${fparse N / cl}'
    Ct0 = ${trapping_fraction}
    trap_per_free = ${trap_per_free}
    trapped_concentration_variables = trapped
  []
  [scaled_empty]
    variable = scaled_empty_sites
    type = NormalizationAux
    normal_factor = ${cl}
    source_variable = empty_sites
  []
  [trapped_sites]
    variable = trapped_sites
    type = NormalizationAux
    normal_factor = ${trap_per_free}
    source_variable = trapped
  []
  [total_sites]
    variable = total_sites
    type = ParsedAux
    expression = 'trapped_sites + empty_sites'
    coupled_variables = 'trapped_sites empty_sites'
  []
[]

[Kernels]
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
  [coupled_time]
    type = ScaledCoupledTimeDerivative
    variable = mobile
    v = trapped
    factor = ${trap_per_free}
    extra_vector_tags = ref
  []
[]

[NodalKernels]
  [time]
    type = TimeDerivativeNodalKernel
    variable = trapped
  []
  [trapping]
    type = TrappingNodalKernel
    variable = trapped
    alpha_t = ${trapping_prefactor}
    N = '${fparse N / cl}'
    Ct0 = ${trapping_fraction}
    mobile_concentration = 'mobile'
    temperature = ${temperature}
    trap_per_free = ${trap_per_free}
    extra_vector_tags = ref
  []
  [release]
    type = ReleasingNodalKernel
    alpha_r = ${release_prefactor}
    temperature = ${temperature}
    detrapping_energy = ${epsilon}
    variable = trapped
  []
[]

[BCs]
  [left]
    type = FunctionDirichletBC
    variable = mobile
    function = 'BC_func'
    boundary = left
  []
  [right]
    type = DirichletBC
    variable = mobile
    value = 0
    boundary = right
  []
[]
[Functions]
  [BC_func]
    type = ParsedFunction
    expression = '${fparse cl / cl}*tanh(1e10 * t)'
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
  [min_trapped]
    type = NodalExtremeValue
    value_type = MIN
    variable = trapped
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
  dtmax = 1e-3
  solve_type = NEWTON
  scheme = BDF2
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = 'none'
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-7
    optimal_iterations = 9
    growth_factor = 1.1
    cutback_factor = 0.909
  []
[]

[Outputs]
  exodus = false
  csv = true
  [dof]
    type = DOFMap
    execute_on = initial
  []
  perf_graph = false
[]
