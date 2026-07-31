# Tungsten adaption of verification Problem #1d from TMAP4/TMAP7 V&V document
# Permeation Problem with Trapping in Trapping-limited Case
# No Soret effect, or solubility included.

# Modeling parameters
node_num = 2000
end_time = 1 #s
trap_per_free = 1e3
# Trapping parameters

!include copper.params

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = ${node_num}
  xmax = ${Cu_thickness}
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
    N = '${fparse Cu_lattice_density / cl}'
    Ct0 = ${Cu_trapping_fraction}
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
    type = ADMatDiffusion
    variable = mobile
    diffusivity = ${Cu_D}
    extra_vector_tags = ref
  []
  [time]
    type = ADTimeDerivative
    variable = mobile
    extra_vector_tags = ref
  []
  [coupled_time]
    type = ADScaledCoupledTimeDerivative
    variable = mobile
    v = trapped
    mat_prop = ${trap_per_free}
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
    alpha_t = ${Cu_trapping_prefactor}
    N = '${fparse Cu_lattice_density / cl}'
    Ct0 = ${Cu_trapping_fraction}
    mobile_concentration = 'mobile'
    temperature = ${T}
    trap_per_free = ${trap_per_free}
    extra_vector_tags = ref
  []
  [release]
    type = ReleasingNodalKernel
    alpha_r = ${Cu_release_prefactor}
    temperature = ${T}
    detrapping_energy = ${Cu_epsilon}
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
    diffusivity = ${Cu_D}
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
  solve_type = NEWTON
  scheme = BDF2
  # petsc_options = '-pc_svd_monitor'
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'svd'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = 'none'
  automatic_scaling = true
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
  perf_graph = false
[]
