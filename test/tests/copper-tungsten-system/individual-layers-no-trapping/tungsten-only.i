# Permeation problem in a single layer (W)
# No trapping
# Adapted from val 1dd
# BCs: Sieverts' (left) and infinite sink (right)
# Numerical parameters
simulation_time = ${units 0.015 s}
!include WCu-atmum3.params
W_probe = ${units 1 mum} #probe in W layer

[Mesh]
  coord_type = 'XYZ' # 3D Cartesian coordinates
  [generated]
    type = CartesianMeshGenerator
    dim = 1
    dx = '${W_thickness}'
    ix = '400' # 400'   # number of elements per subdomain
    subdomain_id = '0' # 2' 
  []
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Materials]
  [diffusivity_solubility]
    type = ADGenericConstantMaterial
    prop_names = 'diffusivity_W' 
    prop_values = '${W_D}' 
  []
[]

[Variables]
  [C_M_W] #mobile concentration of H in W, at/mum^3
    block = '0'
  []
[]

[AuxVariables]
  [pressure]
    family = SCALAR
    initial_condition = ${pressure}
  []
[]

[Kernels]
  [W_diff]
    type = ADMatDiffusion
    block = '0'
    variable = C_M_W
    diffusivity = diffusivity_W
    extra_vector_tags = ref
  []
  [W_time]
    type = ADTimeDerivative
    block = '0'
    variable = C_M_W
    extra_vector_tags = ref
  []
[]

[BCs]
  [left_flux]
    type = EquilibriumBC
    Ko = ${W_S} 
    activation_energy = 0
    boundary = 'left'
    enclosure_var = ${pressure} #H2 domain
    temperature = ${T}
    variable = C_M_W
    p = 0.5 #Sievert's Law
  []
  [right]
    type = ADDirichletBC
    variable = C_M_W
    value = 0
    boundary = 'right'
  []
[]

[VectorPostprocessors]
  [W]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '${W_thickness} 0 0'
    num_points = 200
    sort_by = 'x'
    variable = C_M_W
    outputs = vector
  []
[]

[Postprocessors]
  # Used to obtain varying concentration with time at a point
  [concentration_at_x_W]
    type = PointValue
    variable = C_M_W
    point = '${W_probe} 0 0'
  []

  # Diffusion front verification
  [exact_diffusion_length]   # Analytical Diffusion length (time-independent BC required for this to be correct)
    type = ParsedPostprocessor
    expression = 'sqrt(pi*D*t)'
    constant_names = 'D pi'
    constant_expressions = '${W_D} 3.1415926535897932'
    use_t = true
    outputs = main
  []

  [gradient_left_boundary]
    type = ADSideDiffusiveFluxAverage
    boundary = 'left'
    variable = C_M_W
    diffusivity = 1
    outputs = none
  []

  [left_concentration] 
    type = SideAverageValue
    boundary = 'left'
    variable = C_M_W
    outputs = none
  []

  [simulated_diffusion_length] # x-intercept of tangent line at interface
    type = ParsedPostprocessor
    expression = '-left_concentration/gradient_left_boundary'
    pp_names = 'left_concentration gradient_left_boundary'
    outputs = main
  []

  # conservation of mass check
  [left_influx] # Influx at left boundary
    type = ADSideDiffusiveFluxIntegral
    boundary = 'left'
    variable = C_M_W
    diffusivity = ${W_D}
    outputs = none
  []

  [right_outflux]
    type = ADSideDiffusiveFluxIntegral
    boundary = 'right'
    variable = C_M_W
    diffusivity = ${W_D}
    outputs = none
  []

  [flux_difference]
    type = ParsedPostprocessor
    expression = '-left_influx - right_outflux' # negative sign on influx to account for outward normal vector direction
    pp_names = 'left_influx right_outflux'
    outputs = none
  []
  # mass accumulated as inferred from boundary flux
  [mass_from_flux]
    type = TimeIntegratedPostprocessor
    value = flux_difference
    time_integration_scheme = IMPLICIT-EULER
    outputs = none
  []

  [mass_from_concentration]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_W
    block = '0'
    outputs = none
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  line_search = 'none'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-20
  dtmax = 1e-4
  steady_state_tolerance = 1e-6
  steady_state_detection = true
  steady_state_start_time = 0.01
  automatic_scaling = true
  compute_scaling_once = false
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-7
    optimal_iterations = 6
    growth_factor = 1.05
    cutback_factor = 0.95
  []
[]

[Outputs]
  [exodus]
    type = Exodus
  []
  [main]
    type = CSV
  []
  [vector]
    type = CSV
    sync_times = ${simulation_time}
    sync_only = true
  []
[]
