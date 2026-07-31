# Permeation problem in a single layer (Cu)
# No trapping
# Adapted from val 1dd
# BCs: Sieverts (left) and recombination (right)

!include WCu-atmum3.params
Cu_probe = ${units 1 mum} #probe in Cu layer, measured from left side of Cu

[Mesh]
  coord_type = 'XYZ'
  [generated]
    type = CartesianMeshGenerator
    dim = 1
    dx = '${Cu_thickness}'
    ix = '400'
    subdomain_id = '0' 
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
    prop_names = 'diffusivity_Cu' 
    prop_values = '${Cu_D}' 
  []

  [flux_recombination_surface]
    type = ADParsedMaterial
    block = '0'
    coupled_variables = C_M_Cu
    property_name = 'flux_recombination_surface'
    expression = '- 2 * ${Cu_kr} * C_M_Cu ^ 2' #H2 molecules/mum^2/s
  []
[]

[Variables]
  [C_M_Cu] #mobile concentration of H in Cu, at/mum^3
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
  [Cu_diff]
    type = ADMatDiffusion
    block = '0'
    variable = C_M_Cu
    diffusivity = diffusivity_Cu
    extra_vector_tags = ref
  []
  [Cu_time]
    type = ADTimeDerivative
    block = '0'
    variable = C_M_Cu
    extra_vector_tags = ref
  []
[]

[BCs]
  [left_flux]
    type = EquilibriumBC
    Ko = ${Cu_S} 
    activation_energy = 0
    boundary = 'left'
    enclosure_var = ${pressure} #H2 domain
    temperature = ${T}
    variable = C_M_Cu
    p = 0.5 #Sievert's Law
  []

  [right_recombination]
    type = ADMatNeumannBC
    variable = C_M_Cu
    boundary = 'right'
    value = 1
    boundary_material = flux_recombination_surface
    extra_vector_tags = ref
  []

  # zero concentration on right boundary - use to check diffusion front
  # [right]
  #   type = ADDirichletBC
  #   variable = C_M_Cu
  #   boundary = right
  #   value = 0
  # []
[]

[VectorPostprocessors]
  [Cu]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '${Cu_thickness} 0 0'
    num_points = 200
    sort_by = 'x'
    variable = C_M_Cu
    outputs = vector
  []
[]

[Postprocessors]
  # Used to obtain varying concentration with time at a point
  [concentration_at_x_Cu]
    type = PointValue
    variable = C_M_Cu
    point = '${Cu_probe} 0 0'
  []

  # Diffusion front verification
  [exact_diffusion_length]   # Analytical Diffusion length (time-independent BC required for this to be correct)
    type = ParsedPostprocessor
    expression = 'sqrt(pi*D*t)'
    constant_names = 'D pi'
    constant_expressions = '${Cu_D} 3.1415926535897932'
    use_t = true
    outputs = main
  []

  [gradient_left_boundary]
    type = ADSideDiffusiveFluxAverage
    boundary = 'left'
    variable = C_M_Cu
    diffusivity = 1
    outputs = none
  []

  [left_concentration] 
    type = SideAverageValue
    boundary = 'left'
    variable = C_M_Cu
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
    variable = C_M_Cu
    diffusivity = ${Cu_D}
    outputs = none
  []

  [right_outflux]
    type = ADSideDiffusiveFluxIntegral
    boundary = 'right'
    variable = C_M_Cu
    diffusivity = ${Cu_D}
    outputs = main
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
    outputs = main
  []

  [mass_from_concentration]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_Cu
    block = '0'
    outputs = main
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  line_search = 'bt'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  dtmax = 1e-2
  steady_state_tolerance = 1e-6
  steady_state_detection = true
  steady_state_start_time = 0.01
  automatic_scaling = true
  compute_scaling_once = true
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-8
    optimal_iterations = 6
    growth_factor = 1.15
    cutback_factor = 0.9
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
    execute_on = 'FINAL'
  []
[]
