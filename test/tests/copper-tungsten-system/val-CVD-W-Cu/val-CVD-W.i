# Permeation problem in a single layer (W) to compare against CVD W-Cu experimental study (Sun 2025)

!include val-CVD-W-Cu.params

[Mesh]
  coord_type = 'XYZ'
  [generated]
    type = CartesianMeshGenerator
    dim = 1
    dx = '${W_only_thickness}'
    ix = '400'  # number of elements per subdomain
    subdomain_id = '0'
  []
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Materials]
  [diffusivity]
    type = ADGenericConstantMaterial
    prop_names = 'diffusivity_W' 
    prop_values = '${W_D}' 
  []
  [flux_on_left]
    type = ADParsedMaterial
    block = '0'
    coupled_variables = C_M_W
    property_name = 'flux_on_left'
    expression = '(${W_J0_coef} * ${J0}) - (${W_kr} * C_M_W ^ 2)' #at/mum^2/s
  []
[]

[Variables]
  [C_M_W] #mobile concentration of H in W, at/mum^3
    block = '0'
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
  # recombination-limited upstream
  [left_recombination]
    type = ADMatNeumannBC
    variable = C_M_W
    boundary = 'left'
    value = 1
    boundary_material = flux_on_left
    extra_vector_tags = ref
  []
  # diffusion-limited downstream
  [right]
    type = ADDirichletBC
    variable = C_M_W
    boundary = 'right'
    value = 0
  []
[]

[Postprocessors]
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
    time_integration_scheme = trapezoidal-rule
    outputs = main
    execute_on = "INITIAL TIMESTEP_END"
  []

  [mass_from_concentration]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_W
    block = '0'
    execute_on = "INITIAL TIMESTEP_END"
    outputs = main
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
  scheme = bdf2
  #end_time = ${simulation_time}
  solve_type = NEWTON
  line_search = 'bt'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  dtmax = 5
  steady_state_tolerance = 1e-6
  steady_state_detection = true
  steady_state_start_time = 1.0
  automatic_scaling = true
  compute_scaling_once = true
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-8
    optimal_iterations = 6
    growth_factor = 1.15
    cutback_factor = 0.5
  []
[]

[Outputs]
  file_base = val_W_${T}
  [exodus]
    type = Exodus
  []
  [main]
    type = CSV
    file_base = val_W_${T}_main
  []
[]
