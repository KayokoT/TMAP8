# Permeation problem in a composite layer (W/Cu/W) adapted from Verification Problem #1e
# W faces the H source
# No trapping
# BCs: Sieverts (left) & recombination (right)

# Numerical parameters
!include WCu-atmum3.params
penalty = 1e6 #check conc&sol_ratio to ensure BOTH interface jumps are properly enforced
RUN_ID = "run01" #default, change this in terminal by appending RUN_ID=run**
W1_probe = ${units 1 mum} #probe in W1 layer, measured from left boundary of W1
W2_probe = ${units 1 mum} #probe in W2 layer, measured from left boundary of W2
Cu_probe = ${units 1 mum} #probe in Cu layer, measured from left boundary of Cu
[Mesh]
  coord_type = 'XYZ'
  [generated]
    type = CartesianMeshGenerator
    dim = 1
    dx = '${W_thickness} ${Cu_thickness} ${W_thickness}'
    ix = '800 400 800'   # number of elements per subdomain
    subdomain_id = '0 1 2'   # W1 / Cu / W2
  []
  [interface_W1_to_Cu]
    type = SideSetsBetweenSubdomainsGenerator
    input = generated
    primary_block = '0'   # W1
    paired_block = '1'    # Cu
    new_boundary = 'interface_W1_to_Cu'
  []
  [interface_Cu_to_W1]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface_W1_to_Cu
    primary_block = '1'   # Cu
    paired_block = '0'    # W1
    new_boundary = 'interface_Cu_to_W1'
  []
  [interface_Cu_to_W2]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface_Cu_to_W1
    primary_block = '1'   # Cu
    paired_block = '2'    # W2
    new_boundary = 'interface_Cu_to_W2'
  []
  [interface_W2_to_Cu]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface_Cu_to_W2
    primary_block = '2'   # W2
    paired_block = '1'    # Cu
    new_boundary = 'interface_W2_to_Cu'
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
    prop_names = 'diffusivity_W diffusivity_Cu solubility_W solubility_Cu'
    prop_values = '${W_D} ${Cu_D} ${W_S} ${Cu_S}'
  []

  [converter_to_nonAD]
    type = MaterialADConverter
    ad_props_in = 'solubility_W solubility_Cu'
    reg_props_out = 'solubility_W_nonAD solubility_Cu_nonAD'
    outputs = 'exodus'
  []
  #first interface
  [interface_jump_W1_Cu]
    type = SolubilityRatioMaterial
    solubility_primary = solubility_W
    solubility_secondary = solubility_Cu
    boundary = interface_W1_to_Cu
    concentration_primary = C_M_W1
    concentration_secondary = C_M_Cu
  []

  #second interface
  [interface_jump_Cu_W2]
    type = SolubilityRatioMaterial
    solubility_primary = solubility_Cu
    solubility_secondary = solubility_W
    boundary = interface_Cu_to_W2
    concentration_primary = C_M_Cu
    concentration_secondary = C_M_W2
  []

  # downstream (right) boundary
  [flux_recombination_surface]
    type = ADParsedMaterial
    block = '2'
    coupled_variables = C_M_W2
    property_name = 'flux_recombination_surface'
    expression = '-2 * ${W_kr} * C_M_W2 ^ 2' #at/mum^2/s
  []
[]

[Variables]
  [C_M_W1] #mobile concentration of H in first W layer, at/mum^3
    block = '0'
  []
  [C_M_Cu] #mobile concentration of H in bulk Cu layer, at/mum^3
    block = '1'
  []
  [C_M_W2] #mobile concentration of H in second W layer, at/mum^3
    block = '2'
  []
[]

[AuxVariables]
  [pressure]
    family = SCALAR
    initial_condition = ${pressure}
  []
[]

[Kernels]
  [W1_diff]
    type = ADMatDiffusion
    block = '0'
    variable = C_M_W1
    diffusivity = diffusivity_W
    extra_vector_tags = ref
  []
  [W1_time]
    type = ADTimeDerivative
    block = '0'
    variable = C_M_W1
    extra_vector_tags = ref
  []
  [Cu_diff]
    type = ADMatDiffusion
    block = '1'
    variable = C_M_Cu
    diffusivity = diffusivity_Cu
    extra_vector_tags = ref
  []
  [Cu_time]
    type = ADTimeDerivative
    block = '1'
    variable = C_M_Cu
    extra_vector_tags = ref
  []
  [W2_diff]
    type = ADMatDiffusion
    block = '2'
    variable = C_M_W2
    diffusivity = diffusivity_W
    extra_vector_tags = ref
  []
  [W2_time]
    type = ADTimeDerivative
    block = '2'
    variable = C_M_W2
    extra_vector_tags = ref
  []
[]

[BCs]
  #Sieverts' equilbrium upstream
  [left_flux]
    type = EquilibriumBC
    Ko = ${W_S} 
    boundary = left
    activation_energy = 0
    enclosure_var = ${pressure} #H2 domain
    temperature = ${T}
    variable = C_M_W1
    p = 0.5 #Sievert's Law
  []
  # recombination downstream
  [right_recombination]
    type = ADMatNeumannBC
    variable = C_M_W2
    boundary = 'right'
    value = 1
    boundary_material = flux_recombination_surface
    extra_vector_tags = 'ref'
  []
[]

[InterfaceKernels]
  # Penalized continuity with solubility jump at W1 / Cu interface
  [W1_to_Cu_interface]
    type = ADPenaltyInterfaceDiffusion
    variable = C_M_W1
    neighbor_var = C_M_Cu
    penalty = ${penalty}
    jump_prop_name = solubility_ratio
    boundary = 'interface_W1_to_Cu'
    extra_vector_tags = 'ref'
  []
  # Penalized continuity with solubility jump at Cu / W2 interface
  [Cu_to_W2_interface]
    type = ADPenaltyInterfaceDiffusion
    variable = C_M_Cu
    neighbor_var = C_M_W2
    penalty = ${penalty}
    jump_prop_name = solubility_ratio
    boundary = 'interface_Cu_to_W2'
    extra_vector_tags = 'ref'
  []
[]

[VectorPostprocessors]
  [W1]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '${W_thickness} 0 0'
    num_points = 200
    sort_by = 'x'
    variable = C_M_W1
    outputs = vector
  []
  [Cu]
    type = LineValueSampler
    start_point = '${W_thickness} 0 0'
    end_point = '${fparse ${W_thickness} + ${Cu_thickness}} 0 0'
    num_points = 200
    sort_by = 'x'
    variable = C_M_Cu
    outputs = vector
  []
  [W2]
    type = LineValueSampler
    start_point = '${fparse ${W_thickness} + ${Cu_thickness}} 0 0'
    end_point = '${fparse (2 * ${W_thickness}) + ${Cu_thickness}} 0 0'
    num_points = 200
    sort_by = 'x'
    variable = C_M_W2
    outputs = vector
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  #concentration at points to plot over time
  # Used to obtain varying concentration with time at a point
  [concentration_at_x_W1]
    type = PointValue
    variable = C_M_W1
    point = '${W1_probe} 0 0'
  []

  [concentration_at_x_Cu]
    type = PointValue
    variable = C_M_Cu
    point = '${fparse ${W_thickness} + ${Cu_probe}} 0 0'
  []

  [concentration_at_x_W2]
    type = PointValue
    variable = C_M_W2
    point = '${fparse ${W_thickness} + ${Cu_thickness}+ ${W2_probe}} 0 0'
  []

  # check concentration ratio against solubility to ensure penalty is properly enforced
  # W1/Cu interface
  [W1_S]
    type = SideAverageMaterialProperty
    property = solubility_W_nonAD
    boundary = interface_W1_to_Cu
  []

  [Cu_S_at_W1_interface]
    type = SideAverageMaterialProperty
    property = solubility_Cu_nonAD
    boundary = interface_W1_to_Cu
  []

  [gold_solubility_ratio_1]
    type = ParsedPostprocessor
    pp_names = 'Cu_S_at_W1_interface W1_S'
    expression = 'Cu_S_at_W1_interface / W1_S'
  []

  [Cu_interface_1]
    type = SideAverageValue
    boundary = interface_Cu_to_W1
    variable = C_M_Cu
  []

  [W1_interface]
    type = SideAverageValue
    boundary = interface_W1_to_Cu
    variable = C_M_W1
  []

  [variable_ratio_1]
    type = ParsedPostprocessor
    pp_names = 'Cu_interface_1 W1_interface'
    expression = 'Cu_interface_1 / W1_interface'
  []
  # Cu/W2 interface
  [Cu_S_at_W2_interface]
    type = SideAverageMaterialProperty
    property = solubility_Cu_nonAD
    boundary = interface_Cu_to_W2
  []

  [W2_S]
    type = SideAverageMaterialProperty
    property = solubility_W_nonAD
    boundary = interface_Cu_to_W2
  []

  [gold_solubility_ratio_2]
    type = ParsedPostprocessor
    pp_names = 'W2_S Cu_S_at_W2_interface'
    expression = 'W2_S / Cu_S_at_W2_interface'
  []

  [W2_interface]
    type = SideAverageValue
    boundary = interface_W2_to_Cu
    variable = C_M_W2
  []

  [Cu_interface_2]
    type = SideAverageValue
    boundary = interface_Cu_to_W2
    variable = C_M_Cu
  []

  [variable_ratio_2]
    type = ParsedPostprocessor
    pp_names = 'W2_interface Cu_interface_2'
    expression = 'W2_interface / Cu_interface_2'
  []
  # conservation of mass check
  # inventory (mass currently in the system)
  [mobile_W1_mass]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_W1
    block = '0'
    outputs = none
  []

  [W1_mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'mobile_W1_mass'
    pp_names = 'mobile_W1_mass'
    outputs = main
  []

  [mobile_Cu_mass]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_Cu
    block = '1'
    outputs = none
  []

  [Cu_mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'mobile_Cu_mass'
    pp_names = 'mobile_Cu_mass'
    outputs = main
  []

  [mobile_W2_mass]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_W2
    block = '2'
    outputs = none
  []

  [W2_mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'mobile_W2_mass'
    pp_names = 'mobile_W2_mass'
    outputs = main
  []

  # boundary fluxes
  [left_influx]
    type = ADSideDiffusiveFluxIntegral
    boundary = 'left'
    variable = C_M_W1
    diffusivity = ${W_D}
    outputs = main
  []

  [right_outflux]
    type = ADSideDiffusiveFluxIntegral
    boundary = 'right'
    variable = C_M_W2
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
  []

  #mass accumulated by directly integrating concentration in slab
  [mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'W1_mass_from_concentration + Cu_mass_from_concentration + W2_mass_from_concentration'
    pp_names = 'W1_mass_from_concentration Cu_mass_from_concentration W2_mass_from_concentration'
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
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  dtmax = 100
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
    cutback_factor = 0.9
  []
[]

[Outputs]
  file_base = W_Cu_W_${RUN_ID}
  [exodus]
    type = Exodus
  []
  [main]
    type = CSV
    file_base = W_Cu_W_${RUN_ID}_main
  []
  [vector]
    type = CSV
    execute_on = 'FINAL'
  []
[]