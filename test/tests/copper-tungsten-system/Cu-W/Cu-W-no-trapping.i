# Permeation problem in a composite layer (W & Cu) adapted from Verification Problem #1e
# Cu faces the H source
# No trapping
# BCs: Sieverts (left) & recombination (right)

!include WCu-atmum3.params
W_probe = ${units 1 mum} #probe in W layer, measured from right boundary of Cu
Cu_probe = ${units 1 mum} #probe in Cu layer, measured from left side of Cu
penalty = 1e9 #check conc&sol_ratio to ensure interface jump is properly enforced
RUN_ID = run01 #default, change this when running by adding RUN_ID=run**

[Mesh]
  coord_type = 'XYZ' # 3D Cartesian coordinates
  [generated]
    type = CartesianMeshGenerator
    dim = 1
    dx = '${Cu_thickness} ${W_thickness}'
    ix = '400 400'   # number of elements per subdomain
    subdomain_id = '0 1' 
  []
  [interface_Cu_to_W]
    type = SideSetsBetweenSubdomainsGenerator
    input = generated
    primary_block = '0' # Cu
    paired_block = '1' # W
    new_boundary = 'interface_Cu_to_W'
  []
  [interface_W_to_Cu]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface_Cu_to_W
    primary_block = '1' # W
    paired_block = '0' # Cu
    new_boundary = 'interface_W_to_Cu'
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
    prop_names = 'diffusivity_Cu diffusivity_W solubility_Cu solubility_W'
    prop_values = '${Cu_D} ${W_D} ${Cu_S} ${W_S}'
  []

  [converter_to_nonAD]
    type = MaterialADConverter
    ad_props_in = 'solubility_Cu solubility_W'
    reg_props_out = 'solubility_Cu_nonAD solubility_W_nonAD'
    outputs = 'exodus'
  []

  [interface_jump]
    type = SolubilityRatioMaterial
    solubility_primary = solubility_Cu
    solubility_secondary = solubility_W
    boundary = interface_Cu_to_W
    concentration_primary = C_M_Cu
    concentration_secondary = C_M_W
  []

  # if including dissociation + recombination on the left
  # [flux_dissociation_recombination_surface]
  #   type = ADParsedMaterial
  #   block = '0'
  #   coupled_variables = C_M_Cu
  #   property_name = 'flux_dissociation_recombination_surface'
  #   expression = '(2 * ${Cu_kd} * ${pressure}) - (2 * ${Cu_kr} * C_M_Cu^2)' #at/mum^2/s
  # []

  [flux_recombination_surface]
    type = ADParsedMaterial
    block = '1'
    coupled_variables = C_M_W
    property_name = 'flux_recombination_surface'
    expression = '-2 * ${W_kr} * C_M_W ^ 2' #at/mum^2/s
  []
[]

[Variables]
  [C_M_Cu] #mobile concentration of H in Cu, at/mum^3
    block = '0'
  []
  [C_M_W] #mobile concentration of H in W, at/mum^3 
    block = '1'
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
  [W_diff]
    type = ADMatDiffusion
    block = '1'
    variable = C_M_W
    diffusivity = diffusivity_W
    extra_vector_tags = ref
  []
  [W_time]
    type = ADTimeDerivative
    block = '1'
    variable = C_M_W
    extra_vector_tags = ref
  []
[]

[BCs]
  [left_flux]
    type = EquilibriumBC
    Ko = ${Cu_S} 
    boundary = left
    enclosure_var = ${pressure} #H2 domain
    temperature = ${T}
    variable = C_M_Cu
    p = 0.5 #Sievert's Law
  []
  # if including dissociation + recombination on the left
  # [left_dissociation_recombination]
  #   type = ADMatNeumannBC
  #   variable = C_M_Cu
  #   boundary = 'left'
  #   value = 1
  #   boundary_material = flux_dissociation_recombination_surface
  # []
  [right_recombination]
    type = ADMatNeumannBC
    variable = C_M_W
    boundary = 'right'
    value = 1
    boundary_material = flux_recombination_surface
    extra_vector_tags = 'ref'
  []

[]

[InterfaceKernels]
  # Penalized continuity with solubility jump at Cu/W interface
  [Cu-W-interface]
    type = ADPenaltyInterfaceDiffusion
    variable = C_M_Cu
    neighbor_var = C_M_W
    penalty = ${penalty}
    jump_prop_name = solubility_ratio
    boundary = 'interface_Cu_to_W'
    extra_vector_tags = 'ref'
  []
[]

[VectorPostprocessors]
  [Cu]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '${Cu_thickness} 0 0'
    num_points = 400
    sort_by = 'x'
    variable = C_M_Cu
    outputs = vector
  []
  [W]
    type = LineValueSampler
    start_point = '${Cu_thickness} 0 0'
    end_point = '${fparse ${Cu_thickness} + ${W_thickness}} 0 0'
    num_points = 400
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
    point = '${fparse ${Cu_thickness} + ${W_probe}} 0 0'
  []

  [concentration_at_x_Cu]
    type = PointValue
    variable = C_M_Cu
    point = '${Cu_probe} 0 0'
  []

  [outflux]
    type = SideDiffusiveFluxAverage
    boundary = 'right'
    diffusivity = ${W_D}
    variable = C_M_W
  []

  # check concentration ratio against solubility to ensure penalty is properly enforced
  [W_S]
    type = SideAverageMaterialProperty
    property = solubility_W_nonAD
    boundary = interface_W_to_Cu
  []

  [Cu_S]
    type = SideAverageMaterialProperty
    property = solubility_Cu_nonAD
    boundary = interface_W_to_Cu
  []

  [gold_solubility_ratio]
    type = ParsedPostprocessor
    pp_names = 'Cu_S W_S'
    expression = 'Cu_S / W_S'
  []

  [Cu_interface]
    type = SideAverageValue
    boundary = interface_Cu_to_W
    variable = C_M_Cu
  []

  [W_interface]
    type = SideAverageValue
    boundary = interface_W_to_Cu
    variable = C_M_W
  []

  [variable_ratio]
    type = ParsedPostprocessor
    pp_names = 'Cu_interface W_interface'
    expression = 'Cu_interface / W_interface'
  []

  # inventory (mass currently in the system)
  [mobile_Cu_mass]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_Cu
    block = '0'
    outputs = none
  []

  [Cu_mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'mobile_Cu_mass'
    pp_names = 'mobile_Cu_mass'
    outputs = main
  []

  [mobile_W_mass]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_W
    block = '1'
    outputs = none
  []

  [W_mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'mobile_W_mass'
    pp_names = 'mobile_W_mass'
    outputs = main
  []

  # boundary fluxes
  [left_influx]
    type = ADSideDiffusiveFluxIntegral
    boundary = 'left'
    variable = C_M_Cu
    diffusivity = ${Cu_D}
    outputs = main
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
  []

  #mass accumulated by directly integrating concentration in slab
  [mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'W_mass_from_concentration + Cu_mass_from_concentration'
    pp_names = 'W_mass_from_concentration Cu_mass_from_concentration'
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
  solve_type = NEWTON
  line_search = 'bt'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  dtmax = 1e-2
  steady_state_tolerance = 1e-6
  steady_state_detection = true
  steady_state_start_time = 1.0
  automatic_scaling = true
  compute_scaling_once = false
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-8
    optimal_iterations = 6
    growth_factor = 1.15
    cutback_factor = 0.9
  []
[]

[Outputs]
  file_base = Cu_W_${RUN_ID}
  [exodus]
    type = Exodus
  []
  [main]
    type = CSV
    file_base = Cu_W_${RUN_ID}_main
  []
  [vector]
    type = CSV
    execute_on = 'FINAL'
  []
[]
