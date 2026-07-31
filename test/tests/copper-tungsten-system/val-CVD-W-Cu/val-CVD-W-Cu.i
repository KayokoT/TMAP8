# Permeation problem in a three-layer composite (W / entry-Cu / Cu) to compare against CVD W experimental study (Sun 2025)
# As mentioned by Sun, the entry-Cu is necessary to model between the pure W and pure Cu due to 
# differing solution energies caused by Cu atoms mixing in the W lattice.

# Numerical parameters
!include val-CVD-W-Cu.params
penalty = 1e6 #check conc&sol_ratio to ensure BOTH interface jumps are properly enforced

[Mesh]
  coord_type = 'XYZ'
  [generated]
    type = CartesianMeshGenerator
    dim = 1
    dx = '${W_thickness} ${Cu_entry_thickness} ${Cu_thickness}'
    ix = '800 400 400'   # number of elements per subdomain
    subdomain_id = '0 1 2'   # W /entry-Cu/bulk Cu
  []
  [interface_W_to_CuEntry]
    type = SideSetsBetweenSubdomainsGenerator
    input = generated
    primary_block = '0'   # W
    paired_block = '1'    #  entry-Cu
    new_boundary = 'interface_W_to_CuEntry'
  []
  [interface_CuEntry_to_W]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface_W_to_CuEntry
    primary_block = '1'   # entry-Cu
    paired_block = '0'    # W
    new_boundary = 'interface_CuEntry_to_W'
  []
  [interface_CuEntry_to_Cu]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface_CuEntry_to_W
    primary_block = '1'   # entry-Cu
    paired_block = '2'    # bulk Cu
    new_boundary = 'interface_CuEntry_to_Cu'
  []
  [interface_Cu_to_CuEntry]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface_CuEntry_to_Cu
    primary_block = '2'   # bulk-Cu
    paired_block = '1'    # entry Cu
    new_boundary = 'interface_Cu_to_CuEntry'
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
    prop_names = 'diffusivity_W diffusivity_CuEntry diffusivity_Cu solubility_W solubility_CuEntry solubility_Cu'
    prop_values = '${W_D} ${W_D} ${Cu_D} ${W_S} ${Cu_entry_S} ${Cu_S}' #diffusivity of W = diffusivity of CuEntry
  []

  [converter_to_nonAD]
    type = MaterialADConverter
    ad_props_in = 'solubility_W solubility_CuEntry solubility_Cu'
    reg_props_out = 'solubility_W_nonAD solubility_CuEntry_nonAD solubility_Cu_nonAD'
    outputs = 'exodus'
  []

  [interface_jump_W_CuEntry]
    type = SolubilityRatioMaterial
    solubility_primary = solubility_W
    solubility_secondary = solubility_CuEntry
    boundary = interface_W_to_CuEntry
    concentration_primary = C_M_W
    concentration_secondary = C_M_CuEntry
  []

  [interface_jump_CuEntry_Cu]
    type = SolubilityRatioMaterial
    solubility_primary = solubility_CuEntry
    solubility_secondary = solubility_Cu
    boundary = interface_CuEntry_to_Cu
    concentration_primary = C_M_CuEntry
    concentration_secondary = C_M_Cu
  []

  [flux_left_surface]
    type = ADParsedMaterial
    block = '0'
    coupled_variables = C_M_W
    property_name = 'flux_left_surface'
    expression = '(${W_J0_coef} * ${J0}) - (${WCu_kr} * C_M_W ^ 2)' #at/mum^2/s
  []
[]

[Variables]
  [C_M_W] #mobile concentration of H in W, at/mum^3
    block = '0'
  []
  [C_M_CuEntry] #mobile concentration of H in entry-Cu layer, at/mum^3
    block = '1'
  []
  [C_M_Cu] #mobile concentration of H in bulk Cu, at/mum^3
    block = '2'
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
  [CuEntry_diff]
    type = ADMatDiffusion
    block = '1'
    variable = C_M_CuEntry
    diffusivity = diffusivity_CuEntry
    extra_vector_tags = ref
  []
  [CuEntry_time]
    type = ADTimeDerivative
    block = '1'
    variable = C_M_CuEntry
    extra_vector_tags = ref
  []
  [Cu_diff]
    type = ADMatDiffusion
    block = '2'
    variable = C_M_Cu
    diffusivity = diffusivity_Cu
    extra_vector_tags = ref
  []
  [Cu_time]
    type = ADTimeDerivative
    block = '2'
    variable = C_M_Cu
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
    boundary_material = flux_left_surface
    extra_vector_tags = 'ref'
  []
  # diffusion-limited downstream
  [right]
    type = ADDirichletBC
    variable = C_M_Cu
    boundary = 'right'
    value = 0
  []
[]

[InterfaceKernels]
  # Penalized continuity with solubility jump at W / entry-Cu interface
  [W_to_CuEntry_interface]
    type = ADPenaltyInterfaceDiffusion
    variable = C_M_W
    neighbor_var = C_M_CuEntry
    penalty = ${penalty}
    jump_prop_name = solubility_ratio
    boundary = 'interface_W_to_CuEntry'
    extra_vector_tags = 'ref'
  []
  # Penalized continuity with solubility jump at entry-Cu / bulk-Cu interface
  [CuEntry_to_Cu_interface]
    type = ADPenaltyInterfaceDiffusion
    variable = C_M_CuEntry
    neighbor_var = C_M_Cu
    penalty = ${penalty}
    jump_prop_name = solubility_ratio
    boundary = 'interface_CuEntry_to_Cu'
    extra_vector_tags = 'ref'
  []
[]

[VectorPostprocessors]
  [W]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '${W_thickness} 0 0'
    num_points = 100
    sort_by = 'x'
    variable = C_M_W
    outputs = vector
  []
  [CuEntry]
    type = LineValueSampler
    start_point = '${W_thickness} 0 0'
    end_point = '${fparse ${W_thickness} + ${Cu_entry_thickness}} 0 0'
    num_points = 100
    sort_by = 'x'
    variable = C_M_CuEntry
    outputs = vector
  []
  [Cu]
    type = LineValueSampler
    start_point = '${fparse ${W_thickness} + ${Cu_entry_thickness}} 0 0'
    end_point = '${fparse ${W_thickness} + ${Cu_entry_thickness} + ${Cu_thickness}} 0 0'
    num_points = 100
    sort_by = 'x'
    variable = C_M_Cu
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
  # check concentration ratio against solubility to ensure penalty is properly enforced
  # W/Cu-entry interface
  [W_S]
    type = SideAverageMaterialProperty
    property = solubility_W_nonAD
    boundary = interface_W_to_CuEntry
  []

  [CuEntry_S]
    type = SideAverageMaterialProperty
    property = solubility_CuEntry_nonAD
    boundary = interface_W_to_CuEntry
  []

  [gold_solubility_ratio_1]
    type = ParsedPostprocessor
    pp_names = 'CuEntry_S W_S'
    expression = 'CuEntry_S / W_S'
  []

  [CuEntry_interface_1]
    type = SideAverageValue
    boundary = interface_CuEntry_to_W
    variable = C_M_CuEntry
  []

  [W_interface]
    type = SideAverageValue
    boundary = interface_W_to_CuEntry
    variable = C_M_W
  []

  [variable_ratio_1]
    type = ParsedPostprocessor
    pp_names = 'CuEntry_interface_1 W_interface'
    expression = 'CuEntry_interface_1 / W_interface'
  []
  # Cu-entry/Cu interface
  [Cu_S]
    type = SideAverageMaterialProperty
    property = solubility_Cu_nonAD
    boundary = interface_CuEntry_to_Cu
  []

  [gold_solubility_ratio_2]
    type = ParsedPostprocessor
    pp_names = 'Cu_S CuEntry_S'
    expression = 'Cu_S / CuEntry_S'
  []

  [Cu_interface]
    type = SideAverageValue
    boundary = interface_Cu_to_CuEntry
    variable = C_M_Cu
  []

  [CuEntry_interface_2]
    type = SideAverageValue
    boundary = interface_CuEntry_to_Cu
    variable = C_M_CuEntry
  []

  [variable_ratio_2]
    type = ParsedPostprocessor
    pp_names = 'Cu_interface CuEntry_interface_2'
    expression = 'Cu_interface / CuEntry_interface_2'
  []
  # conservation of mass check
  # inventory (mass currently in the system)
  [mobile_Cu_mass]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_Cu
    block = '2'
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
    block = '0'
    outputs = none
  []

  [W_mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'mobile_W_mass'
    pp_names = 'mobile_W_mass'
    outputs = main
  []
  [mobile_CuEntry_mass]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_CuEntry
    block = '1'
    outputs = none
  []

  [CuEntry_mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'mobile_CuEntry_mass'
    pp_names = 'mobile_CuEntry_mass'
    outputs = main
  []

  # boundary fluxes
  [left_influx]
    type = ADSideDiffusiveFluxIntegral
    boundary = 'left'
    variable = C_M_W
    diffusivity = ${W_D}
    outputs = main
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
    time_integration_scheme = trapezoidal-rule
    outputs = main
  []

  #mass accumulated by directly integrating concentration in slab
  [mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'W_mass_from_concentration + CuEntry_mass_from_concentration + Cu_mass_from_concentration'
    pp_names = 'W_mass_from_concentration CuEntry_mass_from_concentration Cu_mass_from_concentration'
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
  file_base = val_W_entryCu_Cu_${T}
  [exodus]
    type = Exodus
  []
  [main]
    type = CSV
    file_base = val_W_entryCu_Cu_${T}_main
  []
  [vector]
    type = CSV
    execute_on = 'FINAL'
  []
[]