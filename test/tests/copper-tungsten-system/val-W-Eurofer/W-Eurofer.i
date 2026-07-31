# Permeation problem in a composite layer (W & Euro) adapted from Verification Problem #1e

# Numerical parameters
!include W-Eurofer.params
penalty = 1 #check conc&sol_ratio to ensure interface jump is properly enforced
end_time = 6000 #s
Euro_probe = ${units 1 mum -> m}
W_probe =  ${units 1 mum -> m}

[Mesh]
  coord_type = 'XYZ' # 3D Cartesian coordinates
  [generated]
    type = CartesianMeshGenerator
    dim = 1
    dx = '${W_thickness} ${Euro_thickness}'
    ix = '400 400'   # number of elements per subdomain
    subdomain_id = '0 1' 
  []
  [interface_W_to_Euro]
    type = SideSetsBetweenSubdomainsGenerator
    input = generated
    primary_block = '0' # W
    paired_block = '1' # Euro
    new_boundary = 'interface_W_to_Euro'
  []
  [interface_Euro_to_W]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface_W_to_Euro
    primary_block = '1' # Euro
    paired_block = '0' # W
    new_boundary = 'interface_Euro_to_W'
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
    prop_names = 'diffusivity_W diffusivity_Euro solubility_W solubility_Euro'
    prop_values = '${W_D} ${Euro_D} ${W_S} ${Euro_S}'
  []

  [converter_to_nonAD]
    type = MaterialADConverter
    ad_props_in = 'solubility_W solubility_Euro'
    reg_props_out = 'solubility_W_nonAD solubility_Euro_nonAD'
    outputs = 'exodus'
  []

  [interface_jump]
    type = SolubilityRatioMaterial
    solubility_primary = solubility_W
    solubility_secondary = solubility_Euro
    boundary = interface_W_to_Euro
    concentration_primary = C_M_W
    concentration_secondary = C_M_Euro
  []

  # [flux_dissociation_surface]
  #   type = ADParsedMaterial
  #   block = '0'
  #   coupled_variables = C_M_W
  #   property_name = 'flux_dissociation_surface'
  #   expression = '(2 * ${W_kd} * ${pressure}) - (2 * ${W_kr} * C_M_W^2)' #at/mum^2/s, * (1-C_M_W/${W_lattice_density})^2)
  # []

  # [flux_recombination_surface]
  #   type = ADParsedMaterial
  #   block = '1'
  #   coupled_variables = C_M_Euro
  #   property_name = 'flux_recombination_surface'
  #   expression = '- 2 * ${Euro_kr} * C_M_Euro ^ 2' #at/mum^2/s
  # []
[]

[Variables]
  [C_M_W] #mobile concentration of H in W, mol/m^3
    block = '0'
  []
  [C_M_Euro] #mobile concentration of H in Euro, mol/m^3 
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
  [Euro_diff]
    type = ADMatDiffusion
    block = '1'
    variable = C_M_Euro
    diffusivity = diffusivity_Euro
    extra_vector_tags = ref
  []
  [Euro_time]
    type = ADTimeDerivative
    block = '1'
    variable = C_M_Euro
    extra_vector_tags = ref
  []
[]

[BCs]
  [left_concentration]
    type = EquilibriumBC
    Ko = ${W_S} 
    boundary = left
    enclosure_var = ${pressure} #H2 domain
    temperature = ${T}
    variable = C_M_W
    p = 0.5 #Sievert's Law
  []
  # [left_dissociation]
  #   type = ADMatNeumannBC
  #   variable = C_M_W
  #   boundary = 'left'
  #   value = 1
  #   boundary_material = flux_dissociation_surface
  # []
  # [right_recombination]
  #   type = ADMatNeumannBC
  #   variable = C_M_Euro
  #   boundary = 'right'
  #   value = 1
  #   boundary_material = flux_recombination_surface
  # []
  #assume instantaneous recombination
  [right_concentration]
    type = ADDirichletBC
    boundary = right
    value = 0
    variable = C_M_Euro
  []
[]

[InterfaceKernels]
  # Penalized continuity with solubility jump at W/Euro interface
  [W-Euro-interface]
    type = ADPenaltyInterfaceDiffusion
    variable = C_M_W
    neighbor_var = C_M_Euro
    penalty = ${penalty}
    jump_prop_name = solubility_ratio
    boundary = 'interface_W_to_Euro'
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
  [Euro]
    type = LineValueSampler
    start_point = '${W_thickness} 0 0'
    end_point = '${fparse ${W_thickness} + ${Euro_thickness}} 0 0'
    num_points = 100
    sort_by = 'x'
    variable = C_M_Euro
    outputs = vector
  []
[]

[Postprocessors]

  # check concentration ratio against solubility to ensure penalty is properly enforced
  [Euro_S]
    type = SideAverageMaterialProperty
    property = solubility_Euro_nonAD
    boundary = interface_Euro_to_W
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [W_S]
    type = SideAverageMaterialProperty
    property = solubility_W_nonAD
    boundary = interface_Euro_to_W
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [gold_solubility_ratio]
    type = ParsedPostprocessor
    pp_names = 'W_S Euro_S'
    expression = 'W_S / Euro_S'
  []

  [W_interface]
    type = SideAverageValue
    boundary = interface_W_to_Euro
    variable = C_M_W
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [Euro_interface]
    type = SideAverageValue
    boundary = interface_Euro_to_W
    variable = C_M_Euro
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [variable_ratio]
    type = ParsedPostprocessor
    pp_names = 'W_interface Euro_interface'
    expression = 'W_interface / Euro_interface'
  []
  # conservation of mass check
  # mass accumulated in W coating
  [W_mass_from_concentration]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_W
    block = '0'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = main
  []
  # mass accmulated in Eurofer
  [Euro_mass_from_concentration]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_Euro
    block = '1'
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = main
  []

  #conservation of mass check on whole system
   [left_influx] # Influx at left boundary of W slab
    type = ADSideDiffusiveFluxIntegral
    boundary = 'left'
    variable = C_M_W
    diffusivity = ${W_D}
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []

  [right_outflux] #Outflux at right boundary of Euro slab
    type = ADSideDiffusiveFluxIntegral
    boundary = 'right'
    variable = C_M_Euro
    diffusivity = ${Euro_D}
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = main
  []

  [right_outflux_H2] #Outflux at right boundary of Euro slab
    type = ScalePostprocessor
    value = 'right_outflux'
    scaling_factor = 0.5
    outputs = main
  []

  [flux_difference]
    type = ParsedPostprocessor
    expression = '-left_influx - right_outflux' # negative sign on influx to account for outEuroard normal vector direction
    pp_names = 'left_influx right_outflux'
    outputs = none
  []

  # mass accmulated as inferred from boundary flux
  [mass_from_flux]
    type = TimeIntegratedPostprocessor
    value = flux_difference
    time_integration_scheme = IMPLICIT-EULER
    outputs = main
  []

  # mass accmulated as inferred from integrated concentration
  [mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'Euro_mass_from_concentration + W_mass_from_concentration'
    pp_names = 'Euro_mass_from_concentration W_mass_from_concentration'
    outputs = main
  []

   # Used to obtain varying concentration with time at a point
  [concentration_at_x_W]
    type = PointValue
    variable = C_M_W
    point = '${W_probe} 0 0'
  []

  [concentration_at_x_Euro]
    type = PointValue
    variable = C_M_Euro
    point = '${fparse ${W_thickness} + ${Euro_probe}} 0 0'
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
  end_time = ${end_time}
  solve_type = NEWTON
  line_search = 'bt'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-10
  # steady_state_tolerance = 1e-6
  # steady_state_detection = true
  # steady_state_start_time = 1.0
  automatic_scaling = true
  compute_scaling_once = false
  dtmax = 5
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-8
    optimal_iterations = 6
    growth_factor = 1.15
    cutback_factor = 0.9
  []
[]

[Outputs]
  file_base = W_Eurofer_${pressure}
  [exodus]
    type = Exodus
  []
  [main]
    type = CSV
    file_base = W_Eurofer_${pressure}_main
  []
  [vector]
    type = CSV
    execute_on = 'FINAL'
    # sync_times = ${simulation_time}
    # sync_only = true
  []
[]
