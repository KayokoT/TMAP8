# Permeation problem in a composite layer (Cu-W-Ar) adapted from Verification Problem #1e

# Numerical parameters
!include WCuAr-atmum3.params
W_probe = ${units 1 mum} #probe in W layer, measured from right boundary of Cu
Cu_probe = ${units 1 mum} #probe in Cu layer, measured from left side of Cu
Ar_probe = ${units 1 mum} #probe in air, measured from right boundary of W
penalty = 1e4 #check conc&sol_ratio to ensure interface jump is properly enforced
RUN_ID = run01 #default, change this when running by adding RUN_ID=run**
Ar_thickness = '${units 5e-5 m -> mum}'

[Mesh]
  coord_type = 'XYZ'
  [generated]
    type = CartesianMeshGenerator
    dim = 1
    dx = '${Cu_thickness} ${W_thickness} ${Ar_thickness}'
    ix = '400 200 400'   # number of elements per subdomain
    subdomain_id = '0 1 2' # Cu W Ar
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
  [interface_W_to_Ar]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface_W_to_Cu
    primary_block = '1' # W
    paired_block = '2' # Ar
    new_boundary = 'interface_W_to_Ar'
  []
  [interface_Ar_to_W]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface_W_to_Ar
    primary_block = '2'  # Ar
    paired_block = '1'   # W
    new_boundary = 'interface_Ar_to_W'
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
    prop_names = 'diffusivity_Cu diffusivity_W diffusivity_Ar solubility_Cu solubility_W'
    prop_values = '${Cu_D} ${W_D} ${Ar_D} ${Cu_S} ${W_S}'
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
[]

[Variables]
  [C_M_Cu] #mobile concentration of H in Cu, at/mum^3
    block = '0'
  []
  [C_M_W] #mobile concentration of H in W, at/mum^3 
    block = '1'
  []
  [C_M_Ar] #mobile concentration of H in Ar, at/mum^3 
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
  [Ar_diff]
    type = ADMatDiffusion
    block = '2'
    variable = C_M_Ar
    diffusivity = diffusivity_Ar
    extra_vector_tags = ref
  []
  [Ar_time]
    type = ADTimeDerivative
    block = '2'
    variable = C_M_Ar
    extra_vector_tags = ref
  []
[]

[BCs]
  [left_flux]
    type = EquilibriumBC
    Ko = ${Cu_S} 
    activation_energy = 0
    boundary = left
    enclosure_var = ${pressure} #H2 domain
    temperature = ${T}
    variable = C_M_Cu
    p = 0.5 #Sievert's Law
  []
  [right]
    type = ADNeumannBC
    variable = C_M_Ar
    value = 0
    boundary = 'right'
    extra_vector_tags = ref
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
  []
  #sorption law boundary condition
  [W-Ar-interface]
    type = ADInterfaceSorption
    K0 = '${W_S}'
    Ea = 0
    diffusivity = diffusivity_W
    n_sorption = 0.5 #Sievert's Law
    variable = C_M_W
    neighbor_var =  C_M_Ar
    temperature = ${T}
    boundary = 'interface_W_to_Ar'
    unit_scale_neighbor = ${fparse 1e18/ ${NA}} #converts Cg from at/mum^3 -> mol/m^3 to match R
  []
[]

[VectorPostprocessors]
  [Cu]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '${Cu_thickness} 0 0'
    num_points = 100
    sort_by = 'x'
    variable = C_M_Cu
    outputs = vector
  []
  [W]
    type = LineValueSampler
    start_point = '${Cu_thickness} 0 0'
    end_point = '${fparse ${Cu_thickness} + ${W_thickness}} 0 0'
    num_points = 100
    sort_by = 'x'
    variable = C_M_W
    outputs = vector
  []
  [Ar]
    type = LineValueSampler
    start_point = '${fparse ${Cu_thickness} + ${W_thickness}} 0 0'
    end_point = '${fparse ${Cu_thickness} + ${W_thickness} + ${Ar_thickness}} 0 0'
    num_points = 100
    sort_by = 'x'
    variable = C_M_Ar
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

  [concentration_at_x_Ar]
    type = PointValue
    variable = C_M_Ar
    point = '${fparse ${Cu_thickness} + ${W_thickness} + ${Ar_probe}} 0 0'
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

  # conservation of mass check - Cu slab
  [Cu_left_influx] # Influx at left boundary
    type = ADSideDiffusiveFluxIntegral
    boundary = 'left'
    variable = C_M_Cu
    diffusivity = ${Cu_D}
    outputs = none
  []

  [Cu_right_outflux] #Outflux at right boundary of Cu slab
    type = ADSideDiffusiveFluxIntegral
    boundary = interface_Cu_to_W
    variable = C_M_Cu
    diffusivity = ${Cu_D}
    outputs = none
  []

  [Cu_flux_difference]
    type = ParsedPostprocessor
    expression = '-Cu_left_influx - Cu_right_outflux' # negative sign on influx to account for outward normal vector direction
    pp_names = 'Cu_left_influx Cu_right_outflux'
    outputs = none
  []
  # mass accumulated as inferred from boundary flux
  [Cu_mass_from_flux]
    type = TimeIntegratedPostprocessor
    value = Cu_flux_difference
    time_integration_scheme = IMPLICIT-EULER
    outputs = main
  []
  # mass accumulated by directly integrating concentration in slab
  [Cu_mass_from_concentration]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_Cu
    block = '0'
    outputs = main
  []

  # conservation of mass check - W slab
  [W_left_influx] # Influx at left boundary of W slab
    type = ADSideDiffusiveFluxIntegral
    boundary = interface_W_to_Cu
    variable = C_M_W
    diffusivity = ${W_D}
    outputs = main
  []

  [W_right_outflux] #Outflux at right boundary of W slab
    type = ADSideDiffusiveFluxIntegral
    boundary = interface_W_to_Ar
    variable = C_M_W
    diffusivity = ${W_D}
    outputs = none
  []

  [W_flux_difference]
    type = ParsedPostprocessor
    expression = '-W_left_influx - W_right_outflux' # negative sign on influx to account for outward normal vector direction
    pp_names = 'W_left_influx W_right_outflux'
    outputs = none
  []
  # mass accumulation as inferred from boundary flux
  [W_mass_from_flux] 
    type = TimeIntegratedPostprocessor
    value = W_flux_difference
    time_integration_scheme = IMPLICIT-EULER
    outputs = main
  []
  # mass accumulated by directly integrating concentration in slab
  [W_mass_from_concentration]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_W
    block = '1'
    outputs = main
  []

  # conservation of mass check - Ar layer
  [Ar_left_influx] # Influx at left boundary of Ar layer
    type = ADSideDiffusiveFluxIntegral
    boundary = interface_Ar_to_W
    variable = C_M_Ar
    diffusivity = ${Ar_D}
    outputs = none
  []

  [Ar_right_outflux] #Outflux at right boundary of Ar layer
    type = ADSideDiffusiveFluxIntegral
    boundary = 'right'
    variable = C_M_Ar
    diffusivity = ${Ar_D}
    outputs = none
  []

  [Ar_flux_difference]
    type = ParsedPostprocessor
    expression = '-Ar_left_influx - Ar_right_outflux' # negative sign on influx to account for outward normal vector direction
    pp_names = 'Ar_left_influx Ar_right_outflux'
    outputs = none
  []
  # mass accumulation as inferred from boundary flux
  [Ar_mass_from_flux] 
    type = TimeIntegratedPostprocessor
    value = Ar_flux_difference
    time_integration_scheme = IMPLICIT-EULER
    outputs = main
  []
  # mass accumulated by directly integrating concentration in slab
  [Ar_mass_from_concentration]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_Ar
    block = '2'
    outputs = main
  []

  #conservation of mass check on whole system
   [left_influx] # Influx at left boundary of W slab
    type = ADSideDiffusiveFluxIntegral
    boundary = 'left'
    variable = C_M_Cu
    diffusivity = ${Cu_D}
    outputs = none
  []

  [right_outflux] #Outflux at right boundary of W slab
    type = ADSideDiffusiveFluxIntegral
    boundary = 'right'
    variable = C_M_Ar
    diffusivity = ${Ar_D}
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
    outputs = main
  []
  # mass accumulated as inferred from boundary flux
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
  nl_abs_tol = 1e-7
  dtmax = 1e-3
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
  file_base = ${RUN_ID}
  [exodus]
    type = Exodus
  []
  [main]
    type = CSV
    file_base = ${RUN_ID}_main
  []
  [vector]
    type = CSV
    execute_on = 'FINAL'
    # sync_times = ${simulation_time}
    # sync_only = true
  []
[]
