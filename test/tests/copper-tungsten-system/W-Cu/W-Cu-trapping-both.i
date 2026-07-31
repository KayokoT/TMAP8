# Permeation problem in a composite layer (W & Cu) adapted from Verification Problem #1e
# W faces the H source
# 3 traps - 2 in W and 1 in Cu
# Includes commented out code for a 2nd trap in Cu
# BCs: Sieverts (left) & recombination (right)

!include WCu-atmum3.params
Cu_probe = ${units 1 mum} #probe in Cu layer, measured from right boundary of W
W_probe = ${units 1 mum} #probe in W layer, measured from left side of W
penalty = 1e10 #check conc&sol_ratio to ensure interface jump is properly enforced
RUN_ID = run01 #default, change this when running by adding RUN_ID=run**
trap_per_free_1 = 1e3
trap_per_free_2 = 1e3
trap_per_free_3 = 1000

[Mesh]
  coord_type = 'XYZ' 
  [generated]
    type = CartesianMeshGenerator
    dim = 1
    dx = '${W_thickness} ${Cu_thickness}'
    ix = '400 400' 
    subdomain_id = '0 1' 
  []
  [interface_W_to_Cu]
    type = SideSetsBetweenSubdomainsGenerator
    input = generated
    primary_block = '0' # W
    paired_block = '1' # Cu
    new_boundary = 'interface_W_to_Cu'
  []
  [interface_Cu_to_W]
    type = SideSetsBetweenSubdomainsGenerator
    input = interface_W_to_Cu
    primary_block = '1' # Cu
    paired_block = '0' # W
    new_boundary = 'interface_Cu_to_W'
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

  [interface_jump]
    type = SolubilityRatioMaterial
    solubility_primary = solubility_W
    solubility_secondary = solubility_Cu
    boundary = interface_W_to_Cu
    concentration_primary = C_M_W
    concentration_secondary = C_M_Cu
  []

  # if including dissociation + recombination on the left
  # [flux_dissociation_recombination_surface]
  #   type = ADParsedMaterial
  #   block = '0'
  #   coupled_variables = C_M_W
  #   property_name = 'flux_dissociation_recombination_surface'
  #   expression = '(2 * ${W_kd} * ${pressure}) - (2 * ${W_kr} * C_M_W^2)' #at/mum^2/s
  # []

  [flux_recombination_surface]
    type = ADParsedMaterial
    block = '1'
    coupled_variables = C_M_Cu
    property_name = 'flux_recombination_surface'
    expression = '- 2 * ${Cu_kr} * C_M_Cu ^ 2' #H2 molecules/mum^2/s
  []
[]

[Variables]
  [C_M_W] #mobile concentration of H in W, at/mum^3
    block = '0'
  []
  [C_M_Cu] #mobile concentration of H in Cu, at/mum^3 
    block = '1'
  []
  [C_T_W_1] #concentration of H in W in trap 1 scaled by trap_per_free, at/mum^3
    block = '0'
  []
  [C_T_Cu_1] #concentration of H in Cu in trap 1, at/mum^3 
    block = '1'                      
  []
  [C_T_W_2] #concentration of H in W in trap 2, at/mum^3
    block = '0'
  []
  # [C_T_Cu_2] #concentration of H in Cu in trap 2, at/mum^3 
  #   block = '1'                      
  # []
[]

[AuxVariables]
  [pressure]
    family = SCALAR
    initial_condition = ${pressure}
  []
  [W_empty_sites_1]
     block = '0'
  []
  [W_trapped_sites_1]
     block = '0'
  []
  [W_empty_sites_2]
     block = '0'
  []
  [W_trapped_sites_2]
     block = '0'
  []
  [W_total_sites]
     block = '0'
  []
  [Cu_empty_sites_1]
     block = '1'
  []
  [Cu_trapped_sites_1]
     block = '1'
  []
  # [Cu_empty_sites_2]
  #    block = '1'
  # []
  # [Cu_trapped_sites_2]
  #    block = '1'
  # []
  [Cu_total_sites]
     block = '1'
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
  [W_coupled_time_1]
    type = ADScaledCoupledTimeDerivative
    block = '0'
    variable = C_M_W
    v = C_T_W_1
    mat_prop = ${trap_per_free_1}
    extra_vector_tags = ref
  []
  [W_coupled_time_2]
    type = ADScaledCoupledTimeDerivative
    block = '0'
    variable = C_M_W
    v = C_T_W_2
    mat_prop = ${trap_per_free_2}
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
  [Cu_coupled_time_1]
    type = ADScaledCoupledTimeDerivative
    block = '1'
    variable = C_M_Cu
    v = C_T_Cu_1
    mat_prop = ${trap_per_free_3}
    extra_vector_tags = ref
  []
  # [Cu_coupled_time_2]
  #   type = ADScaledCoupledTimeDerivative
  #   block = '1'
  #   variable = C_M_Cu
  #   v = C_T_Cu_1
  #   mat_prop = ${trap_per_free_4}
  #   extra_vector_tags = ref
  # []
[]

[NodalKernels]
  [W_time_1]
    type = TimeDerivativeNodalKernel
    block = '0'
    variable = C_T_W_1
    extra_vector_tags = ref
  []
  [W_trapping_1]
    type = TrappingNodalKernel
    block = '0'
    variable = C_T_W_1
    alpha_t = ${W_trapping_prefactor}
    N = '${W_lattice_density}'
    Ct0 = ${W_trapping_site_fraction_1}
    trapping_energy = ${W_migration_energy}
    mobile_concentration = 'C_M_W'
    temperature = ${T}
    trap_per_free = ${trap_per_free_1}
    extra_vector_tags = ref
  []
  [W_release_1]
    type = ReleasingNodalKernel
    block = '0'
    alpha_r = ${W_detrapping_prefactor}
    temperature = ${T}
    detrapping_energy = ${W_detrapping_energy_1}
    variable = C_T_W_1
    extra_vector_tags = ref
  []
  [W_time_2]
    type = TimeDerivativeNodalKernel
    block = '0'
    variable = C_T_W_2
    extra_vector_tags = ref
  []
  [W_trapping_2]
    type = TrappingNodalKernel
    block = '0'
    variable = C_T_W_2
    alpha_t = ${W_trapping_prefactor}
    N = '${W_lattice_density}'
    Ct0 = ${W_trapping_site_fraction_2}
    trapping_energy = ${W_migration_energy}
    mobile_concentration = 'C_M_W'
    temperature = ${T}
    trap_per_free = ${trap_per_free_2}
    extra_vector_tags = ref
  []
  [W_release_2]
    type = ReleasingNodalKernel
    block = '0'
    alpha_r = ${W_detrapping_prefactor}
    temperature = ${T}
    detrapping_energy = ${W_detrapping_energy_2}
    variable = C_T_W_2
    extra_vector_tags = ref
  []
  [Cu_time_1]
    type = TimeDerivativeNodalKernel
    block = '1'
    variable = C_T_Cu_1
    extra_vector_tags = ref
  []
  [Cu_trapping_1]
    type = TrappingNodalKernel
    block = '1'
    variable = C_T_Cu_1
    alpha_t = ${Cu_trapping_prefactor}
    N = '${Cu_lattice_density}'
    Ct0 = ${Cu_trapping_site_fraction_1}
    trapping_energy = ${Cu_migration_energy}
    mobile_concentration = 'C_M_Cu'
    temperature = ${T}
    trap_per_free = ${trap_per_free_3}
    extra_vector_tags = ref
  []
  [Cu_release_1]
    type = ReleasingNodalKernel
    block = '1'
    alpha_r = ${Cu_detrapping_prefactor}
    temperature = ${T}
    detrapping_energy = ${Cu_detrapping_energy_1}
    variable = C_T_Cu_1
    extra_vector_tags = ref
  []
  # [Cu_time_2]
  #   type = TimeDerivativeNodalKernel
  #   block = '1'
  #   variable = C_T_Cu_2
  # []
  # [Cu_trapping_2]
  #   type = TrappingNodalKernel
  #   block = '1'
  #   variable = C_T_Cu_2
  #   alpha_t = ${Cu_trapping_prefactor}
  #   N = '${Cu_lattice_density}'
  #   Ct0 = ${Cu_trapping_site_fraction_2}
  #   mobile_concentration = 'C_M_Cu'
  #   temperature = ${T}
  #   trap_per_free = ${trap_per_free_4}
  #   extra_vector_tags = ref
  # []
  # [Cu_release_2]
  #   type = ReleasingNodalKernel
  #   block = '1'
  #   alpha_r = ${Cu_detrapping_prefactor}
  #   temperature = ${T}
  #   detrapping_energy = ${Cu_detrapping_energy_2}
  #   variable = C_T_Cu_2
  # []
[]

[AuxKernels]
  [W_empty_sites_1]
    variable = W_empty_sites_1
    type = EmptySitesAux
    N = ${W_lattice_density}
    Ct0 = ${W_trapping_site_fraction_1}
    trap_per_free = ${trap_per_free_1}
    trapped_concentration_variables = C_T_W_1
    block = '0'
  []
  [W_trapped_sites_1]
    variable = W_trapped_sites_1
    type = NormalizationAux
    normal_factor = ${trap_per_free_1}
    source_variable = C_T_W_1
    block = '0'
  []
  [W_empty_sites_2]
    variable = W_empty_sites_2
    type = EmptySitesAux
    N = ${W_lattice_density}
    Ct0 = ${W_trapping_site_fraction_2}
    trap_per_free = ${trap_per_free_2}
    trapped_concentration_variables = C_T_W_2
  []
  [W_trapped_sites_2]
    variable = W_trapped_sites_2
    type = NormalizationAux
    normal_factor = ${trap_per_free_2}
    source_variable = C_T_W_2
  []
  [W_total_sites]
    variable = W_total_sites
    type = ParsedAux
    expression = 'W_trapped_sites_1 + W_empty_sites_1 + W_trapped_sites_2 + W_empty_sites_2'
    coupled_variables = 'W_trapped_sites_1 W_empty_sites_1 W_trapped_sites_2 W_empty_sites_2'
    block = '0'
  []
  [Cu_empty_sites_1]
    variable = Cu_empty_sites_1
    type = EmptySitesAux
    N = ${Cu_lattice_density}
    Ct0 = ${Cu_trapping_site_fraction_1}
    trap_per_free = ${trap_per_free_3}
    trapped_concentration_variables = C_T_Cu_1
  []
  [Cu_trapped_sites_1]
    variable = Cu_trapped_sites_1
    type = NormalizationAux
    normal_factor = ${trap_per_free_3}
    source_variable = C_T_Cu_1
  []
  # [Cu_empty_sites_2]
  #   variable = Cu_empty_sites_2
  #   type = EmptySitesAux
  #   N = ${Cu_lattice_density}
  #   Ct0 = ${Cu_trapping_site_fraction_2}
  #   trap_per_free = ${trap_per_free_4}
  #   trapped_concentration_variables = C_T_Cu_2
  # []
  # [Cu_trapped_sites_2]
  #   variable = Cu_trapped_sites_2
  #   type = NormalizationAux
  #   normal_factor = ${trap_per_free_4}
  #   source_variable = C_T_Cu_2
  # []
  [Cu_total_sites]
    variable = Cu_total_sites
    type = ParsedAux
    expression = 'Cu_trapped_sites_1 + Cu_empty_sites_1' # + Cu_trapped_sites_2 + Cu_empty_sites_2'
    coupled_variables = 'Cu_trapped_sites_1 Cu_empty_sites_1' # Cu_trapped_sites_2 Cu_empty_sites_2'
  []
[]

[BCs]
  [left_flux]
    type = EquilibriumBC
    Ko = ${W_S}
    activation_energy = 0
    boundary = left
    enclosure_var = ${pressure} #H2 domain
    temperature = ${T}
    variable = C_M_W
    p = 0.5 #Sievert's Law
  []
  # if including dissociation + recombination on the left
  # [left_dissociation_recombination]
  #   type = ADMatNeumannBC
  #   variable = C_M_W
  #   boundary = 'left'
  #   value = 1
  #   boundary_material = flux_dissociation_recombination_surface
  # []

  [right_recombination]
    type = ADMatNeumannBC
    variable = C_M_Cu
    boundary = 'right'
    value = 1
    boundary_material = flux_recombination_surface
    extra_vector_tags = ref
  []

[]

[InterfaceKernels]
  # Penalized continuity with solubility jump at W/Cu interface
  [W-Cu-interface]
    type = ADPenaltyInterfaceDiffusion
    variable = C_M_W
    neighbor_var = C_M_Cu
    penalty = ${penalty}
    jump_prop_name = solubility_ratio
    boundary = 'interface_W_to_Cu'
    extra_vector_tags = ref
  []
[]

[VectorPostprocessors]
  [W]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '${W_thickness} 0 0'
    num_points = 400
    sort_by = 'x'
    variable = C_M_W
    outputs = vector
  []
  [Cu]
    type = LineValueSampler
    start_point = '${W_thickness} 0 0'
    end_point = '${fparse ${W_thickness} + ${Cu_thickness}} 0 0'
    num_points = 400
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
    point = '${fparse ${W_thickness} + ${Cu_probe}} 0 0'
  []

  [concentration_at_x_W]
    type = PointValue
    variable = C_M_W
    point = '${W_probe} 0 0'
  []

  # check concentration ratio against solubility to ensure penalty is properly enforced
  [Cu_S]
    type = SideAverageMaterialProperty
    property = solubility_Cu_nonAD
    boundary = interface_Cu_to_W
  []

  [W_S]
    type = SideAverageMaterialProperty
    property = solubility_W_nonAD
    boundary = interface_Cu_to_W
  []

  [gold_solubility_ratio]
    type = ParsedPostprocessor
    pp_names = 'W_S Cu_S'
    expression = 'W_S / Cu_S'
  []

  [W_interface]
    type = SideAverageValue
    boundary = interface_W_to_Cu
    variable = C_M_W
  []

  [Cu_interface]
    type = SideAverageValue
    boundary = interface_Cu_to_W
    variable = C_M_Cu
  []

  [variable_ratio]
    type = ParsedPostprocessor
    pp_names = 'W_interface Cu_interface'
    expression = 'W_interface / Cu_interface'
  []

  # inventory (mass Wrrently in the system)
  [mobile_W_mass]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_W
    block = '0'
    outputs = none
  []

  [trapped_W_mass_1]
    type = ElementIntegralVariablePostprocessor
    variable = C_T_W_1
    block = '0'
    outputs = none
  []

  [trapped_W_mass_1_physical]
    type = ScalePostprocessor
    scaling_factor = '${trap_per_free_1}'
    value = trapped_W_mass_1
    outputs = none
  []

  [trapped_W_mass_2]
    type = ElementIntegralVariablePostprocessor
    variable = C_T_W_2
    block = '0'
    outputs = none
  []

  [trapped_W_mass_2_physical]
    type = ScalePostprocessor
    scaling_factor = '${trap_per_free_2}'
    value = trapped_W_mass_2
    outputs = none
  []

  [W_mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'mobile_W_mass + trapped_W_mass_1_physical + trapped_W_mass_2_physical'
    pp_names = 'mobile_W_mass trapped_W_mass_1_physical trapped_W_mass_2_physical'
    outputs = main
  []

  [mobile_Cu_mass]
    type = ElementIntegralVariablePostprocessor
    variable = C_M_Cu
    block = '1'
    outputs = none
  []

  [trapped_Cu_mass_1]
    type = ElementIntegralVariablePostprocessor
    variable = C_T_Cu_1
    block = '1'
    outputs = none
  []

  [trapped_Cu_mass_1_physical]
    type = ScalePostprocessor
    scaling_factor = '${trap_per_free_3}'
    value = trapped_Cu_mass_1
    outputs = none
  []
  

  [Cu_mass_from_concentration]
    type = ParsedPostprocessor
    expression = 'mobile_Cu_mass + trapped_Cu_mass_1_physical' # trapped_Cu_mass_2'
    pp_names = 'mobile_Cu_mass trapped_Cu_mass_1_physical' #trapped_Cu_mass_2'
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
    expression = 'Cu_mass_from_concentration + W_mass_from_concentration'
    pp_names = 'Cu_mass_from_concentration W_mass_from_concentration'
    outputs = main
  []

  # conservation of mass check through residual calculation
  # another way to calculate conservation of mass, didn't quite work/was unable to determine what the
  # reference concentration should be to contextualize the residual
  # [total_inventory_physical]
  #   type = SumPostprocessor
  #   values = 'mobile_W_mass trapped_W_mass_1_physical mobile_Cu_mass'
  #   execute_on = 'INITIAL TIMESTEP_END'
  #   outputs = main
  # []
  # [net_gain_rate]
  #   type = ParsedPostprocessor
  #   pp_names = 'left_influx right_outflux'
  #   expression = '-left_influx - right_outflux'
  #   exeWte_on = 'INITIAL TIMESTEP_END'
  #   outputs = none
  # []

  # [mass_gained_from_flux_physical]
  #   type = TimeIntegratedPostprocessor
  #   value = net_gain_rate
  #   time_integration_scheme = trapezoidal-rule
  #   execute_on = 'INITIAL TIMESTEP_END'
  #   outputs = main
  # []

  # [inventory_change]
  #   type = ChangeOverTimePostprocessor
  #   postprocessor = total_inventory_physical
  #   change_with_respect_to_initial = false
  #   execute_on = 'INITIAL TIMESTEP_END'
  #   outputs = none
  # []

  # # residual: should be ~0 if mass is conserved
  # [mass_conservation_residual]
  #   type = ParsedPostprocessor
  #   pp_names = 'inventory_change mass_gained_from_flux_physical'
  #   expression = 'inventory_change + mass_gained_from_flux_physical'
  #   execute_on = 'INITIAL TIMESTEP_END'
  #   outputs = main
  # []

  # [mass_conservation_residual_relative]
  #   type = ParsedPostprocessor
  #   pp_names = 'mass_conservation_residual total_inventory_physical'
  #   expression = 'if(abs(total_inventory_physical) > 1e-30, mass_conservation_residual / total_inventory_physical, 0)'
  #   execute_on = 'INITIAL TIMESTEP_END'
  #   outputs = main
  # []

  # trapping sites
  # [W_trapped_sites_integral]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = W_trapped_sites_1
  #   block = '0'
  #   outputs = main
  # []

  # [W_total_sites_integral]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = W_total_sites
  #   block = '0'
  #   outputs = main
  # []

  # [W_trap_occupancy_integral]
  #   type = ParsedPostprocessor
  #   expression = 'W_trapped_sites_integral / W_total_sites_integral'
  #   pp_names = 'W_trapped_sites_integral W_total_sites_integral'
  #   outputs = main
  # []

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
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  dtmax = 1
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
  file_base = W_Cu_${RUN_ID}
  [exodus]
    type = Exodus
  []
  [main]
    type = CSV
    file_base = W_Cu_${RUN_ID}_main
  []
  [vector]
    type = CSV
    execute_on = 'FINAL'
  []
[]
