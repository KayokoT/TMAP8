# Verification Problem #1dd from TMAP7 V&V document
# Permeation Problem without Trapping sites
# No Soret effect or solubility included.

# modeling parameters
nx_num = 200 # (-)
diffusivity = '${units 2.72e-8 m^2/s}' #at 1500 K
simulation_time = ${units 0.01 s}
interval_time_min = ${units 1e-5 s}
interval_time = ${units 1e-5 s}
cl = '${units 6.115e24 at/m^3}' 

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = ${nx_num}
  xmax = '${units 2.5e-5 m}'#-> mum}
[]

[Variables]
  [mobile]
  []
[]

[Kernels]
  [diff]
    type = MatDiffusion
    variable = mobile
    diffusivity = ${diffusivity}
  []
  [time]
    type = TimeDerivative
    variable = mobile
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = mobile
    value = '${fparse cl / cl}'
    boundary = left
  []
  [right]
    type = DirichletBC
    variable = mobile
    value = 0
    boundary = right
  []
[]

[Postprocessors]
  [outflux]
    type = SideDiffusiveFluxAverage
    boundary = 'right'
    diffusivity = ${diffusivity}
    variable = mobile
  []
  [scaled_outflux]
    type = ScalePostprocessor
    value = outflux
    scaling_factor = ${cl} 
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
  end_time = ${simulation_time}
  dt = ${interval_time}
  dtmin = ${interval_time_min}
  solve_type = NEWTON
  scheme = BDF2
  nl_abs_tol = 1e-13
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  verbose = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  csv = true
  [dof]
    type = DOFMap
    execute_on = initial
  []
  perf_graph = true
[]
