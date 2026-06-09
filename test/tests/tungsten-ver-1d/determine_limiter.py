# determines whether hydrogen diffusion through tungsten is trapping-limited or diffusion-limited
import math

lam = 3.198e-10   #m, lattice constant [1, Table 1] at 1500
nu = 1e13        #s^-1, Debye frequency
rho = 0.01      #-, trapping site fraction [4]
D0 = 7.44e-8    #m^2/s, diffusivity pre-exponential [1]
Ed = 0.13       #eV, diffusion activation energy [1]
ep = 1.2      #eV, trap energy, natural traps are 0.5-0.85, vacancies are 1.2-1.5 [3]
k = 8.617333262e-5    #eV/K Boltzmann's constant
T = 1500        #K, temperature
c = 1e-4        #-, dissolved gas atom fraction. Assume DGA concentration is around 1e-4*lattice density [2]
release_energy = ep/k
print(f'The release energy is: {release_energy} K')
lattice_density = 2/(lam**3) #since BCC, so 2 atoms per unit cell
print(f'The lattice density is: {lattice_density} at/m^3')
num = lam**2*nu*math.exp((Ed-ep)/(k*T))
denom = rho*D0
zeta = num/denom + c/rho

factor = 10 #arbitrary factor to determine whether zeta >> c/rho
threshold = c/rho
print(f"Zeta equals {zeta}")
print(f"c/rho equals {c/rho}")
if zeta > threshold*factor:
    #effective diffusivity limit; diffusion is the rate-limiting step
    print("strongly diffusion-limited")

else:
    #deep-trapping limit; trapping is the rate-limiting step
    print("mixed or trapping-limited")

#Key takeaways: Diffusion is diffusion-limited when the trapping energy is low; 
#It is trapping-limited when it is high. As temperatures rise, it becomes more diffusion-limited.
# [1] Kong, X., et al. First-principles calculations of hydrogen solution and diffusion in tungsten: Temperature and defect-trapping effects. Acta Materialia., v. 84, pp. 426-435. 2015. doi.org/10.1016/j.actamat.2014.10.039
# [2] Ambrosek, J., & Ambrosek, J. Verification and Validation of TMAP7. Idaho National Laboratory. p. 469. 2008.
# [3] Alivaliollahi, et al. “Hydrogen interaction with vacancy defects in tungsten: Unraveling the influence on diffusion mechanisms and mechanical properties.” Int. J. Hyd. En., v. 57, pp. 889-903, 2024.
# [4] Roth, J.; Schmid, K. “Hydrogen in tungsten as plasma-facing material” Phys. Scr., v. T145, pp. 014031, 2011. 

