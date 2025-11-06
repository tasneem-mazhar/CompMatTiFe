#!/usr/bin/env python3
"""
Simple Birch-Murnaghan EOS fitting for TiFe
Usage: python fit_eos.py ev_data.txt
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys

def birch_murnaghan(V, E0, V0, B0, B0_prime):
    """Birch-Murnaghan equation of state"""
    # Conversion: 1 Ry/Ang^3 = 14710.513 / (0.529177)^3 GPa
    conversion = 14710.513 / (0.529177**3)
    eta = (V0 / V)**(2.0/3.0)
    E = E0 + (9.0 * V0 * B0 / 16.0 / conversion) * (
        (eta - 1.0)**2 * (6.0 - 4.0 * eta + B0_prime * (eta - 1.0))
    )
    return E

def pressure_bm(V, V0, B0, B0_prime):
    """Pressure from BM EOS"""
    eta = (V0 / V)**(1.0/3.0)
    P = (3.0 * B0 / 2.0) * (eta**7 - eta**5) * (
        1.0 + 0.75 * (B0_prime - 4.0) * (eta**2 - 1.0)
    )
    return P

def solve_volume_at_pressure(V0, B0, B0_prime, P_target):
    """Find volume at given pressure"""
    from scipy.optimize import fsolve
    def equation(V):
        return pressure_bm(V, V0, B0, B0_prime) - P_target
    V_guess = V0 * (1.0 - P_target / B0)
    return fsolve(equation, V_guess)[0]

# Read data
if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = 'ev_data.txt'

print("Reading:", filename)
data = np.loadtxt(filename)
volumes = data[:, 0]
energies = data[:, 1]

print("\nData points:", len(volumes))
for i, (v, e) in enumerate(zip(volumes, energies), 1):
    print(f"  {i}. V={v:.4f} Ang^3, E={e:.8f} Ry")

# Initial guess
E0_guess = np.min(energies)
V0_guess = volumes[np.argmin(energies)]
B0_guess = 150.0
B0_prime_guess = 4.0

# Fit
popt, pcov = curve_fit(birch_murnaghan, volumes, energies,
                       p0=[E0_guess, V0_guess, B0_guess, B0_prime_guess])
perr = np.sqrt(np.diag(pcov))

E0, V0, B0, B0_prime = popt

# Calculate residuals and quality metrics
E_fit = birch_murnaghan(volumes, *popt)
residuals = energies - E_fit
rmse = np.sqrt(np.mean(residuals**2))
max_residual = np.max(np.abs(residuals))

# Calculate R-squared
ss_res = np.sum(residuals**2)
ss_tot = np.sum((energies - np.mean(energies))**2)
r_squared = 1 - (ss_res / ss_tot)

# Convert bulk modulus to other units
B0_kbar = B0 * 10.0  # GPa to kbar
B0_eV_Ang3 = B0 * 0.00624150913  # GPa to eV/Ang^3

# Calculate lattice parameter at equilibrium
a0 = V0**(1.0/3.0)  # for cubic

# Print results
print("\n" + "="*70)
print("BIRCH-MURNAGHAN EQUATION OF STATE FIT RESULTS")
print("="*70)
print("\nFitted Parameters:")
print("-" * 70)
print(f"E0 (equilibrium energy):     {E0:.10f} +/- {perr[0]:.10f} Ry")
print(f"V0 (equilibrium volume):     {V0:.6f} +/- {perr[1]:.6f} Ang^3")
print(f"a0 (lattice parameter):      {a0:.6f} Ang (for cubic cell)")
print(f"B0 (bulk modulus):           {B0:.2f} +/- {perr[2]:.2f} GPa")
print(f"                             {B0_kbar:.1f} kbar")
print(f"                             {B0_eV_Ang3:.6f} eV/Ang^3")
print(f"B0' (pressure derivative):   {B0_prime:.3f} +/- {perr[3]:.3f}")
print("\nFit Quality:")
print("-" * 70)
print(f"R-squared:                   {r_squared:.8f}")
print(f"RMSE:                        {rmse:.8f} Ry ({rmse*1000:.4f} mRy)")
print(f"Max residual:                {max_residual:.8f} Ry ({max_residual*1000:.4f} mRy)")
print("="*70)

# Volume at different pressures
print("\nVolume at Different Pressures:")
print("-" * 70)
print("  Pressure         Volume          Change")
print("  (GPa)           (Ang^3)          (%)")
print("-" * 70)
for P in [0, 5, 10, 20, 50, 100]:
    V_P = solve_volume_at_pressure(V0, B0, B0_prime, P)
    dV_p = ((V_P - V0) / V0) * 100
    print(f"  {P:3d}            {V_P:8.4f}         {dV_p:+6.2f}")
print("-" * 70)

# Highlight 10 GPa result
V_10GPa = solve_volume_at_pressure(V0, B0, B0_prime, 10.0)
dV = V_10GPa - V0
dV_percent = (dV / V0) * 100

print(f"\nAt P = 10 GPa (100 kbar):")
print(f"  Volume:                      {V_10GPa:.6f} Ang^3")
print(f"  Change:                      {dV:.6f} Ang^3 ({dV_percent:.2f}%)")
print(f"  Compression:                 {-dV_percent:.2f}%")
print("="*70)

# Plot
V_fit = np.linspace(min(volumes)*0.95, max(volumes)*1.05, 200)
E_fit = birch_murnaghan(V_fit, *popt)

plt.figure(figsize=(10, 6))
plt.plot(volumes, energies, 'ro', markersize=10, label='QE data')
plt.plot(V_fit, E_fit, 'b-', linewidth=2, label='BM EOS fit')
plt.axvline(V0, color='g', linestyle='--', alpha=0.5, 
           label=f'$V_0$ = {V0:.3f} $\AA^3$')
plt.xlabel('Volume ($\AA^3$)', fontsize=14)
plt.ylabel('Energy (Ry)', fontsize=14)
plt.title('TiFe: E(V) Curve and Birch-Murnaghan Fit', fontsize=14, fontweight='bold')
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)

# Add text box with results
textstr = f'$V_0$ = {V0:.3f} $\AA^3$\n$B_0$ = {B0:.1f} GPa\n$B_0\'$ = {B0_prime:.2f}\n$R^2$ = {r_squared:.6f}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, 
        fontsize=11, verticalalignment='top', bbox=props)

plt.tight_layout()
plt.savefig('eos_fit.png', dpi=300)
print("\nPlot saved: eos_fit.png")

# Save results
with open('eos_results.txt', 'w') as f:
    f.write("="*70 + "\n")
    f.write("BIRCH-MURNAGHAN EQUATION OF STATE FIT RESULTS\n")
    f.write("="*70 + "\n\n")
    f.write("Fitted Parameters:\n")
    f.write("-"*70 + "\n")
    f.write(f"E0 (equilibrium energy):     {E0:.10f} +/- {perr[0]:.10f} Ry\n")
    f.write(f"V0 (equilibrium volume):     {V0:.6f} +/- {perr[1]:.6f} Ang^3\n")
    f.write(f"a0 (lattice parameter):      {a0:.6f} Ang (for cubic cell)\n")
    f.write(f"B0 (bulk modulus):           {B0:.2f} +/- {perr[2]:.2f} GPa\n")
    f.write(f"                             {B0_kbar:.1f} kbar\n")
    f.write(f"B0' (pressure derivative):   {B0_prime:.3f} +/- {perr[3]:.3f}\n")
    f.write(f"\nFit Quality:\n")
    f.write("-"*70 + "\n")
    f.write(f"R-squared:                   {r_squared:.8f}\n")
    f.write(f"RMSE:                        {rmse:.8f} Ry\n")
    f.write(f"Max residual:                {max_residual:.8f} Ry\n")
    f.write(f"\nAt P = 10 GPa (100 kbar):\n")
    f.write(f"  Volume:                    {V_10GPa:.6f} Ang^3\n")
    f.write(f"  Change:                    {dV:.6f} Ang^3 ({dV_percent:.2f}%)\n")
    f.write(f"  Compression:               {-dV_percent:.2f}%\n")

print("Results saved: eos_results.txt")
