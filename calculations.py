import numpy as np
from ngsolve import *
import matplotlib.pyplot as plt

def Total_energy_loss_cylinder(r_sample, k_bag, eps, T_surface):
    thickness = r_sample #må endres
    D1 = 0.135
    D2 = D1 + 2*thickness
    r1 = D1/2
    r2 = r1 + thickness
    L = 0.3
    A = np.pi*D1*L
    A2 = 2*np.pi*r2*L 

    #temps
    K = 273.15
    T1 = 34 + K
    Ts = T_surface + K # må endres
    T2 = 4.7 + K
    Tf = (Ts+T2)/2

    #General
    g = 9.81
    pa = 1
    sigma = 5.67e-8

    #air
    beta = 1/Tf
    alfa = 1.944e-5
    rho = 1.246
    nu = 1.426e-5
    mu = 1.728e-5
    Ma = 28.97
    rho = 1.246
    cp = 1006
    kappa = 0.02439

    #convection
    Pr = mu*cp/kappa
    Ra =(g*beta*(Ts-T2)*D2**3)/(alfa*nu)
    Nussel = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
    h_cv = (kappa/D2)*Nussel


    h_rad = eps*sigma*(Ts**2 + T2**2)*(Ts+T2)

    R_cv = 1/(h_cv*A2)
    R_cd = np.log(r2/r1)/(2*np.pi*L*k_bag)
    R_rad = 1/(h_rad*A2)

    R_tot = R_cd + R_cv*R_rad/(R_rad+R_cv)

    P = (T1-T2)/R_tot


    print("Area inner cyl: ", A)
    print("Total Area: ", A2)
    print("film_temperature:", Tf)
    print("Prandtl:", Pr)
    print("Raynolds number: ", Ra)
    print("Nussel number: ",Nussel)
    print("h_cv:", h_cv)

    print("R_cd: ",R_cd)
    print("R_cv: ",R_cv)
    print("R_rad: ", R_rad)
    print("Energy: ", P, "W")
    print("P = ", P/A, "W/m^2")


def analytical_cylinder_temp_profile(d, N, T1, T2):
    r1 = 0.0675
    r2 = d+r1
    rs = np.linspace(r1,r2,N)
    return  T1 +(T2-T1)*(np.log(rs)-np.log(r1))/(np.log(r2)-np.log(r1)), rs


def analytical_cylinder_conduction_energy_loss(r, N, kappa, T1, T2): ####HMMMMM very sketch. Don't think this is the way to do it
    L = 0.30
    r1 = 0.0675
    A1 = 2*np.pi*r1*L
    T, rs = analytical_cylinder_temp_profile(r, N, T1, T2)
    q = - 2*np.pi*kappa*L*np.gradient(T)/(A1*rs)
    return q

def cylinder_conduction_flux(T1,T2, d, k):
    L = 0.3
    r1 = 0.0675
    r2 = r1+d
    A = 2*np.pi*r1*L 
    p = 2*np.pi*L*k*(T1-T2)/np.log(r2/r1)
    j = p/A
    return p, j

def get_line_temp_from_gfu_3d(mesh,r_inner,d_sleepingbag, z, gfu, N):
    r_outer = r_inner + d_sleepingbag
    r_vals = np.linspace(r_inner, r_outer, N)

    points = [(r, 0, z) for r in r_vals]
    r_plot = np.linspace(0,d_sleepingbag,N)
    temperatures = [gfu(mesh(x, y, z)) for (x, y, z) in points]
    return temperatures, r_plot


def get_line_temp_from_gfu_2d(mesh,r_inner,d_sleepingbag, gfu, N):
    r_outer = r_inner + d_sleepingbag
    r_vals = np.linspace(r_inner, r_outer, N)

    points = [(r, 0) for r in r_vals]
    r_plot = np.linspace(0,d_sleepingbag,N)
    temperatures = [gfu(mesh(x, y)) for (x, y) in points]
    return temperatures, r_plot

def get_heat_loss(mesh, gfu, h_rad,h_conv, T_ambient):
    W = Integrate((h_rad + h_conv) * (gfu - T_ambient), mesh, BND, definedon=mesh.Boundaries("bag"))
    Area_inner = Integrate(1, mesh, BND, definedon=mesh.Boundaries("baby"))

    W_per_m2 = W / Area_inner

    return W, Area_inner, W_per_m2



def get_h_conv(D, Ts, T2):
    g = 9.81
    K = 273.15
    Ts = Ts + K
    T2 = T2 + K
    Tf = (Ts + T2) / 2
    beta = 1 / Tf
    alfa = 1.944e-5
    nu = 1.426e-5
    mu = 1.728e-5
    cp = 1006
    kappa = 0.02439
    Pr = mu * cp / kappa
    Ra = (g * beta * (Ts - T2) * D**3) / (alfa * nu)
    Nussel = (0.6 + (0.387 * Ra**(1/6)) / (1 + (0.559/Pr)**(9/16))**(8/27))**2
    h_cv = (kappa / D) * Nussel
    return h_cv

def get_h_rad(Ts, T2, eps):
    sigma = 5.67e-8
    K = 273.15
    Ts = Ts + K
    T2 = T2 + K
    h_rad = eps * sigma * (Ts**2 + T2**2) * (Ts + T2)
    return h_rad

def get_h_array(D, T_ambient, eps):
    Temp_array = np.linspace(T_ambient, 100, 1000)
    h_r = get_h_rad(eps,Temp_array, T_ambient)
    h_cv = get_h_conv(D, Temp_array, T_ambient)
    h = h_r + h_cv
    return Temp_array, h

def plot_3d_temperatures(mesh,r_inner,d_sleepingbag, z, gfu, N):
    temp, r = get_line_temp_from_gfu_3d(mesh,r_inner,d_sleepingbag, z, gfu, N)

    plt.plot(r, temp)
    plt.xlabel("Radius (m)")
    plt.ylabel("Temperature (°C)")
    plt.title("Temperature Profile from Inner to Outer Wall")
    plt.grid(True)
    plt.show()

if __name__ == "__main ":
    Total_energy_loss_cylinder(0.0075,0.031,0.8,14)

