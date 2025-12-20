import numpy as np
import sympy as sp
from sympy import var
import matplotlib.pyplot as plt

'''
Please refer to the readme file for the relevant
background and summary of the program.
'''

#Asks the user to input the initial conditions for V(t) and A(t)
t_0_V = input("Enter initial time condition for volume function V(t_0_V) and the time derivative of the volume function V_dot(t_0_V): ")
V_1 = input("Enter initial volume condition V(t_0_V): ")
V_2 = input("Enter initial volume derivative condition V_dot(t_0_V): ")
t_0_A = input("Enter initial time condition for area function A(t_0_A) and the time derivative of the area function A_dot(t_0_A): ")
A_1 = input("Enter initial area condition A(t_0_A): ")
A_2 = input("Enter initial area derivative condition A_dot(t_0_A): ")

#Defines the energy density and momentum of the brane constants, as input from the user
p_1 = float(input("Input the momentum of the brane along the extra dimension 1 (p_1): "))
p_2 = float(input("Input the momentum of the brane along the extra dimension 2 (p_2): "))
T = float(input("Input the energy density of the brane (T): "))

#Defines the momentum of the light ray along one of the brane dimensions and the two extra dimensions, as input from the user
a_1 = float(input("Input the momentum of the light ray along one of the brane dimensions (a_1): "))
a_4 = float(input("Input the momentum of the light ray along the extra dimension 1 (a_4): "))
a_5 = float(input("Input the momentum of the light ray along the extra dimension 2 (a_5): "))

#Sets the initial positions for the brane at time t_ini; Y_1_h_1(t_ini), Y_1_h_2(t_ini), Y_2_h_1(t_ini), and Y_2_h_2(t_ini)
brane_initial_positions = [0, 0, 0, 0]

#Sets the initial positions for the light ray at time t_ini; x_1_h_1(t_ini), x_1_h_2(t_ini), x_4_h_1(t_ini), x_4_h_2(t_ini), x_5_h_1(t_ini), and x_5_h_2(t_ini)
light_ray_initial_positions = [0, 0, 0, 0, 0, 0]

#Allows the user to input the initial time and the final time to specify in the ODE solvers for both the brane and light ray motion
t_ini = float(input("Input initial time: "))
t_end = float(input("Input final time: "))

#Allows the user to input the time step for Euler's and fourth-order Runge-Kutta methods
t_step_euler = float(input("Input the time step for Euler's method: "))
t_step_rk4 = float(input("Input the time step for the fourth-order Runge-Kutta method: "))

'''
This function solves the integration constants p and q for V(t)
and integration constants r and s for A(t) based on the initial 
conditions (t_0_V, t_0_A, V_1, V_2, A_1, and A_2). The function
returns integration constants p, q, r and s.
'''

def init_conditions_volume_area(t_0_V, t_0_A, V_1, V_2, A_1, A_2):
 
    t = var('t')

    p, q, r, s = sp.symbols('p q r s')

    #Creates two Sympy expressions for V(t) and V'(t)
    volume_time = sp.exp(p*t + q)
    volume_prime_time = p * sp.exp(p*t + q)

    #Substitutes t = 0 into V(t) and V'(t) to create equations to solve for p and q
    v_vals = [sp.Eq(volume_time.subs(t, t_0_V), V_1), sp.Eq(volume_prime_time.subs(t, t_0_V), V_2)]
    sol_p_q = sp.solve(v_vals, [p,q], dict=True)[0] 

    #Creates two Sympy expressions for A(t) and A'(t)
    area_time = sp.exp(r*t + s)
    area_prime_time = r * sp.exp(r*t + s)

    #Substitutes t = 0 into A(t) and A'(t) to create equations to solve for r and s
    a_vals = [sp.Eq(area_time.subs(t, t_0_A), A_1), sp.Eq(area_prime_time.subs(t, t_0_A), A_2)]
    sol_r_s = sp.solve(a_vals, [r,s], dict=True)[0] 

    #Enforces condition that p > 0 and that the constraint equation does not yield imaginary solutions
    while sol_p_q[p] < 0 or (2 * ((4/5) * sol_p_q[p] ** 2 - (5/6) * sol_r_s[r] ** 2)) < 0:

        #If the conditions are not met, the user is asked to input the initial volume conditions until the conditions are met
        V_1 = float(input("Please enter a new value for the initial volume condition V(t_0_V): "))
        V_2 = float(input("Please enter a new value for the initial volume derivative condition V_dot(t_0_V): "))

        #Solves for p and q
        v_vals = [sp.Eq(volume_time.subs(t, t_0_V), V_1), sp.Eq(volume_prime_time.subs(t, t_0_V), V_2)]
        sol_p_q = sp.solve(v_vals, [p,q], dict=True)[0] #first solution only

    return float(sol_p_q[p]), float(sol_p_q[q]), float(sol_r_s[r]), float(sol_r_s[s])
    

#Calls the function to solve for the constants of integration for V(t) = e^(pt + q) and A(t) = e^(rt + s) based on the initial conditions
p, q, r, s = init_conditions_volume_area(float(t_0_V), float(t_0_A), float(V_1), float(V_2), float(A_1), float(A_2))

'''
This function defines the volume function V(t) (eq 4 in documentation) and takes in a time t,
integration constant p, and integration constant q. The function returns
the value of the function at time t.
'''

def V_func(t, p, q):
    return np.exp(p*t + q)

'''
This function defines the area function A(t) (eq 5 in documentation) and takes in a time t,
integration constant r, and integration constant s. The function returns
the value of the function at time t.
'''

def A_func(t, r, s):
    return np.exp(r*t + s)

'''
This function defines the constraint equation (eq 9 in documentation) and takes in integration constants p
and r. The function solves the principal value h_1 of the constraint equation, and defines
another constant h_2, which is the negative of h_1. The constraint equation is quadratic, so there
are 2 solutions (h_1 and h_2). The function returns both solutions.
'''

def constraint_eq(p, r):
    
    h_1 = np.sqrt( 2 * ((4/5) * p**2 - (5/6) * r**2))
    h_2 = - h_1
    return h_1, h_2

#Calls the constraint equation function and captures the two solutions 
h_1, h_2 = constraint_eq(p, r)

'''
This function defines the tau_1 function tau_1(t) (eq 7 in documentation) and takes in a time t,
and one of the constants h_1 or h_2. The function returns the special solution
for when alpha = delta = 1 and beta = gamma = 0 for tau_1.
'''

def tau_1_func(h, t):
    return 0

'''
This function defines the tau_2 function tau_2(t) (eq 8 in documentation) and takes in a time t,
and one of the constants h_1 or h_2. The function returns the special solution
for when alpha = delta = 1 and beta = gamma = 0 for tau_2.
'''

def tau_2_func(h, t):
    return np.exp(h * t)

'''
This function defines the ODEs related to the brane motion. One ODE describing
the motion along the extra dimension 1 (Y_1) and one ODE describing
the motion the extra dimension 2 (Y_2) for constant h_1. Another two ODEs 
describing the brane motion for constant h_2. This function takes in the Y_pos, t, 
h_1, h_2, T (energy density of the brane), p_1 (brane momentum 1), and p_2 (brane momentum 2) and
returns the values of the four ODEs at pos Y and time t, as a numpy array. 
'''

def brane_motion_eqs(Y, t, h_1, h_2, T, p_1, p_2):

    #Calls the volume and area functions
    Vt = V_func(t, p, q)
    At = A_func(t, r, s)
    
    #Calls the tau_1 and tau_2 for both h_1 and h_2
    tau_1_h_1 = tau_1_func(h_1, t)
    tau_2_h_1 = tau_2_func(h_1, t)
    tau_1_h_2 = tau_1_func(h_2, t)
    tau_2_h_2 = tau_2_func(h_2, t)

    #Defines the two denominators for h_1 and h_2
    denominator_h_1 = np.sqrt(T**2 * Vt**1.6 * At**(-1) + ((tau_1_h_1**2 + tau_2_h_1**2) / tau_2_h_1) * p_1**2 + (p_2**2 - 2 * tau_1_h_1 * p_1 * p_2) / tau_2_h_1)
    denominator_h_2 = np.sqrt(T**2 * Vt**1.6 * At**(-1) + ((tau_1_h_2**2 + tau_2_h_2**2) / tau_2_h_2) * p_1**2 + (p_2**2 - 2 * tau_1_h_2 * p_1 * p_2) / tau_2_h_2)  

    #Defines the 2 ODEs describing the motion of the brane for h_1 (eqs 2 and 3 in documentation)
    brane_ODE_1_h_1 = ( (Vt**0.8) * (At**(-0.5)) * ((tau_1_h_1**2 + tau_2_h_1**2 - tau_1_h_1)/tau_2_h_1) ) / denominator_h_1
    brane_ODE_2_h_1 = ( (Vt**0.8) * (At**(-0.5)) * (p_2 - tau_1_h_1 * p_1) / (tau_2_h_1) ) / denominator_h_1

    #Defines the 2 ODEs describing the motion of the brane for h_2 (eqs 2 and 3 in documentation)
    brane_ODE_1_h_2 = ( (Vt**0.8) * (At**(-0.5)) * ((tau_1_h_2**2 + tau_2_h_2**2 - tau_1_h_2)/tau_2_h_2) ) / denominator_h_2
    brane_ODE_2_h_2 = ( (Vt**0.8) * (At**(-0.5)) * (p_2 - tau_1_h_2 * p_1) / (tau_2_h_2) ) / denominator_h_2

    return np.array([brane_ODE_1_h_1, brane_ODE_1_h_2, brane_ODE_2_h_1, brane_ODE_2_h_2])

'''
This function defines the ODEs related to the light ray motion. One ODE describing
the motion along one of the brane dimensions (x_1) and two ODEs for the motion along 
the two extra dimensions (x_4 and x_5) for constant h_1. Another three ODEs 
describing the light ray motion for constant h_2. This function takes in the X_pos, t, 
h_1, h_2, a_1 (light ray momentum 1), a_4 (light ray momentum 4), and a_5 (light ray momentum 5) and
returns the values of the six ODEs at pos X and time t, as a numpy array. 
'''

def null_geodesic_motion_eqs(x, t, h_1, h_2, a_1, a_4, a_5):
    
    #Calls the volume and area functions
    Vt = V_func(t, p, q)
    At = A_func(t, r, s)

    #Calls the tau_1 and tau_2 for both h_1 and h_2
    tau_1_h_1 = tau_1_func(h_1, t)
    tau_2_h_1 = tau_2_func(h_1, t)
    tau_1_h_2 = tau_1_func(h_2, t)
    tau_2_h_2 = tau_2_func(h_2, t)

    #Defines the two denominators for h_1 and h_2
    denominator_h_1 = np.sqrt(At ** (2/3) * a_1 ** 2 + (((tau_1_h_1 ** 2 + tau_2_h_1 ** 2) * a_4 ** 2))/(At * tau_2_h_1) + (a_5 ** 2)/(At * tau_2_h_1) - (2 * tau_1_h_1)/(At * tau_2_h_1))
    denominator_h_2 = np.sqrt(At ** (2/3) * a_1 ** 2 + (((tau_1_h_2 ** 2 + tau_2_h_2 ** 2) * a_4 ** 2))/(At * tau_2_h_2) + (a_5 ** 2)/(At * tau_2_h_2) - (2 * tau_1_h_2)/(At * tau_2_h_2)) 

    #Defines the 3 ODEs describing the motion of the light ray for h_1 (eqs 2 and 3 in documentation)
    null_geodesic_eq_1_h_1 = ((Vt ** (4/5)) * (At ** (2/3)) * a_1)/denominator_h_1
    null_geodesic_eq_4_h_1 = ((Vt ** (4/5)) * ((((tau_1_h_1 ** 2 + tau_2_h_1 ** 2) * a_4 - tau_1_h_1 * a_5))/(At * tau_2_h_1)))/denominator_h_1
    null_geodesic_eq_5_h_1 = ((Vt ** (4/5)) * ((-tau_1_h_1 * a_4 + a_5)/(At * tau_2_h_1)))/denominator_h_1

    #Defines the 3 ODEs describing the motion of the light ray for h_2
    null_geodesic_eq_1_h_2 = ((Vt ** (4/5)) * (At ** (2/3)) * a_1)/denominator_h_2
    null_geodesic_eq_4_h_2 = ((Vt ** (4/5)) * ((((tau_1_h_2 ** 2 + tau_2_h_2 ** 2) * a_4 - tau_1_h_2 * a_5))/(At * tau_2_h_2)))/denominator_h_2
    null_geodesic_eq_5_h_2 = ((Vt ** (4/5)) * ((-tau_1_h_2 * a_4 + a_5)/(At * tau_2_h_2)))/denominator_h_2

    return np.array([null_geodesic_eq_1_h_1, null_geodesic_eq_1_h_2, null_geodesic_eq_4_h_1, null_geodesic_eq_4_h_2, null_geodesic_eq_5_h_1, null_geodesic_eq_5_h_2])

'''
This function defines the ODE solver using Euler's method for the ODEs desribing the brane motion. 
This function takes in the Y_pos, t_ini, t_end, time_step, h_1, h_2, T (energy density of the brane),
p_1 (brane momentum 1), and p_2 (brane momentum 2). The function returns the time values and the solution
of the four ODEs.
'''
def brane_evolution_euler(Y, t_ini, t_end, time_step, h_1, h_2, T, p_1, p_2):

    #Defines the number of steps based on the initial time, final time, and time step
    num_steps = int((t_end - t_ini) / time_step) + 1 #The +1 ensures we include the last time value

    #Defines the time array based on the initial time, final time, and number of steps
    t_vals = np.linspace(t_ini, t_end, num_steps)

    #Defines a num_steps by 4 zeros array to store the brane positions for each ODE
    brane_pos = np.zeros((num_steps, 4))

    '''
    brane_pos array is defined as such, where each column corresponds to each ODE

    t_ini                     
    brane_initial_positions[0] 0 0 0 0 0 0 0 0 0 0 0 0 0 ... t_end
    brane_initial_positions[1] 0 0 0 0 0 0 0 0 0 0 0 0 0 ... t_end
    brane_initial_positions[2] 0 0 0 0 0 0 0 0 0 0 0 0 0 ... t_end
    brane_initial_positions[3] 0 0 0 0 0 0 0 0 0 0 0 0 0 ... t_end
    '''
    
    #Sets the initial positions for the brane at time t_ini and includes in the first element of the brane_pos array
    brane_pos[0, :] = brane_initial_positions

    #Defines a function to find the value of each ODE at current time t to find the next step for Euler's method
    def f(Y, t_val):
        return brane_motion_eqs(Y, t_val, h_1, h_2, T, p_1, p_2)
    
    #Calculates the next time step for Euler's method and updates the 4 positions
    for i in range(num_steps - 1):
        
        #Stores time at step i
        t = t_vals[i]

        #Stores the positions of the brane at time step i
        y = brane_pos[i, :]

        #Calculates one step of Euler's method
        k1 = time_step * f(Y, t)
        
        #Updates the four positions
        brane_pos[i+1, :] = y + k1

    return t_vals, brane_pos

'''
This function defines the ODE solver using the fourth-order Runge-Kutta method for the ODEs desribing the brane motion. 
This function takes in the Y_pos, t_ini, t_end, time_step, h_1, h_2, T (energy density of the brane),
p_1 (brane momentum 1), and p_2 (brane momentum 2). The function returns the time values and the solution
of the four ODEs.
'''
def brane_evolution_RK4(Y, t_ini, t_end, time_step, h_1, h_2, T, p_1, p_2):
    
   #Defines the number of steps based on the initial time, final time, and time step
    num_steps = int((t_end - t_ini) / time_step) + 1 #The +1 ensures we include the last time value

    #Defines the time array based on the initial time, final time, and number of steps
    t_vals = np.linspace(t_ini, t_end, num_steps)

    #Defines a num_steps by 4 zeros array to store the brane positions for each ODE
    brane_pos = np.zeros((num_steps, 4))
    
    #Sets the initial positions for the brane at time t_ini and includes in the brane_pos array
    brane_pos[0, :] = brane_initial_positions

    #Defines a function to find the value of each ODE at current time t to find the next step for fourth-order Runge-Kutta method
    def f(Y, t_val):
        return brane_motion_eqs(Y, t_val, h_1, h_2, T, p_1, p_2)

    for i in range(num_steps - 1):
        t = t_vals[i]
        y = brane_pos[i, :] # Current positions (will be updated)
              

        #Calculates one step of fourth-order Runge-Kutta method
        k1 = time_step * f(Y, t)
        k2 = time_step * f(Y + 0.5 * k1, t + 0.5 * time_step)
        k3 = time_step * f(Y + 0.5 * k2, t + 0.5 * time_step)
        k4 = time_step * f(Y + k3, t + time_step)
        
        #Updates the four positions
        brane_pos[i+1, :] = y + (1/6) * (k1 + 2*k2 + 2*k3 + k4)
        
    return t_vals, brane_pos

'''
This function defines the ODE solver using Euler's method for the ODEs desribing the light ray motion. 
This function takes in the x_pos, t_ini, t_end, time_step, h_1, h_2, a_1 (light ray momentum 1), 
a_4 (light ray momentum 4), and a_5 (light ray momentum 5). The function returns the time values and the solution
of the six ODEs.
'''
def null_geodesic_evolution_euler(x, t_ini, t_end, time_step, h_1, h_2, a_1, a_4, a_5):

    #Defines the number of steps based on the initial time, final time, and time step
    num_steps = int((t_end - t_ini) / time_step) + 1 #The +1 ensures we include the last time value

    #Defines the time array based on the initial time, final time, and number of steps
    t_vals = np.linspace(t_ini, t_end, num_steps)

    #Defines a num_steps x 6 zeros array to store the light ray positions for each ODE
    light_ray_pos = np.zeros((num_steps, 6))

    #Sets the initial positions for the light ray at time t_ini and includes in the light_ray_pos array
    light_ray_pos[0, :] = light_ray_initial_positions

    #Defines a function to find the value of each ODE at current time t to find the next step for Euler's method
    def f(x, t_val):
        return null_geodesic_motion_eqs(x, t_val, h_1, h_2, a_1, a_4, a_5)
    
    #Calculates the next time step for Euler's method and updates the 4 positions
    for i in range(num_steps - 1):
        
        #Stores time at step i
        t = t_vals[i]

        #Stores the positions of the light ray at time step i
        pos = light_ray_pos[i, :]

        #Calculates one step of Euler's method
        k1 = time_step * f(x, t)

        #Updates the four positions
        light_ray_pos[i+1, :] = pos + k1

    return t_vals, light_ray_pos

'''
This function defines the ODE solver using the fourth-order Runge-Kutta method for the ODEs desribing the light ray motion. 
This function takes in the x_pos, t_ini, t_end, time_step, h_1, h_2, a_1 (light ray momentum 1), 
a_4 (light ray momentum 4), and a_5 (light ray momentum 5). The function returns the time values and the solution
of the six ODEs.
'''
def null_geodesic_evolution_RK4(x, t_ini, t_end, time_step, h_1, h_2, a_1, a_4, a_5):

    #Defines the number of steps based on the initial time, final time, and time step
    num_steps = int((t_end - t_ini) / time_step) + 1 #The +1 ensures we include the last time value

    #Defines the time array based on the initial time, final time, and number of steps
    t_vals = np.linspace(t_ini, t_end, num_steps)

    #Defines a num_steps by 6 zeros array to store the light ray positions for each ODE
    light_ray_pos = np.zeros((num_steps, 6))

    #Sets the initial positions for the light ray at time t_ini and includes in the light_ray_pos array
    light_ray_pos[0, :] = light_ray_initial_positions

    #Defines a function to find the value of each ODE at current time t to find the next step for Euler's method
    def f(x, t_val):
        return null_geodesic_motion_eqs(x, t_val, h_1, h_2, a_1, a_4, a_5)
    
    for i in range(num_steps - 1):
        
        #Stores time at step i
        t = t_vals[i]

        #Stores the positions of the light ray at time step i
        y = light_ray_pos[i, :] 
              
        #Calculates one step of fourth-order Runge-Kutta method
        k1 = time_step * f(x, t)
        k2 = time_step * f(x + 0.5 * k1, t + 0.5 * time_step)
        k3 = time_step * f(x + 0.5 * k2, t + 0.5 * time_step)
        k4 = time_step * f(x + k3, t + time_step)
        
        #Updates the four positions
        light_ray_pos[i+1, :] = y + (1/6) * (k1 + 2*k2 + 2*k3 + k4)
        
    return t_vals, light_ray_pos
    
'''
Calls both the brane_evolution_euler and brane_evolution_RK4 functions and captures
the time values and the brane positions (the solution). The first argument in the function
can be anything as the ODE has no dependence on the Y value.
'''
t_brane_euler, brane_pos_euler = brane_evolution_euler(0, t_ini, t_end, t_step_euler, h_1, h_2, T, p_1, p_2)
t_brane_RK4, brane_pos_RK4 = brane_evolution_RK4(0, t_ini, t_end, t_step_rk4, h_1, h_2, T, p_1, p_2)

'''
Calls both the null_geodesic_evolution_euler and null_geodesic_evolution_RK4 functions and captures
the time values and the light ray positions (the solution). The first argument in the function
can be anything as the ODE has no dependence on the X value.
'''
t_light_ray_euler, light_ray_pos_euler = null_geodesic_evolution_euler(0, t_ini, t_end, t_step_euler, h_1, h_2, a_1, a_4, a_5)
t_light_ray_RK4, light_ray_pos_RK4 = null_geodesic_evolution_RK4(0, t_ini, t_end, t_step_rk4, h_1, h_2, a_1, a_4, a_5)

#Creates a 2x2 grid of subplots to plot the brane motion plots (total of 4 plots for two ODEs for h_1 and two ODEs h_2)
fig, axes = plt.subplots(2, 2, figsize=(10, 6))

'''
Populates the subplots for each solution to the brane
motion ODEs
'''
axes[0, 0].plot(t_brane_euler, brane_pos_euler[:, 0], 'orange')
axes[0, 0].plot(t_brane_RK4, brane_pos_RK4[:, 0], 'b:', lw = 2)
axes[0, 0].set_xlabel('Time (t)')
axes[0, 0].set_ylabel('Position ($Y^{1}$)')
axes[0, 0].set_title(f'Brane Position along $Y_{1}$ for $h_1 = {h_1:.3f}$')
axes[0, 0].legend(['Euler', 'RK4'], loc = 'best', ncol = 1)

axes[0, 1].plot(t_brane_euler, brane_pos_euler[:, 1], 'gray') 
axes[0, 1].plot(t_brane_RK4, brane_pos_RK4[:, 1], 'b:', lw = 2)
axes[0, 1].set_xlabel('Time (t)')
axes[0, 1].set_ylabel('Position ($Y^{1}$)')
axes[0, 1].set_title(f'Brane Position along $Y_{1}$ for $h_2 = {h_2:.3f}$')
axes[0, 1].legend(['Euler', 'RK4'], loc = 'best', ncol = 1)

axes[1, 0].plot(t_brane_euler, brane_pos_euler[:, 2], 'orange') 
axes[1, 0].plot(t_brane_RK4, brane_pos_RK4[:, 2], 'b:', lw = 2)
axes[1, 0].set_xlabel('Time (t)')
axes[1, 0].set_ylabel('Position ($Y^{2}$)')
axes[1, 0].set_title(f'Brane Position along $Y_{2}$ for $h_1 = {h_1:.3f}$')
axes[1, 0].legend(['Euler', 'RK4'], loc = 'best', ncol = 1)

axes[1, 1].plot(t_brane_euler, brane_pos_euler[:, 3], 'gray') 
axes[1, 1].plot(t_brane_RK4, brane_pos_RK4[:, 3], 'b:', lw = 2)
axes[1, 1].set_xlabel('Time (t)')
axes[1, 1].set_ylabel('Position ($Y^{2}$)')
axes[1, 1].set_title(f'Brane Position along $Y_{2}$ for $h_2 = {h_2:.3f}$')
axes[1, 1].legend(['Euler', 'RK4'], loc = 'best', ncol = 1)

plt.tight_layout()
plt.show()

#Creates a 2x3 grid of subplots to plot the light ray motion plots (total of 6 plots for three ODEs for h_1 and three ODEs h_2)
fig, axes = plt.subplots(3, 2, figsize=(10, 6)) 

axes[0, 0].plot(t_light_ray_euler, light_ray_pos_euler[:, 0], 'orange')
axes[0, 0].plot(t_light_ray_RK4, light_ray_pos_RK4[:, 0], 'b:', lw = 2)
axes[0, 0].set_xlabel('Time (t)')
axes[0, 0].set_ylabel('Position ($x^{1}$)')
axes[0, 0].set_title(f'Light Ray Position along $x^{1}$ for $h_1 = {h_1:.3f}$')
axes[0, 0].legend(['Euler', 'RK4'], loc = 'best', ncol = 1)

axes[0, 1].plot(t_light_ray_euler, light_ray_pos_euler[:, 1], 'gray')
axes[0, 1].plot(t_light_ray_RK4, light_ray_pos_RK4[:, 1], 'b:', lw = 2)
axes[0, 1].set_xlabel('Time (t)')
axes[0, 1].set_ylabel('Position ($x^{1}$)')
axes[0, 1].set_title(f'Light Ray Position along $x^{1}$ for $h_2 = {h_2:.3f}$')
axes[0, 1].legend(['Euler', 'RK4'], loc = 'best', ncol = 1)

axes[1, 0].plot(t_light_ray_euler, light_ray_pos_euler[:, 2], 'orange')
axes[1, 0].plot(t_light_ray_RK4, light_ray_pos_RK4[:, 2], 'b:', lw = 2)
axes[1, 0].set_xlabel('Time (t)')
axes[1, 0].set_ylabel('Position ($x^{4}$)')
axes[1, 0].set_title(f'Light Ray Position along $x^{4}$ for $h_1 = {h_1:.3f}$')
axes[1, 0].legend(['Euler', 'RK4'], loc = 'best', ncol = 1)

axes[1, 1].plot(t_light_ray_euler, light_ray_pos_euler[:, 3], 'gray')
axes[1, 1].plot(t_light_ray_RK4, light_ray_pos_RK4[:, 3], 'b:', lw = 2)
axes[1, 1].set_xlabel('Time (t)')
axes[1, 1].set_ylabel('Position ($x^{4}$)')
axes[1, 1].set_title(f'Light Ray Position along $x^{4}$ for $h_2 = {h_2:.3f}$')
axes[1, 1].legend(['Euler', 'RK4'], loc = 'best', ncol = 1)

axes[2, 0].plot(t_light_ray_euler, light_ray_pos_euler[:, 4], 'orange')
axes[2, 0].plot(t_light_ray_RK4, light_ray_pos_RK4[:, 4], 'b:', lw = 2)
axes[2, 0].set_xlabel('Time (t)')
axes[2, 0].set_ylabel('Position ($x^{5}$)')
axes[2, 0].set_title(f'Light Ray Position along $x^{5}$ for $h_1 = {h_1:.3f}$')
axes[2, 0].legend(['Euler', 'RK4'], loc = 'best', ncol = 1)

axes[2, 1].plot(t_light_ray_euler, light_ray_pos_euler[:, 5], 'gray')
axes[2, 1].plot(t_light_ray_RK4, light_ray_pos_RK4[:, 5], 'b:', lw = 2)
axes[2, 1].set_xlabel('Time (t)')
axes[2, 1].set_ylabel('Position ($x^{5}$)')
axes[2, 1].set_title(f'Light Ray Position along $x^{5}$ for $h_2 = {h_2:.3f}$')
axes[2, 1].legend(['Euler', 'RK4'], loc = 'best', ncol = 1)

plt.tight_layout()
plt.show()

plt.tight_layout()
plt.show()
