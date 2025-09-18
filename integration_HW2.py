import numpy as np
import matplotlib.pyplot as plt
from sympy import var
from sympy import sympify
import argparse

'''
This code uses the sympy package and both the trapezoid and simpson rules to approximate the function E_x for values
of x from a to b in steps whichever the user chooses. The user inputs the function they want to integrate with the
integration bounds and slice size, and the code calculates the integral using both the trapezoid and simpson rules and plots
the results. For the HW task, please run this on the command line: python integration_HW2.py 0 3 30 1000 0.1 and input exp(-1 * t**2)
when prompted. 

-- 09/18/25
'''

#Creates parser
parser = argparse.ArgumentParser()

#Creates arguments for a, b, num_points, N, and increment_step
parser.add_argument('a', type= int, help= 'The lower bound of integration.')
parser.add_argument('b', type= int, help= 'The upper bound of integration.')
parser.add_argument('num_points', type= int, help= 'Number of points for function E_x.')
parser.add_argument('N', type= int, help= 'Integration slice size.')
parser.add_argument('increment_step', type= float, help= 'Step size to determine points within integration bounds a to b.')

#Parses the arguments
args = parser.parse_args()

#Defines function for trapezoid rule
def trapezoid_rule(a, b, num_points, N, increment_step, integration_func):
    
    
    t = var('t')
    
    #Converts user inputted function (a string) into a SymPy object
    expr = sympify(integration_func)

    #Defines a linearly spaced array with x points bounded between a and b
    x_vals = np.linspace(a, b, num_points)

    #Defines an empty array to store each value calculated from the trapezoid rule for a given x value
    E_x_trapezoid = np.array([])
    
    #Loops through each x-val and calculates integral value from a to x-val using trapezoid rule
    for i in x_vals:
        
        #Defines delta x
        delta_x = (i - a)/N

        #Calculates the summation in the trapezoid rule by using the SymPy method subs
        area = 0
        for j in range(N):
            area += expr.subs(t, 0 + j*delta_x) 

        #Calculates the constant terms in the trapezoid rule by using the SymPy method subs
        end_point_1 = 0.5*expr.subs(t, a)
        end_point_2 = 0.5*expr.subs(t, i)

        #Sums up the values from the trapezoid rule and stores it in the empty array defined outside the main loop 
        integral_val_trap_rule = delta_x*(end_point_1 + end_point_2 +area)
        E_x_trapezoid = np.append(E_x_trapezoid, integral_val_trap_rule) 


    return E_x_trapezoid, x_vals

#Defines function for simpson rule
def simpson_rule(a, b, num_points, N, increment_step, integration_func):
    
    t = var('t')

    #Converts user inputted function (a string) into a SymPy object
    expr = sympify(integration_func)

    #Defines a linearly spaced array with x points bounded between a and b
    x_vals = np.linspace(a, b, num_points)

    #Defines an empty array to store each value calculated from the trapezoid rule for a given x value
    E_x_simpson = np.array([])
    
    for i in x_vals:
        
        #Defines delta x
        delta_x = (i - a)/N

        #Calculates the constant terms in the simpson rule by using the SymPy method subs 
        end_point_1 = expr.subs(t, a)
        end_point_2 = expr.subs(t, i)

        #Calculates both summations in the simpson rule by using the SymPy method subs
        area_1 = 0
        area_2 = 0

        for j in range(N//2):
            area_1 += 4*expr.subs(t, (a + (2*j -1) )*delta_x)

        for k in range(N//2 - 1):
            area_2 += 2*expr.subs(t, a+ 2*k*delta_x)

        #Sums up the values from the simpson rule and stores it in the empty array defined outside the main loop 
        integral_val_simp_rule = (delta_x/3)*(end_point_1+end_point_2+area_1+area_2)
        E_x_simpson = np.append(E_x_simpson, integral_val_simp_rule)
    
    return E_x_simpson, x_vals

#Prompts the user to input a function for integration
func = input("Please input the function you want to integrate (use t as your variable). For example, if you want to integrate e^t, input exp(t), or if the function is sint, input sin(t): ")

#Calls both the trapezoid and simpson functions based on the parameters the user inputs when running the .py file on the terminal, and stores the x and E_x values in two arrays
E_x_trapezoid, x_vals = trapezoid_rule(args.a, args.b, args.num_points, args.N, args.increment_step, func)
E_x_simpson, x_vals = simpson_rule(args.a, args.b, args.num_points, args.N, args.increment_step, func)

#Defines two subplots to visualize the function obtained from both the trapezoid and simpson rules
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 4)) 

#Plots the function obtained from the trapezoid rule, and sets the axes labels and title
ax1.plot(x_vals, E_x_trapezoid, c = "b", linestyle = 'dotted')
ax1.set_xlabel("x")
ax1.set_ylabel("E_x")
ax1.set_title("E_x using Trapezoid Rule")

#Plots the function obtained from the simpson rule, and sets the axes labels and title
ax2.plot(x_vals, E_x_simpson, c = "g", linestyle = 'dashdot')
ax2.set_xlabel("x")
ax2.set_ylabel("E_x")
ax2.set_title("E_x using Simpson Rule")

plt.show()