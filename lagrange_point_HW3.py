import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const

'''
This code solves for the root of the quintic equation given in the homework
directions by using both the Newton and Secant methods, respectively. 
The program then plots the points obtained from both the methods
and displays to the user.

-- 09/24/25
'''

#Defines constants and parameters used in the equation

m = 7.348e22 #Mass of the moon in kg
M = const.M_earth.value #Mass of the earth in kg

omega = 2.662e-6 #Angular velocity of both the Moon and the satellite in Hz
R = 3.844e8 #Earth-moon distance in m
G = const.G.value #Newton's Universal Law of Gravitation constant

#Defines a function to represent the equation we are trying to solve
def func(x):
    return (G*M)/(x**2) - (G*m)/((R-x)**2) - x * omega**2

#Defines a function to represent the derivative of the equation we are solving
def deriv(x):
    return (-2*G*M)/(x**3) - (2*G*m)/((R-x)**3) - omega**2

#Defines a function to calculate the slope of the secant line
def secant(x_1, x_2):
    return (func(x_2) - func(x_1))/(x_2 - x_1)


#Plots a certain range of r values to estimate the root of the equation to use the Newton and/or Secant method

r = np.linspace(2e8, 3.4e8, 1000) #3.262e8 looks like the zero
plt.plot(r, func(r))
plt.xlabel('r: distance from the center of Earth to L1 lagrange point (m)')
plt.ylabel('f(r): equation we are solving')
plt.title('f(r) vs. r; graph to estimate root(s)')
plt.show()

#Defines array with the initial guess of the solution for the Newton method; obtained from previous graph
r_guess = [3.262e8]

#For loop for Newton's method

for i in range(100):
    
    #Defines a variable to calculate the next point for the Newton method
    val = r_guess[i] - func(r_guess[i])/deriv(r_guess[i])
    
    #Checks if the difference between the new value and previous value are below a tolerance value, if they are, then they are roughly equal and we have found a solution
    if abs(val - r_guess[i]) < 1.e-10:
        break

    else:

        #Appends new point to r_guess array
        r_guess.append(val)


#Displays root in meters and kilometers
print(f'From Newtons method, the solution or root of the equation is r = {r_guess[-1]:.4f} m = {r_guess[-1]/1000:.4f} km. ')       

#Defines array with the initial two guesses of the solution for the Secant method
r_guess_1 = [3.262e8, 3.35e8]


#For loop for Secant method        
for i in range(100):

    #Defines a variable to calculate the next point for the Secant method
    val = r_guess_1[i + 1] - func(r_guess_1[i + 1])/secant(r_guess_1[i], r_guess_1[i + 1])
    

    #Checks if the difference between the new value and previous value are below a tolerance value, if they are, then they are roughly equal and we have found a solution
    if abs(val - r_guess_1[i + 1]) < 1.e-10:
        break

    else:

        #Appends new point to r_guess_2 array
        r_guess_1.append(val)


#Displays root in meters and kilometers           
print(f'From the Secant method, the solution or root of the equation is r = {r_guess_1[-1]:.4f} m = {r_guess_1[-1]/1000:.4f} km. ')        

#Defines two subplots to visualize the points obtained from both the Newton and Secant methods
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 4)) 

#Redefines the arrays with points to numpy arrays for plotting
r_guess = np.array(r_guess)
r_guess_1 = np.array(r_guess_1)


#Plots the points obtained from the Newton method against the function value at specific point
ax1.plot(r_guess, func(r_guess), c = "b", linestyle = 'dotted')
ax1.set_xlabel("r")
ax1.set_ylabel("f(r)")
ax1.set_title("Plot showing f(r) vs. r for Newton method ")

#Plots the root of the equation and creates legend
ax1.plot(r_guess[-1], func(r_guess[-1]), 'r*', markersize=10, label="Root (Newton)")
ax1.legend()

#Plots the points obtained from the Secant method against the function value at specific point
ax2.plot(r_guess_1, func(r_guess_1), c = "g", linestyle = 'dashdot')
ax2.set_xlabel("r")
ax2.set_ylabel("f(r)")
ax2.set_title("Plot showing f(r) vs. r for Secant method")   

#Plots the root of the equation and creates legend
ax2.plot(r_guess_1[-1], func(r_guess_1[-1]), 'b*', markersize=10, label = "Root (Secant)")
ax2.legend()


plt.show()

