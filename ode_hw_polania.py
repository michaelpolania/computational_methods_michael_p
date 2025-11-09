import numpy as np
import argparse
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots   

'''
This code simulates the orbit of a ball orbiting around a rod with mass M and length L 
using both the fourth-order Runge-Kutta method and the Leapfrog method. In addition, 
the code uses argparse, so the user can input the initial conditions, mass, length, 
and time step on the terminal. The code uses plotly to plot the orbits using both methods. 
Based on both plots, we note that both methods yield similar results for the orbit of the ball.
For the assignment, please run this command on the terminal: python ode_hw_polania.py 1 0 0 1 10 2 0.0001 

-- 11/08/25
'''
#Creates parser
parser = argparse.ArgumentParser()

#Creates arguments for the initial conditions, mass, length, and time step
parser.add_argument('x0', type= float, help= 'The initial x position.')
parser.add_argument('y0', type= float, help= 'The initial y position.')
parser.add_argument('vx0', type= float, help= 'The initial x velocity.')
parser.add_argument('vy0', type= float, help= 'The initial y velocity.')
parser.add_argument('M', type= float, help= 'The mass of the rod.')
parser.add_argument('L', type= float, help= 'The rod length.')
parser.add_argument('h', type= float, help= 'The time step.')

#Parses the arguments
args = parser.parse_args()

#Defines initial conditions
x0 = args.x0 
y0 = args.y0
vx0 = args.vx0
vy0 = args.vy0

#Vectorizes initial conditions
r0 = np.array([x0, y0, vx0, vy0])

#Defines constants for problem
G = 1
M = args.M
L = args.L

#Defines a function for the four linear equations of motion
def fgrav(r, t):

    x_dot = r[2]
    y_dot = r[3]

    vx_dot = -G*M*r[0]/((r[0] ** 2 + r[1] ** 2) *np.sqrt(r[0] ** 2 + r[1]** 2 + 0.25*L**2))
    vy_dot = -G*M*r[1]/((r[0] ** 2 + r[1] ** 2) *np.sqrt(r[0] ** 2 + r[1]** 2 + 0.25*L**2))

    return np.array([x_dot, y_dot, vx_dot, vy_dot])

#Defines a function to calculate one step of a fourth-order Runge-Kutta method
def rk4_step(r, t, h, f):
    k1 = h * f(r, t)
    k2 = h * f(r + 0.5 * k1, t + 0.5 * h)
    k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
    k4 = h * f(r + k3, t + h)
    return r + (k1 + 2 * k2 + 2 * k3 + k4) / 6

#Defines a function to calculate one step of the leapfrog method
def leapfrog_step(r, t, h, f):
    k1 = r + .5 * h * f(r, t)
    k2 = r + h * f(k1, t + .5 * h)
    k3 = k1 + h * f(k2, t + h)
    return k2 + h * f(k3, t + 1.5 * h)

#Defines a function to evolve the orbit over time using the RK4 method
def evolve_orbit_rk4(h):

    t = np.arange(0, 10, h)
    x = np.zeros(len(t))
    y = np.zeros(len(t))
    r = r0

    for time in range(len(t)):
        x[time] = r[0]
        y[time] = r[1]
        r = rk4_step(r, time, h, fgrav)
    return x, y, t

#Defines a function to evolve the orbit over time using the leapfrog method
def evolve_orbit_leapfrog(h):

    t = np.arange(0, 10, h)
    x = np.zeros(len(t))
    y = np.zeros(len(t))
    r = r0

    for time in range(len(t)):
        x[time] = r[0]
        y[time] = r[1]
        r = leapfrog_step(r, time, h, fgrav)
    return x, y, t

#Calls the functions to evolve the orbit using both the RK4 and leapfrog methods
x_pos_rk4, y_pos_rk4, t  = evolve_orbit_rk4(args.h)
x_pos_leapfrog, y_pos_leapfrog, t  = evolve_orbit_leapfrog(args.h)

#Create a subplot using Plotly with 1 row and 2 columns
fig = make_subplots(rows=1, cols=2, subplot_titles=('Orbit using RK4', 'Orbit using Leapfrog'))

# Adds the traces for both plots to the subplot figure
fig.add_trace(
    go.Scatter(x=x_pos_rk4, y=y_pos_rk4, mode='lines', name='RK4 Method'),
    row=1, col=1
)

fig.add_trace(
    go.Scatter(x=x_pos_leapfrog, y=y_pos_leapfrog, mode='lines', name='Leapfrog Method'),
    row=1, col=2
)

#Creates axis titles for both plots
fig.update_xaxes(title_text="x position", row=1, col=1)
fig.update_yaxes(title_text="y position", row=1, col=1)
fig.update_xaxes(title_text="x position", row=1, col=2)
fig.update_yaxes(title_text="y position", row=1, col=2)

fig.show()