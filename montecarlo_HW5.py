import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
import numpy as np
import argparse

'''
This code does three simulations, including a photon/multiple photons traveling in
a slab and a photon traveling in the Sun. The code produces an animation for a
single photon in a slab and for the photon in the Sun. For the parameters 
specified in the HW, run the following command in the terminal:
python montecarlo_HW5.py 1000 150 0 0 0 100

-- 10/26/25
'''

#Creates parser
parser = argparse.ArgumentParser()

#Creates arguments for a, b, num_points, N, and increment_step
parser.add_argument('slab_width', type= int, help= 'The slab width the user wants to consider in m.')
parser.add_argument('l', type= int, help= 'The mean free path in m.')
parser.add_argument('init_pos_x', type= int, help= 'The initial x position of the photon in m.')
parser.add_argument('init_pos_y', type= int, help= 'The initial y position of the photon in m.')
parser.add_argument('init_pos_r', type= int, help= 'The initial r position of the photon in the Sun in cm.')
parser.add_argument('num_photons', type= int, help= 'The initial r position of the photon in the Sun in cm.')

#Parses the arguments
args = parser.parse_args()

#Defines a function to simulate a single photon in a slab w/ electrons
def one_photon_slab(slab_width, init_pos_x, init_pos_y, l):

    #Defines initial x and y positions of the photon
    init_x = init_pos_x
    init_y = init_pos_y

    #Defines total dist in m
    tot_dist = 0

    #Stores x and y postions in an array
    x_pos = [init_x]
    y_pos = [init_y]

    #Defines counter, a Boolean to track if the photon escaped, and the maximum number of scatters
    i = 0
    photon_alive = True
    max_scatters = 10000

    #Loops while the photon is still in the slab and the counter is less than the number of scatters
    while (photon_alive and i < max_scatters): 

     #Defines scattering distance from an exponential distribution; the distance the photon travels after a scattering event
     distance = -l * np.log(np.random.rand())

     #Updates the total distance   
     tot_dist += distance

     #Defines the angle the photon scatters from -pi to pi from a uniform distribution
     theta = np.pi * np.random.uniform(-1, 1)

     #Updates the x and y positions, respectively based on the scattering angle    
     new_x = x_pos[i] + distance * np.cos(theta)
     new_y = y_pos[i] + distance * np.sin(theta)

     #Checks if the photon is within the slab and appends new position to x_pos and y_pos arrays     
     if (init_pos_x <= new_x <= slab_width and init_pos_y <= new_y <= slab_width):
        x_pos.append(new_x)
        y_pos.append(new_y)
        i += 1

    #If it is not, the simulation stops as the photon has escaped the slab
    else:
        photon_alive = False

    #Checks to see if the photon got reflected back to the original position, if so, the simulation restarts
    if new_x == init_x and new_y == init_y:
        x_pos = [init_x]
        y_pos = [init_y]

        new_x = init_x
        new_y = init_y
        
        tot_dist = 0

        i = 0

    return x_pos, y_pos, tot_dist 

#Defines a function to simulate multiple photons in a slab w/ electrons
def multiple_photons_slab(num_photons, slab_width, init_pos_x, init_pos_y, l):

    num_scatter = [] #Number of scatterings per photon
    escaped = 0 #Counter for photons that transmit through slab
    reflected = 0 #Counter for photons that reflect back

    
    #for loop to iterate through all photons
    for p in range(num_photons):

        #Defines initial x and y positions of the photon
        init_x = init_pos_x
        init_y = init_pos_y

        #Defines total dist in m
        tot_dist = 0

        #Stores x and y postions in an array
        x_pos = [init_x]
        y_pos = [init_y]

        #Defines counter, a Boolean to track if the photon escaped, and the maximum number of scatters
        i = 0
        photon_alive = True
        max_scatters = 10000

        #Loops while the photon is still in the slab and the counter is less than the number of scatters
        while (i < max_scatters and photon_alive): 

            #Checks if the last position is greater than the slab width, if so, goes to the next photon
            if x_pos[-1] > slab_width:
                escaped += 1
                photon_alive = False
                break

            #Checks if the last position is less than the slab width, if so, goes to the next photon        
            if x_pos[-1] < init_pos_x:
                reflected += 1
                photon_alive = False
                break

            #Defines scattering distance from an exponential distribution; the distance the photon travels after a scattering event
            distance = -l * np.log(np.random.rand())

            #Updates the total distance     
            tot_dist += distance

            #Defines the angle the photon scatters from -pi to pi from a uniform distribution
            theta = np.pi * np.random.uniform(-1, 1)

           #Updates the x and y positions, respectively based on the scattering angle    
            new_x = x_pos[i] + distance * np.cos(theta)
            new_y = y_pos[i] + distance * np.sin(theta)

            #Checks if the photon is within the slab and appends new position to x_pos and y_pos arrays  
            if (init_pos_x <= new_x <= slab_width and init_pos_y <= new_y <= slab_width):
                x_pos.append(new_x)
                y_pos.append(new_y)
                i += 1
            
            #If the photon is not within the slab, it must have escaped or reflected
            else:

                #Checks whether the new x position is greater than the slab width, if so, the photon has escaped
                if new_x > slab_width:
                    escaped += 1
                    photon_alive = False
                    
                #Checks whether the new x position is less than the slab width, if so, the photon has reflected
                elif new_x < init_pos_x:
                    reflected += 1
                    photon_alive = False
                
            #Adds the number of scatters for specific photon to num scatter array
            num_scatter.append(i)

    #Calculates the probability that a photon escaped or reflected
    prob_escaped = escaped/num_photons
    prob_refleted = reflected/num_photons

    #Displays probabilities
    print(f'The probability that a photon escapes is {prob_escaped*100}% for a smaple of {num_photons}.')
    print(f'The probability that a photon reflects is {prob_refleted*100}% for a smaple of {num_photons}.')

#Defines a function to simulate a photon in the Sun
def one_photon_sun(init_pos_r):

    #Defines the radius of the Sun and Thompson scattering in units of cm
    R_sun = 696340 * 10**5
    sigma_T = 6.652 * 10**-25

    #Defines an array with initial radius
    r = [init_pos_r]

    #Defines the electron number density in cm^-3
    n_e = (2.5 * 10**26) * np.exp(-r[-1]/(0.096 * R_sun)) 

    #Defines the mean free path 
    l = (n_e * sigma_T) ** -1

    #Initializes the total distance traveled
    tot_dist = 0

    #Defines counter and the maximum number of scatters
    i = 0
    max_scatters = 10000

    #Loops while counter is less than the number of max scatters
    while (i < max_scatters):

        #Defines scattering distance from an exponential distribution; the distance the photon travels after a scattering event
        distance = -l * np.log(np.random.rand())

        #Updates the total distance  
        tot_dist += distance

        #Defines the angle the photon scatters from 0 to pi from a uniform distribution
        theta = np.random.uniform(0, np.pi)

        #Calculates the change in r
        del_r = distance * np.cos(theta)
        

        #Checks if the change in r is negative
        if del_r < 0:
            del_r = np.abs(del_r)

        #Calculates new position r    
        new_R = r[-1] + distance * np.cos(theta)
        
        r.append(new_R)

        i += 1
    
    return r, tot_dist

#Calls each function
pos_x, pos_y, dist = one_photon_slab(args.slab_width, args.init_pos_x, args.init_pos_y, args.l)
multiple_photons_slab(args.num_photons, args.slab_width, args.init_pos_x, args.init_pos_y, args.l)
r_pos, r_tot_dist = one_photon_sun(args.init_pos_r)

#Creates a figure for the one photon slab simulation, used claude.AI to debug section
fig, ax = plt.subplots(figsize=(10, 6))

#Initializes slab by starting the rectangle at coordinate (0,0)
slab_rect = Rectangle((0, 0), args.slab_width, args.slab_width, facecolor='lightblue', alpha=0.2, edgecolor='blue', linewidth=2)

ax.add_patch(slab_rect)

# Initializes the line to track the photon path and point to track the position of the photon
line, = ax.plot([], [], 'g-', linewidth=1, alpha=0.7, label='Photon path')
point, = ax.plot([], [], 'ro', markersize=8, label='Current position')

# Defines axis properties
ax.set_xlim(min(pos_x), max(pos_x))
ax.set_ylim(min(pos_y), max(pos_y))
ax.set_xlabel('x position (m)')
ax.set_ylabel('y position (m)')
ax.set_title(f'Photon Random Walk in Slab (Mean Free Path = {args.l} m)')
ax.legend()
ax.grid(True, alpha=0.3)

# Generates textbox at the top left of the slab to keep track of the scatter count 
scatter_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, verticalalignment='top', fontsize=10, bbox=dict(boxstyle='round', facecolor='salmon', alpha=0.5))

# Initializes function
def init():
    line.set_data([], [])
    point.set_data([], [])
    scatter_text.set_text('')
    return line, point, scatter_text

# Animation function
def animate(i):
    
    x = pos_x[:i+1] #Saves x points, so for i = 0 array (hypothetical points) the array would be [0], and for i = 1 the array would be [0, 1], saves previous point in array
    y = pos_y[:i+1] #Same logic as previous line
    line.set_data(x, y)
    
    # Shows the current position and updates textbox to show the scatter number
    if i < len(pos_x):
        point.set_data([pos_x[i]], [pos_y[i]])
        scatter_text.set_text(f'Scatters: {i}/{len(pos_x)-1}')
    
    return line, point, scatter_text

# Creates the animation
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=len(pos_x), interval=50, blit=True, repeat=False)

plt.show()

#Creates a figure for the one photon Sun simulation, used claude.AI to debug section
fig, ax = plt.subplots(figsize=(10, 6))

#Defines speed of light in m/s
c = 3 * 10**8

#Finds the total time (in s) the photon traveled for the simulation
total_time = (r_tot_dist/100) / c

#Defines a time array
t = np.linspace(0, total_time, len(r_pos))

#Converts each element in the r_pos array to meters
for i in range(len(r_pos)):
    r_pos[i] = r_pos[i]/100

# Initializes the line to track the photon path and point to track the position of the photon
line, = ax.plot([], [], 'g-', linewidth=1, alpha=0.7, label='Photon path')
point, = ax.plot([], [], 'ro', markersize=8, label='Current position')

# Defines axis properties
ax.set_ylim(min(r_pos), max(r_pos))
ax.set_xlim(0, total_time)
ax.set_xlabel('Time t (s)')
ax.set_ylabel('Radial Position (m)')
ax.set_title(f'Photon Random Walk in the Sun (Starting from {r_pos[0]} m)')
ax.legend()
ax.grid(True, alpha=0.3)

# Generates textbox at the top left of the figure to keep track of the scatter count 
scatter_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, verticalalignment='top', fontsize=10, bbox=dict(boxstyle='round', facecolor='salmon', alpha=0.5))

# Initializes function
def init():
    line.set_data([], [])
    point.set_data([], [])
    scatter_text.set_text('')
    return line, point, scatter_text

# Animation function
def animate(i):
    
    time = t[:i+1] #Saves time t points, so for i = 0 array (hypothetical points) the array would be [0], and for i = 1 the array would be [0, 1], saves previous point in array
    radial = r_pos[:i+1] #Same logic as previous line
    line.set_data(time, radial)
    
    # Shows the current radial position and updates textbox to show the scatter number
    if i < len(t):
        point.set_data([t[i]], [r_pos[i]])
        scatter_text.set_text(f'Scatters: {i}/{len(r_pos)-1}')
    
    return line, point, scatter_text

# Creates the animation
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=len(r_pos), interval=50, blit=True, repeat=False)

plt.show()
