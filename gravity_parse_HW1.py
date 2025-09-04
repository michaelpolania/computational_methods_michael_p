import argparse

'''
This code uses argparse to calculate the time it takes for an object to reach the ground, given the initial height and gravitational acceleration.
The user will input the height (h) and gravitational acceleration (g) in the terminal when running the .py file. For instance, if h = 1 and g = 9.8,
then the user will run python gravity_parse_HW1.py 1 9.8. It is important to note that h must go first then g (order matters).

-- 09/03/25
'''

#Creates parser
parser = argparse.ArgumentParser()

#Creates arguments for the height and gravitational acceleration
parser.add_argument('h', type= float, help= 'The height of the object above ground in meters.')
parser.add_argument('g', type= float, help= 'The gravitational acceleration constant in m/s^2.')

#Parses the arguments
args = parser.parse_args()

#Defines a function to calculate the time for an object to reach the ground in a gravational field, given height above the ground
def calc_time (h, g):
    t = ((2*h)/g)**.5 #Equation to solve for time
    print(f'The ball takes {t:.2f} seconds to reach the ground.') #Prints the time up to two sig figs
    return t

#Calls the function based on the height and gravitational acceleration the user inputs when running the .py file on the terminal
time(args.h, args.g)
