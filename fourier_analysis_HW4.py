from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import argparse
from numpy import zeros
from cmath import exp,pi

'''
This code reads in the TIC0010756751.fits file and computes the Fourier transform using two approaches.
In the first approach, we compute the full Discrete Fourier Transform, and create a power spectrum. 
We then compute the inverse Discrete Fourier Transform and compare to our original light curve. 
In the second aproach, we use a linear interpolation and use the numpy fft method to compute
the fast Fourier transform to create a power spectrum. We then use the numpy ifft method to compute the
inverse fast Fourier transform and compare to original light curve. We note that both methods yield similar
plots. 

Please run both of these commands in the terminal to analyze both time epochs studied: 
python fourier_analysis_HW4.py 1380 1410 tic0010756751, python fourier_analysis_HW4.py 2440 2453 tic0010756751 

-- 10/12/25
'''

#Creates parser
parser = argparse.ArgumentParser()

#Creates arguments for a, b, num_points, N, and increment_step
parser.add_argument('lower_bound', type= float, help= 'The lower bound of time epoch of interest.')
parser.add_argument('upper_bound', type= float, help= 'The upper bound of time epoch of interest.')
parser.add_argument('dat', type= str, help= 'The system_name, all lower case w/o .fits ')

#Parses the arguments
args = parser.parse_args()

#Stores filename as a string
filename = f'{args.dat}.fits'

#Opens the fits file
hdul = fits.open(filename)

#Creates an astropy table of data taken in from the .fits file 
dat_table = Table(hdul[1].data)

#Defines a function to calculate the Fourier transform (hard-coded way) and inverse Fourier transform w/o interpolation, and plots the power spectrum and lightcurve
def fourier_transform(lower_bound, upper_bound, dat, sys_name):

    #Defines three subplots to visualize the original lightcurve, the power spectrum, and the inverse transform 
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(18, 5)) 

    #Plots original lightcurve (flux vs. time)
    ax1.plot(dat_table['times'], dat_table['fluxes'])
    ax1.set_xlim(lower_bound, upper_bound) 
    ax1.set_xlabel('Time (JD)')
    ax1.set_ylabel('Flux')
    ax1.set_title(f'Light Curve for {sys_name.upper()}')
    

    #Checks the times column to determine where the times are larger than the user-inputted lower bound and smaller than the upper-bound
    condition_dat_table = (dat['times'] > lower_bound) & (dat['times'] < upper_bound)

    #Filters original astropy table to only include points from epoch of interest
    filtered_dat_table = dat[condition_dat_table]

    #Defines function to calculate the full Discrete Fourier transform of the filtered flux points
    def dft(fluxes):

        N = len(fluxes)

        c = zeros(N, complex)

        for k in range(N):
            for n in range(N):
                c[k] += fluxes[n] * exp(-2j*pi*k*n/N)
        
        return c, N

    #Calls the dft function to calculate the discrete fourier transform of the filtered fluxes 
    discrete_fourier_transform, N = dft(filtered_dat_table['fluxes'])
    
    discrete_fourier_transform[0] = 0  #Sets k = 0 term to 0

    #Defines k values array for first 100 values
    k_vals = np.arange(0, 100, 1)
    
    #Calculates the coefficients (c_k) squared from the DFT
    dft_squared = discrete_fourier_transform * np.conjugate(discrete_fourier_transform)

    #Filters and stores the first 100 squared coefficients
    filt_dft_squared = dft_squared[:100]

    #Plots the Power Spectrum
    ax2.bar(k_vals, filt_dft_squared)
    ax2.set_xlim(0, 100)
    ax2.set_xlabel('$k$')
    ax2.set_ylabel('$c_k^2$')
    ax2.set_title('Power Spectrum')
        
    #Defines function to calculate the full inverse Discrete Fourier transform of the filtered flux points
    def inverse_dft(c, N):

        y = zeros(N, complex)

        for k in range(N):
            for n in range(N):
                y[n] += c[k] * exp(2j*pi*k*n/N)

        return y/N
    
    #Calls the dft inverse function to calculate the discrete fourier transform of the filtered fluxes 
    flux_vals = inverse_dft(discrete_fourier_transform, len(filtered_dat_table['times']))

    #Plots the Inverse Discrete Fourier Transform 
    ax3.plot(filtered_dat_table['times'], flux_vals)
    ax3.set_xlim(lower_bound, upper_bound)
    ax3.set_xlabel("Time (JD)")
    ax3.set_ylabel("Flux")
    ax3.set_title("Inverse Discrete Fourier Transform ($c_k$ to flux)")

    #Save the entire figure as a .png file
    fig.savefig(f'dft_range_{lower_bound}_{upper_bound}.png', dpi=300)

    plt.show()

#Defines a function to calculate the Fast Fourier transform (using the numpy fft method) and inverse Fast Fourier transform w/ linear interpolation, and plots the power spectrum and lightcurve
def interp_dft(lower_bound, upper_bound, dat, sys_name):

    #Defines three subplots to visualize the original lightcurve, the power spectrum, and the inverse transform 
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(18, 5)) 
    
    #Checks the times column to determine where the times are larger than the user-inputted lower bound and smaller than the upper-bound
    condition_dat_table = (dat['times'] > lower_bound) & (dat['times'] < upper_bound)

    #Filters original astropy table to only include points from epoch of interest
    filtered_dat_table = dat[condition_dat_table]

    #Save times and fluxes
    times = filtered_dat_table['times']
    fluxes = filtered_dat_table['fluxes']

    #Create a linearized time grid
    uniform_times = np.linspace(times.min(), times.max(), len(times))

    #Linearly interpolate the fluxes relative to the linearized time grid
    interpolated_fluxes = np.interp(uniform_times, times, fluxes) 

    #Plots original lightcurve (flux vs. time) w/ interpolated fluxes
    ax1.plot(uniform_times, fluxes, color = "green")
    ax1.set_xlim(lower_bound, upper_bound) 
    ax1.set_xlabel('Time (JD)')
    ax1.set_ylabel('Flux')
    ax1.set_title(f'Light Curve for {sys_name.upper()}')

    #Compute the fast fourier transform of the interpolated fluxes
    fast_fourier = np.fft.fft(interpolated_fluxes)
    
    #Calculates the coefficients (c_k) squared from the FFT
    fft_squared = fast_fourier * np.conjugate(fast_fourier)
    fft_squared[0] = 0  #Sets k = 0 term to 0

    #Defines k values array for first 100 values
    k_vals = np.arange(0, 100, 1)

    #Filters and stores the first 100 squared coefficients
    filt_fft_squared = fft_squared[:100]

    #Plots the Power Spectrum
    ax2.bar(k_vals, filt_fft_squared, color = "green")
    ax2.set_xlim(0, 100)
    ax2.set_xlabel('$k$')
    ax2.set_ylabel('$c_k^2$')
    ax2.set_title('Power Spectrum')

    #Attempted to find the max index and filter out points to do inverse FFT but was not sure how to shift times
    '''
    index = filt_fft_squared.argmax()

    
    keep = filt_fft_squared > 0.1*filt_fft_squared[index]

    c_new = np.copy(fast_fourier)
    c_new[~keep] = 0.0  #set the indices not to be kept to 0.0
    '''
    
    #Compute the inverse FFT 
    inverse_fast_fourier = np.fft.ifft(fast_fourier)

    #Plots the Inverse Discrete Fourier Transform 
    ax3.plot(uniform_times, inverse_fast_fourier, color = "green")
    ax3.set_xlim(lower_bound, upper_bound)
    ax3.set_xlabel("Time (JD)")
    ax3.set_ylabel("Flux")
    ax3.set_title("Inverse Fast Fourier Transform ($c_k$ to flux)")

    #Save the entire figure as a .png file
    fig.savefig(f'fft_range_{lower_bound}_{upper_bound}.png', dpi=300)
    plt.show()

fourier_transform(args.lower_bound, args.upper_bound, dat_table, args.dat) 
interp_dft(args.lower_bound, args.upper_bound, dat_table, args.dat)












