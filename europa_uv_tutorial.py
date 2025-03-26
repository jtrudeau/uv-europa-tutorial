#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% [markdown]
# # Europa UV Spectral Analysis Tutorial
# 
# ## Introduction
# 
# Welcome to this hands-on tutorial on analyzing ultraviolet (UV) spectra of Jupiter's moon Europa! In this tutorial, you will:
# 1. Learn how to work with astronomical data
# 2. Analyze real UV spectra from the Hubble Space Telescope (HST)
# 3. Discover what these spectra tell us about Europa
# 
# ### Why Study Europa's UV Emissions?
# Europa is one of the most intriguing moons in our solar system:
# - It likely has a subsurface ocean of liquid water
# - It might have conditions suitable for life
# - It shows signs of current activity (water vapor plumes)
# 
# By studying UV emissions, we can learn about:
# - The composition of Europa's tenuous atmosphere
# - Potential plume activity
# - Interactions with Jupiter's magnetic field

# %% [markdown]
# ## Setup
# 
# First, let's install and import the packages we need:

# %%
# Install required packages
# !pip install astropy astroquery numpy matplotlib pandas

# Import packages
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astroquery.mast import Observations
import pandas as pd
import os
from astropy.time import Time

# Set up plotting style
plt.style.use('ggplot')  # Using a built-in style instead of seaborn
# %matplotlib inline

# %% [markdown]
# ## Part 1: Understanding FITS Files
# 
# FITS (Flexible Image Transport System) is the standard file format in astronomy. Let's create a simple example to understand its structure:

# %%
# Create a simple example spectrum
wavelength = np.linspace(100, 200, 1000)  # wavelength in nm
flux = np.exp(-(wavelength - 150)**2 / 100)  # example spectral shape
noise = np.random.normal(0, 0.05, len(wavelength))
flux += noise

# Plot the example spectrum
plt.figure(figsize=(10, 6))
plt.plot(wavelength, flux, 'b-', label='Example Spectrum')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Flux (arbitrary units)')
plt.title('Example UV Spectrum')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

# %% [markdown]
# ### Exercise 1: Create and Examine a FITS File
# 
# Let's save our example spectrum as a FITS file and examine its structure:

# %%
# Create a FITS file
def create_example_fits():
    # Create primary HDU
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header['OBSERVER'] = 'Student'  # replace with your name
    primary_hdu.header['DATE-OBS'] = '2025-03-06'
    
    # Create data table - fixing the array shape issue
    col1 = fits.Column(name='WAVELENGTH', format='E', array=wavelength)
    col2 = fits.Column(name='FLUX', format='E', array=flux)
    table_hdu = fits.BinTableHDU.from_columns([col1, col2])
    
    # Create HDU list and write file
    hdul = fits.HDUList([primary_hdu, table_hdu])
    hdul.writeto('example_spectrum.fits', overwrite=True)

# Create the file
create_example_fits()

# Examine the file
with fits.open('example_spectrum.fits') as hdul:
    print("File structure:")
    hdul.info()
    
    print("\nHeader information:")
    for key, value in hdul[0].header.items():
        print(f"{key}: {value}")

# %% [markdown]
# ## Part 2: Accessing Real HST Data
# 
# Now let's search for real HST observations of Europa:

# %%
# Search for Europa observations
obs_table = Observations.query_criteria(target_name="Europa", obs_collection="HST")
print(f"Found {len(obs_table)} HST observations of Europa")

# Convert to pandas DataFrame for easier analysis
df = obs_table.to_pandas()

# Show summary of instruments used
print("\nInstruments used to observe Europa:")
print(df['instrument_name'].value_counts())

# %% [markdown]
# ### Exercise 2: Find UV Observations
# 
# Let's filter for STIS UV observations:

# %%
# Filter for STIS observations
stis_obs = df[df['instrument_name'].str.contains('STIS', na=False)]
print(f"Found {len(stis_obs)} STIS observations")

# Look for UV observations (FUV or NUV)
uv_obs = stis_obs[stis_obs['instrument_name'].str.contains('FUV|NUV', na=False)]
print(f"Found {len(uv_obs)} UV observations")

# Display details of UV observations
print("\nUV Observation Details:")
print(uv_obs[['obs_id', 'instrument_name', 'filters', 't_min']].head())

# %% [markdown]
# ## Part 3: Analyzing Real HST Data
# 
# Now let's analyze the real HST data we've downloaded:

# %%
def mjd_to_readable_date(mjd):
    """Convert Modified Julian Date to a human-readable date string"""
    try:
        t = Time(mjd, format='mjd')
        return t.iso
    except:
        return f"MJD {mjd}"

def analyze_spectrum(filename):
    """Analyze a UV spectrum"""
    with fits.open(filename) as hdul:
        # Get observation details from primary and extension headers
        primary_header = hdul[0].header
        ext_header = hdul[1].header
        
        # For HST STIS data:
        # - EXPSTART in extension header: MJD (Modified Julian Date) of exposure start
        # - TEXPTIME in extension header: Total exposure time in seconds
        
        # Get observation date (using EXPSTART from extension header if available)
        if 'EXPSTART' in ext_header:
            # Get MJD value
            mjd = ext_header['EXPSTART']
            # Convert to readable date
            obs_date = mjd_to_readable_date(mjd)
        else:
            # Fallback to DATE from primary header (file creation date)
            obs_date = primary_header.get('DATE', 'Unknown')
        
        # Get exposure time (using TEXPTIME from extension header if available)
        if 'TEXPTIME' in ext_header:
            exp_time = ext_header['TEXPTIME']
        elif 'EXPTIME' in ext_header:
            exp_time = ext_header['EXPTIME']
        else:
            # Fallback to EXPTIME from primary header
            exp_time = primary_header.get('EXPTIME', 0)
        
        # Print observation details for this file
        print(f"\nFile: {os.path.basename(filename)}")
        print(f"Observation date: {obs_date}")
        print(f"Exposure time: {exp_time} seconds")
        
        # Get the observation ID from the filename
        obs_id = os.path.basename(filename).split('_')[0]
        
        # Get spectrum data
        spec_data = hdul[1].data
        wavelength = np.array(spec_data['WAVELENGTH'][0])  # Angstroms
        flux = np.array(spec_data['FLUX'][0])
        error = np.array(spec_data['ERROR'][0])
        
        # Convert wavelength to nm
        wavelength_nm = wavelength / 10.0
        
        return wavelength_nm, flux, error, obs_date, exp_time, obs_id

# Define paths to our real HST data
data_dir = "data/hst/8224/calibrated/x1d"
x1d_files = [
    os.path.join(data_dir, "o5d601010_x1d.fits"),  # First observation
    os.path.join(data_dir, "o5d601080_x1d.fits"),  # Middle observation
    os.path.join(data_dir, "o5d6010a0_x1d.fits")   # Last observation
]

# Check if the real data exists
if os.path.exists(x1d_files[0]):
    # Analyze first observation
    wave, flux, error, date, exp_time, obs_id = analyze_spectrum(x1d_files[0])

    # Create plot
    plt.figure(figsize=(12, 6))
    plt.plot(wave, flux, 'b-', label='Flux')
    plt.fill_between(wave, flux-error, flux+error, color='blue', alpha=0.3, label='Error')

    # Add emission lines
    lines = {
        'H Lyman-α': 121.6,
        'O I': 130.4,
        'C II': 133.5,
        'O I (2)': 135.6,
        'N V': 124.0
    }

    for name, wav in lines.items():
        plt.axvline(x=wav, color='r', linestyle='--', label=name)

    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Flux (erg/s/cm²/Å)')
    plt.title(f'Europa UV Spectrum\nObservation ID: {obs_id}, Exposure: {exp_time}s')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
else:
    print(f"Real HST data not found at {x1d_files[0]}")
    print("Please make sure you've downloaded the data to the correct location.")

# %% [markdown]
# ### Exercise 3: Analyze Emission Features
# 
# Let's measure the strength of emission features in the real data:

# %%
def measure_emission_line(wavelength, flux, error, line_center, window=1.0):
    """Measure properties of an emission line"""
    # Select region around line
    mask = np.abs(wavelength - line_center) < window
    line_wave = wavelength[mask]
    line_flux = flux[mask]
    line_error = error[mask]
    
    # Calculate properties
    peak_flux = np.max(line_flux)
    total_flux = np.sum(line_flux)
    snr = peak_flux / np.mean(line_error)
    
    return {
        'peak_flux': peak_flux,
        'total_flux': total_flux,
        'snr': snr
    }

# Check if real data exists
if os.path.exists(x1d_files[0]):
    # Measure emission lines
    for name, wav in lines.items():
        results = measure_emission_line(wave, flux, error, wav)
        print(f"\n{name} ({wav} nm):")
        print(f"Peak flux: {results['peak_flux']:.2e}")
        print(f"Total flux: {results['total_flux']:.2e}")
        print(f"Signal-to-noise: {results['snr']:.1f}")

# %% [markdown]
# ## Part 4: Time Series Analysis
# 
# Let's analyze how the spectrum changes across our three observations:

# %%
def analyze_time_series(file_list):
    """Analyze changes in emission lines across observations"""
    # Check if all files exist
    if not all(os.path.exists(f) for f in file_list):
        print("Some data files not found. Please check paths.")
        return
    
    obs_ids = []
    dates = []
    mjd_values = []
    lya_flux = []
    oi_flux = []
    
    # Analyze each observation
    for i, filename in enumerate(file_list):
        wave, flux, error, date, exp_time, obs_id = analyze_spectrum(filename)
        
        # Measure emission lines
        lya = measure_emission_line(wave, flux, error, 121.6)
        oi = measure_emission_line(wave, flux, error, 130.4)
        
        # Extract MJD value if available
        mjd = None
        if "MJD" in date:
            try:
                mjd = float(date.split(" ")[1])
            except:
                pass
        
        obs_ids.append(obs_id)
        dates.append(date)
        if mjd is not None:
            mjd_values.append(mjd)
        lya_flux.append(lya['peak_flux'])
        oi_flux.append(oi['peak_flux'])
    
    # Create time series plot
    plt.figure(figsize=(10, 6))
    
    # If we have MJD values, use them for x-axis
    if len(mjd_values) == len(lya_flux):
        plt.plot(mjd_values, lya_flux, 'ro-', label='Lyman-α')
        plt.plot(mjd_values, oi_flux, 'bo-', label='O I')
        plt.xlabel('Modified Julian Date (MJD)')
        
        # Add a second x-axis with human-readable dates
        ax1 = plt.gca()
        ax2 = ax1.twiny()
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(mjd_values)
        ax2.set_xticklabels([d.split('T')[0] for d in dates], rotation=45)
    else:
        # Otherwise use observation IDs
        plt.plot(range(len(obs_ids)), lya_flux, 'ro-', label='Lyman-α')
        plt.plot(range(len(obs_ids)), oi_flux, 'bo-', label='O I')
        plt.xlabel('Observation')
        plt.xticks(range(len(obs_ids)), obs_ids, rotation=45)
    
    plt.ylabel('Peak Flux (erg/s/cm²/Å)')
    plt.title('Emission Line Strength vs. Time')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Run time series analysis if data exists
if all(os.path.exists(f) for f in x1d_files):
    analyze_time_series(x1d_files)
else:
    print("Not all observation files found. Time series analysis skipped.")

# %% [markdown]
# ## Advanced Analysis: Statistical Methods
# 
# Let's apply some statistical methods to our real HST data:

# %%
from scipy.stats import linregress
from scipy.optimize import curve_fit

def gaussian(x, amplitude, center, width):
    """Gaussian function for fitting emission lines"""
    return amplitude * np.exp(-(x - center)**2 / (2 * width**2))

def advanced_line_analysis(wavelength, flux, error, line_center, window=2.0):
    """Perform advanced analysis of an emission line"""
    # Select data around the line
    mask = np.abs(wavelength - line_center) < window
    wave_subset = wavelength[mask]
    flux_subset = flux[mask]
    error_subset = error[mask]
    
    try:
        # Fit Gaussian to the line
        p0 = [np.max(flux_subset), line_center, 0.5]  # Initial guess
        popt, pcov = curve_fit(gaussian, wave_subset, flux_subset, p0=p0)
        
        # Calculate line properties
        amplitude, center, width = popt
        
        # Calculate integrated flux
        integrated_flux = amplitude * width * np.sqrt(2 * np.pi)
        
        # Plot the results
        plt.figure(figsize=(10, 6))
        plt.errorbar(wave_subset, flux_subset, yerr=error_subset, fmt='o', label='Data')
        
        # Plot the fit
        x_fit = np.linspace(min(wave_subset), max(wave_subset), 100)
        y_fit = gaussian(x_fit, *popt)
        plt.plot(x_fit, y_fit, 'r-', label='Gaussian Fit')
        
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Flux')
        plt.title(f'Emission Line Analysis: {line_center} nm')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()
        
        return {
            'amplitude': amplitude,
            'center': center,
            'width': width,
            'integrated_flux': integrated_flux
        }
    except Exception as e:
        print(f"Error fitting line: {e}")
        return None

# Check if real data exists
if os.path.exists(x1d_files[0]):
    # Perform advanced analysis on Lyman-alpha line
    wave, flux, error, date, exp_time, obs_id = analyze_spectrum(x1d_files[0])
    results = advanced_line_analysis(wave, flux, error, 121.6)
    
    if results:
        print("\nAdvanced Lyman-α Analysis:")
        print(f"Center: {results['center']:.2f} nm")
        print(f"Width: {results['width']:.2f} nm")
        print(f"Integrated flux: {results['integrated_flux']:.2e}")

# %% [markdown]
# ## Comparing Multiple Observations
# 
# Let's compare the spectra from all three observations:

# %%
def compare_observations(file_list):
    """Compare spectra from multiple observations"""
    plt.figure(figsize=(12, 8))
    
    for i, filename in enumerate(file_list):
        if os.path.exists(filename):
            wave, flux, error, date, exp_time, obs_id = analyze_spectrum(filename)
            
            # Plot with different colors
            colors = ['b', 'g', 'r']
            plt.plot(wave, flux, color=colors[i], label=f'Obs: {obs_id} (Exp: {exp_time:.0f}s)')
    
    # Add emission lines
    for name, wav in lines.items():
        plt.axvline(x=wav, color='k', linestyle='--', alpha=0.5)
        plt.text(wav, plt.ylim()[1]*0.9, name, rotation=90, ha='right')
    
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Flux (erg/s/cm²/Å)')
    plt.title('Comparison of Europa UV Spectra')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

# Compare observations if data exists
if all(os.path.exists(f) for f in x1d_files):
    compare_observations(x1d_files)
else:
    print("Not all observation files found. Comparison skipped.")

# %% [markdown]
# ## Part 5: Scientific Interpretation
# 
# Now let's interpret what we've found:
# 
# 1. **Emission Lines**:
#    - Lyman-α (121.6 nm): Indicates hydrogen atoms
#    - O I (130.4 nm): Shows presence of oxygen
#    - C II (133.5 nm): Carbon in the atmosphere
# 
# 2. **Time Variations**:
#    - Changes might indicate:
#      * Plume activity
#      * Magnetospheric effects
#      * Changes in solar illumination
# 
# 3. **Implications**:
#    - Presence of water (H and O)
#    - Active processes on Europa
#    - Interaction with Jupiter's environment
# 
# ### Exercise 4: Write a Scientific Summary
# 
# Based on your analysis:
# 1. What emission features did you detect?
# 2. How do they vary with time?
# 3. What might cause these variations?
# 4. What does this tell us about Europa?

# %% [markdown]
# ## Part 6: Connection to Future Missions
# 
# The Europa Clipper mission (launching in 2024) will study Europa in unprecedented detail. Let's explore how our analysis connects to this mission:
# 
# 1. **Europa-UVS Instrument**:
#    - Will observe similar UV emissions
#    - Much closer to Europa (better resolution)
#    - Can track changes over multiple flybys
# 
# 2. **Science Goals**:
#    - Search for plume activity
#    - Study atmospheric composition
#    - Investigate surface composition
# 
# 3. **Predictions**:
#    - Based on our HST analysis, we might expect:
#      * Stronger emission features
#      * Better spatial resolution
#      * Temporal variations related to orbital position

# %% [markdown]
# ## Next Steps
# 
# You can extend this analysis by:
# 1. Analyzing more observations
# 2. Studying correlations with Europa's orbital position
# 3. Comparing with other satellites
# 4. Preparing for Europa Clipper observations
# 
# ### Resources
# - [HST Data Handbook](https://hst-docs.stsci.edu/)
# - [STIS Instrument Handbook](https://www.stsci.edu/hst/stis)
# - [Europa Clipper Mission](https://europa.nasa.gov/)
# - [Recent Europa Publications](https://ui.adsabs.harvard.edu/)

# %% [markdown]
# ## Conclusion
# 
# In this tutorial, you've learned:
# 1. How to work with FITS files
# 2. How to access HST data
# 3. How to analyze UV spectra
# 4. How to interpret the results scientifically
# 
# These skills are valuable for studying not just Europa, but many other astronomical objects as well. 

# %% [markdown]
# ## Understanding HST Data Headers
# 
# HST data uses Modified Julian Date (MJD) for timestamps:
# - MJD = JD - 2400000.5
# - MJD 51456 corresponds to October 5, 1999
# - The decimal part represents the fraction of the day
# 
# For example:
# - MJD 51456.36073172 is approximately 8:39 UT on October 5, 1999
# - MJD 51456.63515321 is approximately 15:14 UT on October 5, 1999
# 
# These observations were taken over about 6.5 hours, showing Europa at different orbital positions. 