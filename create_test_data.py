import numpy as np
from astropy.io import fits
import os

def create_test_spectrum(filename, time_offset=0):
    """Create a test spectrum with realistic features"""
    # Create wavelength array (1150-1700 Å)
    wavelength = np.linspace(1150, 1700, 1024)
    
    # Create baseline flux
    flux = np.exp(-(wavelength - 1425)**2 / 100000)
    
    # Add emission lines
    emission_lines = {
        'Ly-α': 1216,
        'NV': 1240,
        'OI': 1304,
        'CII': 1335,
        'OI_2': 1356
    }
    
    for wav in emission_lines.values():
        # Add emission line with time-variable strength
        strength = 1.0 + 0.2 * np.sin(time_offset)
        flux += strength * np.exp(-(wavelength - wav)**2 / 2)
    
    # Add noise
    noise = np.random.normal(0, 0.05, len(wavelength))
    flux += noise
    
    # Create error array
    error = np.sqrt(np.abs(flux)) + 0.1
    
    # Create primary HDU
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header['DATE-OBS'] = f'2024-03-06T{10+time_offset:02d}:00:00'
    primary_hdu.header['TELESCOP'] = 'HST'
    primary_hdu.header['INSTRUME'] = 'STIS'
    
    # Create data table
    col1 = fits.Column(name='WAVELENGTH', format='E', array=[wavelength])
    col2 = fits.Column(name='FLUX', format='E', array=[flux])
    col3 = fits.Column(name='ERROR', format='E', array=[error])
    
    table_hdu = fits.BinTableHDU.from_columns([col1, col2, col3])
    
    # Create HDU list and write file
    hdul = fits.HDUList([primary_hdu, table_hdu])
    hdul.writeto(filename, overwrite=True)

def create_test_dataset():
    """Create a set of test files"""
    # Create output directory if it doesn't exist
    os.makedirs('test_data/europa_uv_tutorial/raw', exist_ok=True)
    
    # Create three observations at different times
    for i, offset in enumerate([0, 2, 4]):
        filename = f'test_data/europa_uv_tutorial/raw/test_obs_{i+1}_x1d.fits'
        create_test_spectrum(filename, offset)
        print(f"Created {filename}")

if __name__ == '__main__':
    create_test_dataset() 