from astroquery.mast import Observations
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
from datetime import datetime

# First, let's analyze the downloaded file if it exists
download_path = "downloads"
fits_file = os.path.join(download_path, "o5d601040_asn.fits")

# List of files we should look for based on the association file
member_files = []

if os.path.exists(fits_file):
    print(f"\nAnalyzing downloaded file: {fits_file}")
    try:
        with fits.open(fits_file) as hdul:
            # Print file structure
            print("\nFITS file structure:")
            hdul.info()
            
            # Print header information from the primary HDU
            print("\nHeader information:")
            for key, value in hdul[0].header.items():
                if key not in ['COMMENT', 'HISTORY']:
                    print(f"{key}: {value}")
            
            # Check if there's data in the primary HDU
            if hdul[0].data is not None:
                print("\nPrimary HDU contains data with shape:", hdul[0].data.shape)
            
            # If there are extensions, examine them
            for i in range(1, len(hdul)):
                print(f"\nExtension {i} information:")
                if hasattr(hdul[i], 'columns') and hdul[i].columns is not None:
                    print("Column names:", hdul[i].columns.names)
                if hdul[i].data is not None:
                    if hasattr(hdul[i].data, 'shape'):
                        print(f"Data shape: {hdul[i].data.shape}")
                    else:
                        print("Data is present but has no shape attribute")
                        
            # ASN files are association files that point to other files
            # Let's check if this is an association file
            if 'ASN_TAB' in hdul[0].header:
                print("\nThis is an association file that points to these member files:")
                if len(hdul) > 1 and hasattr(hdul[1], 'data'):
                    for row in hdul[1].data:
                        if 'MEMNAME' in hdul[1].columns.names:
                            member = row['MEMNAME'].strip()
                            print(f"- {member}")
                            member_files.append(member)
                        elif 'MEMNAME' in hdul[1].data.names:
                            member = row['MEMNAME'].strip()
                            print(f"- {member}")
                            member_files.append(member)
    except Exception as e:
        print(f"Error analyzing FITS file: {e}")

# Now let's look for any spectral data files in the downloads directory
print("\nSearching for spectral data files in downloads directory...")
all_fits_files = glob.glob(os.path.join(download_path, "*.fits"))
print(f"Found {len(all_fits_files)} FITS files in downloads directory:")
for file in all_fits_files:
    print(f"- {os.path.basename(file)}")

# Look for specific file types that might contain spectral data
spectral_files = []
for file in all_fits_files:
    basename = os.path.basename(file)
    # Look for x1d, x2d, tag files which might contain spectral data
    if any(ext in basename for ext in ['_x1d', '_x2d', '_tag', '_flt', '_raw']):
        spectral_files.append(file)

if spectral_files:
    print(f"\nFound {len(spectral_files)} potential spectral data files:")
    for file in spectral_files:
        print(f"- {os.path.basename(file)}")
    
    # Try to analyze the first spectral file
    print(f"\nAnalyzing first spectral file: {os.path.basename(spectral_files[0])}")
    try:
        with fits.open(spectral_files[0]) as hdul:
            # Print file structure
            print("\nFITS file structure:")
            hdul.info()
            
            # Determine what kind of file this is
            file_type = None
            if '_x1d' in spectral_files[0]:
                file_type = 'x1d'
            elif '_x2d' in spectral_files[0]:
                file_type = 'x2d'
            elif '_tag' in spectral_files[0]:
                file_type = 'tag'
            elif '_flt' in spectral_files[0]:
                file_type = 'flt'
            elif '_raw' in spectral_files[0]:
                file_type = 'raw'
            
            print(f"\nFile type: {file_type}")
            
            # Extract and plot data based on file type
            if file_type == 'x1d':
                # For 1D extracted spectra
                if len(hdul) > 1 and hdul[1].data is not None:
                    spec_data = hdul[1].data
                    print("\nSpectrum data columns:", spec_data.names)
                    
                    # Extract wavelength and flux
                    wavelength = spec_data['WAVELENGTH']  # Typically in Angstroms
                    flux = spec_data['FLUX']
                    error = spec_data['ERROR'] if 'ERROR' in spec_data.names else np.sqrt(np.abs(flux))
                    
                    # Convert wavelength to nanometers for better readability
                    wavelength_nm = wavelength / 10.0
                    
                    # Create the plot
                    plt.figure(figsize=(12, 6))
                    plt.plot(wavelength_nm, flux, 'b-', label='Flux')
                    plt.fill_between(wavelength_nm, flux-error, flux+error, color='blue', alpha=0.3, label='Error')
                    
                    # Add vertical lines for key emission features
                    lyman_alpha = 121.6  # nm
                    oxygen_line = 130.4  # nm
                    plt.axvline(x=lyman_alpha, color='r', linestyle='--', label='H Lyman-α (121.6 nm)')
                    plt.axvline(x=oxygen_line, color='g', linestyle='--', label='O I (130.4 nm)')
                    
                    # Add labels and title
                    plt.xlabel('Wavelength (nm)')
                    plt.ylabel('Flux (erg/s/cm²/Å)')
                    plt.title('HST STIS UV Spectrum of Europa')
                    plt.legend()
                    plt.grid(True, alpha=0.3)
                    
                    # Save the plot
                    plt.savefig(os.path.join(download_path, 'europa_uv_spectrum.png'))
                    print(f"\nSpectrum plot saved to {os.path.join(download_path, 'europa_uv_spectrum.png')}")
                    
                    # Analyze the spectrum
                    print("\nSpectrum Analysis:")
                    print(f"Wavelength range: {wavelength_nm.min():.1f} - {wavelength_nm.max():.1f} nm")
                    print(f"Mean flux: {np.mean(flux):.2e} erg/s/cm²/Å")
                    print(f"Max flux: {np.max(flux):.2e} erg/s/cm²/Å")
                    
                    # Check for emission features near Lyman-alpha and oxygen line
                    lyman_alpha_idx = np.abs(wavelength_nm - lyman_alpha).argmin()
                    oxygen_idx = np.abs(wavelength_nm - oxygen_line).argmin()
                    
                    lyman_alpha_flux = flux[lyman_alpha_idx]
                    oxygen_flux = flux[oxygen_idx]
                    
                    print(f"Flux near H Lyman-α (121.6 nm): {lyman_alpha_flux:.2e} erg/s/cm²/Å")
                    print(f"Flux near O I (130.4 nm): {oxygen_flux:.2e} erg/s/cm²/Å")
                    
                    # Determine if we're seeing emission or absorption
                    local_bg_lyman = np.median(flux[max(0, lyman_alpha_idx-5):min(len(flux), lyman_alpha_idx+6)])
                    local_bg_oxygen = np.median(flux[max(0, oxygen_idx-5):min(len(flux), oxygen_idx+6)])
                    
                    if lyman_alpha_flux > local_bg_lyman * 1.2:
                        print(f"Lyman-α appears to show emission (signal is {lyman_alpha_flux/local_bg_lyman:.1f}x local background)")
                    elif lyman_alpha_flux < local_bg_lyman * 0.8:
                        print(f"Lyman-α appears to show absorption (signal is {lyman_alpha_flux/local_bg_lyman:.1f}x local background)")
                    
                    if oxygen_flux > local_bg_oxygen * 1.2:
                        print(f"O I appears to show emission (signal is {oxygen_flux/local_bg_oxygen:.1f}x local background)")
                    elif oxygen_flux < local_bg_oxygen * 0.8:
                        print(f"O I appears to show absorption (signal is {oxygen_flux/local_bg_oxygen:.1f}x local background)")
                    
                    # Check for overall UV reflectance trend
                    if wavelength_nm.max() - wavelength_nm.min() > 50:  # If we have a wide wavelength range
                        # Smooth the spectrum for trend analysis
                        window_size = max(5, len(flux) // 20)
                        smoothed_flux = np.convolve(flux, np.ones(window_size)/window_size, mode='valid')
                        smoothed_wave = wavelength_nm[window_size-1:]
                        
                        # Check if flux decreases toward shorter wavelengths (UV downturn)
                        if len(smoothed_flux) > 10:
                            first_quarter = smoothed_flux[:len(smoothed_flux)//4]
                            last_quarter = smoothed_flux[-len(smoothed_flux)//4:]
                            
                            if np.median(first_quarter) < np.median(last_quarter) * 0.8:
                                print("\nThe spectrum shows a UV downturn (lower reflectance at shorter wavelengths),")
                                print("which is consistent with the presence of sulfur compounds or other UV absorbers on Europa's surface.")
            
            elif file_type == 'x2d':
                # For 2D rectified spectral images
                if len(hdul) > 1 and hdul[1].data is not None:
                    image_data = hdul[1].data
                    print(f"\nImage data shape: {image_data.shape}")
                    
                    # Extract wavelength information from header
                    if 'CRVAL1' in hdul[1].header and 'CD1_1' in hdul[1].header and 'CRPIX1' in hdul[1].header:
                        crval1 = hdul[1].header['CRVAL1']
                        cd1_1 = hdul[1].header['CD1_1']
                        crpix1 = hdul[1].header['CRPIX1']
                        
                        # Calculate wavelength array
                        wavelength = crval1 + cd1_1 * (np.arange(image_data.shape[1]) - crpix1)
                        wavelength_nm = wavelength / 10.0  # Convert to nm
                        
                        # Sum along spatial axis to get 1D spectrum
                        flux = np.sum(image_data, axis=0)
                        
                        # Create the plot
                        plt.figure(figsize=(12, 6))
                        plt.plot(wavelength_nm, flux, 'b-', label='Flux')
                        
                        # Add vertical lines for key emission features
                        lyman_alpha = 121.6  # nm
                        oxygen_line = 130.4  # nm
                        plt.axvline(x=lyman_alpha, color='r', linestyle='--', label='H Lyman-α (121.6 nm)')
                        plt.axvline(x=oxygen_line, color='g', linestyle='--', label='O I (130.4 nm)')
                        
                        # Add labels and title
                        plt.xlabel('Wavelength (nm)')
                        plt.ylabel('Summed Flux')
                        plt.title('HST STIS UV Spectrum of Europa (from 2D image)')
                        plt.legend()
                        plt.grid(True, alpha=0.3)
                        
                        # Save the plot
                        plt.savefig(os.path.join(download_path, 'europa_uv_spectrum_2d.png'))
                        print(f"\nSpectrum plot saved to {os.path.join(download_path, 'europa_uv_spectrum_2d.png')}")
                        
                        # Also create a 2D image plot
                        plt.figure(figsize=(12, 6))
                        plt.imshow(image_data, aspect='auto', origin='lower', 
                                  extent=[wavelength_nm.min(), wavelength_nm.max(), 0, image_data.shape[0]])
                        plt.colorbar(label='Flux')
                        plt.xlabel('Wavelength (nm)')
                        plt.ylabel('Spatial Position')
                        plt.title('HST STIS UV 2D Spectral Image of Europa')
                        
                        # Add vertical lines for key emission features
                        plt.axvline(x=lyman_alpha, color='r', linestyle='--', label='H Lyman-α (121.6 nm)')
                        plt.axvline(x=oxygen_line, color='g', linestyle='--', label='O I (130.4 nm)')
                        plt.legend()
                        
                        # Save the plot
                        plt.savefig(os.path.join(download_path, 'europa_uv_2d_image.png'))
                        print(f"\n2D image plot saved to {os.path.join(download_path, 'europa_uv_2d_image.png')}")
            
            elif file_type == 'tag':
                # For time-tagged photon events
                print("\nTAG files contain time-tagged photon events.")
                print("These files require specialized processing to extract spectral information.")
                print("They are useful for studying time-variable phenomena like auroral emissions.")
                
                # Print information about the extensions
                for i in range(len(hdul)):
                    if hdul[i].data is not None:
                        if hasattr(hdul[i].data, 'names'):
                            print(f"\nExtension {i} column names: {hdul[i].data.names}")
                        if hasattr(hdul[i].data, 'shape'):
                            print(f"Extension {i} data shape: {hdul[i].data.shape}")
            
            else:
                # For other file types
                print("\nThis file type requires specialized processing.")
                print("It may contain raw or partially processed data.")
                
                # Print information about the extensions
                for i in range(len(hdul)):
                    if hdul[i].data is not None:
                        if hasattr(hdul[i].data, 'shape'):
                            print(f"\nExtension {i} data shape: {hdul[i].data.shape}")
                            
                            # If it's a 2D image, plot it
                            if len(hdul[i].data.shape) == 2:
                                plt.figure(figsize=(10, 8))
                                plt.imshow(hdul[i].data, origin='lower', cmap='viridis')
                                plt.colorbar(label='Counts')
                                plt.title(f'HST STIS Data - {os.path.basename(spectral_files[0])} - Ext {i}')
                                plt.xlabel('Pixel')
                                plt.ylabel('Pixel')
                                
                                # Save the plot
                                plot_file = os.path.join(download_path, f'europa_ext{i}_image.png')
                                plt.savefig(plot_file)
                                print(f"\nImage plot saved to {plot_file}")
    
    except Exception as e:
        print(f"Error analyzing spectral file: {e}")
else:
    print("\nNo spectral data files found in downloads directory.")
    print("You may need to download the actual data files from MAST.")

# Print information about Europa's UV features
print("\n\nExpected features in Europa's UV spectrum:")
print("1. Hydrogen Lyman-α emission at 121.6 nm - Indicates presence of hydrogen in Europa's tenuous atmosphere")
print("2. Oxygen emission at 130.4 nm - Evidence of oxygen in the atmosphere, likely from water ice sputtering")
print("3. Overall low UV reflectance - Europa's surface is generally dark in UV")
print("4. Possible UV downturn at shorter wavelengths - Could indicate sulfur compounds on the surface")
print("5. Potential absorption features - May reveal composition of surface materials")

print("\nInterpreting the data:")
print("- If we see strong emission lines at 121.6 nm (H) and 130.4 nm (O), this indicates auroral activity")
print("  and confirms the presence of a tenuous atmosphere around Europa.")
print("- The presence of these emissions could be evidence of water molecules being dissociated")
print("  into hydrogen and oxygen, possibly from subsurface water escaping through cracks or plumes.")
print("- If we see a reflectance spectrum (no emission lines), the overall shape can tell us about")
print("  surface composition, particularly the presence of sulfur compounds from Io's plasma torus.")
print("- The G140L grating used in these observations covers approximately 115-170 nm, which includes")
print("  the key emission lines and the spectral range where sulfur compounds show distinctive absorption.")

# Query HST observations of Europa
print("\n\nQuerying HST observations of Europa...")
obs_table = Observations.query_criteria(target_name="Europa", obs_collection="HST")
print(f"\nFound {len(obs_table)} total HST observations of Europa.")

# Convert to pandas DataFrame for easier exploration
df = obs_table.to_pandas()

# Display summary of all instruments used to observe Europa
print("\nAll instruments used to observe Europa:")
instruments = df['instrument_name'].value_counts()
for instrument, count in instruments.items():
    print(f"- {instrument}: {count} observations")

# Filter for STIS observations
stis_obs = df[df['instrument_name'].str.contains('STIS', na=False)]
print(f"\nFound {len(stis_obs)} STIS observations of Europa.")

# Show STIS observation details
print("\nSTIS Observation Summary:")
print(f"Date range: {stis_obs['t_min'].min()} to {stis_obs['t_max'].max()}")
print(f"Filters/Gratings used: {', '.join(stis_obs['filters'].dropna().unique())}")

# Find UV observations (focusing on FUV and NUV)
stis_uv = stis_obs[stis_obs['instrument_name'].str.contains('FUV|NUV', na=False)]
print(f"\nFound {len(stis_uv)} STIS UV observations that might show UV features.")

# Display details about the UV observations
if len(stis_uv) > 0:
    print("\nSTIS UV Observations:")
    print(stis_uv[['obs_id', 'instrument_name', 'filters', 't_min', 'proposal_id', 'proposal_pi']].head().to_string())
    
    # Get products for the first UV observation
    print("\nGetting products for the first STIS UV observation...")
    first_obs = stis_uv.iloc[0]
    uv_products = Observations.get_product_list(first_obs['obsid'])
    uv_products_df = uv_products.to_pandas()
    
    # Print all available file types to understand what's available
    print("\nAll available file types:")
    file_types = uv_products_df['productFilename'].str.extract(r'_([a-z0-9]{3})\.fits', expand=False).dropna().unique()
    print(", ".join(file_types))
    
    # Look for spectral files in priority order
    # 1. First check for x1d files (calibrated 1D extracted spectra)
    x1d_files = uv_products_df[uv_products_df['productFilename'].str.contains(r'_x1d\.fits', case=False, regex=True)]
    
    # 2. Check for x2d files (2D rectified spectral images)
    x2d_files = uv_products_df[uv_products_df['productFilename'].str.contains(r'_x2d\.fits', case=False, regex=True)]
    
    # 3. Check for tag files (time-tagged photon events)
    tag_files = uv_products_df[uv_products_df['productFilename'].str.contains(r'_tag\.fits', case=False, regex=True)]
    
    # 4. Check for sx1/sx2 files (summed extractions)
    sx_files = uv_products_df[uv_products_df['productFilename'].str.contains(r'_sx[12]\.fits', case=False, regex=True)]
    
    print(f"\nFound {len(x1d_files)} x1d files, {len(x2d_files)} x2d files, {len(tag_files)} tag files, and {len(sx_files)} sx files.")
    
    # Determine which file type to use based on availability and priority
    target_files = None
    file_type_description = ""
    
    if len(x1d_files) > 0:
        target_files = x1d_files
        file_type_description = "x1d (1D extracted spectra) - Best for analyzing spectral features like emission lines"
    elif len(x2d_files) > 0:
        target_files = x2d_files
        file_type_description = "x2d (2D rectified spectral images) - Good for seeing spatial distribution of spectral features"
    elif len(sx_files) > 0:
        target_files = sx_files
        file_type_description = "sx1/sx2 (summed extractions) - Alternative calibrated spectra"
    elif len(tag_files) > 0:
        target_files = tag_files
        file_type_description = "tag (time-tagged photon events) - Good for time-variable phenomena like auroral emissions"
    
    if target_files is not None:
        print(f"\nUsing {file_type_description}")
        print("\nAvailable files:")
        print(target_files[['productFilename', 'description', 'size']].to_string())
        
        # Print information about how to access the data
        print("\nTo download and analyze these files:")
        print("1. Visit the MAST Portal: https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html")
        print(f"2. Search for Europa observations with proposal ID {first_obs['proposal_id']}")
        print(f"3. Look for observation ID {first_obs['obs_id']}")
        print("4. Download the files and use astropy.io.fits to analyze them")
        
        # If we have G140L grating data (as shown in the observations)
        if 'G140L' in first_obs['filters']:
            print("\nThe G140L grating covers approximately 115-170 nm, which includes:")
            print("- The key Lyman-α (121.6 nm) and OI (130.4 nm) emission lines")
            print("- The spectral range where auroral emissions would be detected")
            print("- The UV region where sulfur compounds show distinctive absorption")

        def analyze_x1d_spectrum(fits_file, save_dir="plots"):
            """Analyze a single x1d spectrum file and create plots."""
            print(f"\nAnalyzing spectrum file: {os.path.basename(fits_file)}")
            
            # Create plots directory if it doesn't exist
            os.makedirs(save_dir, exist_ok=True)
            
            with fits.open(fits_file) as hdul:
                # Get observation details from header
                header = hdul[0].header
                obs_date = header.get('DATE-OBS', 'Unknown')
                exp_time = header.get('EXPTIME', 0)
                
                # Get spectrum data
                spec_data = hdul[1].data
                wavelength = np.array(spec_data['WAVELENGTH'][0])  # Convert to 1D array
                flux = np.array(spec_data['FLUX'][0])  # Convert to 1D array
                error = np.array(spec_data['ERROR'][0])  # Convert to 1D array
                
                # Convert wavelength to nm for better readability
                wavelength_nm = wavelength / 10.0
                
                # Create the plot
                plt.figure(figsize=(15, 10))
                
                # Plot spectrum with error region
                plt.plot(wavelength_nm, flux, 'b-', label='Flux', linewidth=1)
                plt.fill_between(wavelength_nm, flux-error, flux+error, 
                                color='blue', alpha=0.3, label='Error')
                
                # Add vertical lines for key emission features
                features = {
                    'H Lyman-α': 121.6,
                    'O I': 130.4,
                    'O I': 135.6,
                    'C II': 133.5,
                    'N V': 124.0
                }
                
                colors = ['r', 'g', 'c', 'm', 'y']
                for (feature, wav), color in zip(features.items(), colors):
                    plt.axvline(x=wav, color=color, linestyle='--', 
                               label=f'{feature} ({wav} nm)', alpha=0.7)
                
                # Customize plot
                plt.xlabel('Wavelength (nm)', fontsize=12)
                plt.ylabel('Flux (erg/s/cm²/Å)', fontsize=12)
                plt.title(f'HST STIS UV Spectrum of Europa\nDate: {obs_date}, Exposure: {exp_time}s', 
                         fontsize=14)
                plt.legend(fontsize=10)
                plt.grid(True, alpha=0.3)
                
                # Add text box with key statistics
                stats_text = f"""Observation Statistics:
Mean Flux: {np.mean(flux):.2e}
Max Flux: {np.max(flux):.2e}
Wavelength Range: {wavelength_nm.min():.1f} - {wavelength_nm.max():.1f} nm"""
                
                plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
                        fontsize=10, verticalalignment='top', 
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                # Save plot
                base_name = os.path.splitext(os.path.basename(fits_file))[0]
                plot_file = os.path.join(save_dir, f'{base_name}_spectrum.png')
                plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                plt.close()
                
                # Analyze emission features
                print("\nEmission Feature Analysis:")
                for feature, wav in features.items():
                    idx = np.abs(wavelength_nm - wav).argmin()
                    local_bg = np.median(flux[max(0, idx-5):min(len(flux), idx+6)])
                    feature_flux = flux[idx]
                    ratio = feature_flux / local_bg
                    
                    print(f"\n{feature} ({wav} nm):")
                    print(f"Flux: {feature_flux:.2e} erg/s/cm²/Å")
                    print(f"Local background: {local_bg:.2e} erg/s/cm²/Å")
                    print(f"Signal/Background ratio: {ratio:.2f}")
                    
                    if ratio > 1.2:
                        print("Status: EMISSION detected")
                    elif ratio < 0.8:
                        print("Status: ABSORPTION detected")
                    else:
                        print("Status: No significant feature detected")

        # Analyze all three x1d files
        data_dir = "data/hst/8224/calibrated/x1d"
        x1d_files = [
            "o5d601010_x1d.fits",  # First observation
            "o5d601080_x1d.fits",  # Middle observation
            "o5d6010a0_x1d.fits"   # Last observation
        ]

        print("Starting analysis of HST STIS UV spectra of Europa...")
        print("=" * 80)

        for fits_file in x1d_files:
            full_path = os.path.join(data_dir, fits_file)
            if os.path.exists(full_path):
                analyze_x1d_spectrum(full_path)
                print("=" * 80)
            else:
                print(f"Warning: File not found - {fits_file}")

        print("\nAnalysis complete. Check the 'plots' directory for spectrum visualizations.")

# If no STIS UV data is found, look for other instruments that might have UV data
if len(stis_uv) == 0:
    # Check for COS (Cosmic Origins Spectrograph) data
    cos_obs = df[df['instrument_name'].str.contains('COS', na=False)]
    if len(cos_obs) > 0:
        print(f"\nFound {len(cos_obs)} COS observations that might have UV data.")
        print("Consider modifying the script to analyze COS data instead.")
    
    # Check for ACS/SBC (Advanced Camera for Surveys / Solar Blind Channel) data
    acs_sbc_obs = df[(df['instrument_name'].str.contains('ACS', na=False)) & 
                     (df['instrument_name'].str.contains('SBC', na=False))]
    if len(acs_sbc_obs) > 0:
        print(f"\nFound {len(acs_sbc_obs)} ACS/SBC observations that might have UV imaging data.")
        print("Consider modifying the script to analyze ACS/SBC data instead.")