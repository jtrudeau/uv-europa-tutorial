# Student Guide: Analyzing UV Spectra of Europa

## Introduction to Europa Science

Europa is one of Jupiter's most fascinating moons and a prime target for astrobiology research due to its subsurface ocean. With a diameter of 3,121 kilometers, it's slightly smaller than Earth's Moon but contains more water than all of Earth's oceans combined. Europa's icy surface is relatively young (40-90 million years) and marked by a complex network of cracks and ridges, suggesting dynamic geological processes.

### Why Study Europa in the UV?

Ultraviolet observations reveal several key features of Europa:

1. **Atmospheric Composition**: UV spectra can detect emissions from Europa's tenuous atmosphere, primarily oxygen and hydrogen.
2. **Plume Activity**: Water vapor plumes, similar to those on Saturn's moon Enceladus, may be detectable in the UV.
3. **Magnetospheric Interactions**: Europa orbits within Jupiter's intense magnetic field, causing aurora-like emissions.
4. **Surface Composition**: UV reflectance spectra provide information about Europa's surface materials.

## The MAST Archive: Your Gateway to HST Data

### What is MAST?

The Mikulski Archive for Space Telescopes (MAST) is NASA's primary repository for astronomical data from space-based observatories. MAST is hosted by the Space Telescope Science Institute (STScI) and provides access to data from:

- Hubble Space Telescope (HST)
- James Webb Space Telescope (JWST)
- Transiting Exoplanet Survey Satellite (TESS)
- Kepler/K2
- Many other missions

### How to Access MAST Data

1. **Portal Interface**: Visit the [MAST Portal](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html) for a user-friendly search interface.
2. **API Access**: Use `astroquery.mast` in Python to programmatically access data. See the [astroquery documentation](https://astroquery.readthedocs.io/en/latest/mast/mast.html) and [MAST API documentation](https://mast.stsci.edu/api/v0/) for detailed guides.
3. **Direct Download**: Use direct URLs for specific datasets.

### Step-by-Step Guide to Using the MAST Portal

1. **Basic Search**:
   - Go to the MAST Portal homepage
   - In the search box, enter "Europa"
   - Under "Mission," select "HST" to filter for Hubble observations
   - Click "Search" to see all HST observations of Europa

2. **Finding Specific Observations**:
   - Look for Proposal ID "8224" in the results (our tutorial data is from this project)
   - Click on an observation to see its details
   - The "Preview" tab shows quick-look images or spectra
   - The "Files" tab lists all available data products

3. **Filtering Data**:
   - You can use the left sidebar to filter by:
     * Instrument (e.g., "STIS")
     * Data type (e.g., "SPECTROSCOPIC")
     * Observation date
     * Wavelength range

4. **Downloading Data**:
   - Select files by checking the boxes next to them
   - Click "Download" to get selected files
   - Choose between direct download or a download script
   - For large datasets, use the download script option
  - ***For the tutorial data has already been provided for you***

5. **Advanced Search Tips**:
   - Use the "Advanced Search" form for complex queries
   - Combine multiple filters (e.g., date range AND instrument)
   - Search by coordinates using the "Resolve" button
   - Save your searches for future use

6. **Working with Data Products**:
   - Look for files ending in:
     * _raw.fits: Raw telescope data
     * _flt.fits: Calibrated, flat-fielded images
     * _x1d.fits: Extracted 1D spectra
     [**These are the ones we work with in the tutorial]**
     * _x2d.fits: 2D rectified spectra
   - Download associated calibration files if needed

### Understanding FITS Files

FITS (Flexible Image Transport System) is the standard file format for astronomical data:

- **Headers**: Contain metadata about the observation (date, exposure time, coordinates)
- **Data Extensions**: Contain the actual data (images, spectra, tables)
- **Multiple HDUs**: A single FITS file can contain multiple Headers+Data Units

HST data typically includes multiple file types:
- **_raw.fits**: Raw, uncalibrated data
- **_flt.fits**: Flat-fielded, calibrated data
- **_x1d.fits**: Extracted 1D spectra
- **_x2d.fits**: Rectified 2D spectra

## HST Proposal 8224: Exploring Europa's Aurora

Our tutorial focuses on data from HST Proposal 8224, "UV Imaging of Europa & Ganymede: Unveiling Satellite Aurora & Electrodynamical Interactions," led by Principal Investigator Melissa A. McGrath (NASA/MSFC).

### Proposal Overview

This 1999 study aimed to:
1. Detect and characterize auroral emissions around Europa
2. Study the interaction between Europa and Jupiter's magnetosphere
3. Compare Europa's emissions with those of Ganymede
4. Investigate the tenuous atmospheres of these icy moons

### Observation Strategy

The observations were conducted using the STIS (Space Telescope Imaging Spectrograph) instrument with the FUV-MAMA detector and G140L grating, covering wavelengths from 114.0-173.0 nm. This spectral range includes:

- Hydrogen Lyman-α emission (121.6nm)
- Oxygen emissions (130.4nm, 135.6nm)
- Various other emission lines that could indicate atmospheric composition

The observations were taken on October 5, 1999, spanning approximately 6.5 hours (8:39 UT to 15:14 UT), allowing the team to observe Europa at different orbital positions around Jupiter.

### Scientific Context

These observations were groundbreaking because:
1. They provided some of the first UV spectra of Europa taken with high sensitivity
2. They helped establish the presence of an oxygen atmosphere around Europa
3. They studied time variability that could be related to plume activity
4. They improved our understanding of moon-magnetosphere interactions

## Time Information in HST Data

Astronomical observations typically use Modified Julian Date (MJD) for timestamps:
- MJD = JD - 2400000.5 (where JD is Julian Date)
- The integer part represents the day
- The decimal part represents the fraction of a day

For HST Proposal 8224:
- MJD 51456 corresponds to October 5, 1999
- The observations span approximately MJD 51456.36 to 51456.64
- This represents Europa at different orbital positions around Jupiter

## Analyzing the Data

In our tutorial, we're working with:
1. **x1d files**: 1D extracted spectra showing flux vs. wavelength
2. **Header information**: Contains critical metadata about the observations
3. **Multiple observations**: Allow us to study time variability

Key features to look for in the spectra:
- **Emission lines**: Peaks at specific wavelengths that indicate atomic/molecular emissions
- **Temporal variations**: Changes in emission strength over time
- **Correlation with orbital position**: How emissions vary with Europa's position around Jupiter

## References and Further Reading

1. McGrath, M. A., et al. (2004). "Satellite atmospheres." Jupiter: The Planet, Satellites and Magnetosphere, 1, 457-483. [link.](https://lasp.colorado.edu/mop/files/2015/08/jupiter_ch19-1.pdf)

2. Hall, D. T., et al. (1995). "Detection of an oxygen atmosphere on Jupiter's moon Europa." Nature, 373(6516), 677-679. [link.](https://www.nature.com/articles/373677a0)

3. Roth, L., et al. (2014). "Transient water vapor at Europa's south pole." Science, 343(6167), 171-174. [link.](https://pubmed.ncbi.nlm.nih.gov/24336567/)

4. HST Proposal Information: [MAST Holdings for Proposal 8224](https://archive.stsci.edu/proposal_search.php?id=8224&mission=hst)

5. Space Telescope Science Institute Documentation: [STIS Data Handbook](https://hst-docs.stsci.edu/stisdhb)

6. NASA's Europa Mission: [Europa Clipper](https://europa.nasa.gov/)

## Connection to Current Research

The data you're analyzing in this tutorial is connected to ongoing research about Europa and several upcoming missions that will revolutionize our understanding of this fascinating moon.

### Europa Clipper Mission (NASA)

**Launch**: October 2024
**Arrival**: 2030

Key Science Goals:
1. **Ice Shell & Ocean**:
   - Confirm the presence of subsurface water
   - Determine ice shell thickness
   - Map potential landing sites for future missions

2. **Composition**:
   - Analyze surface chemistry
   - Search for organic compounds
   - Study the composition of any active plumes

3. **Current Activity**:
   - Monitor plume activity
   - Study cryovolcanic features
   - Observe surface changes

Instruments relevant to UV studies:
- **Europa-UVS**: Ultraviolet Spectrograph
  * Will observe in the 55-210 nm range
  * Higher resolution than HST STIS
  * Can detect plume activity
  * Will study atmospheric composition

### JUICE Mission (ESA)

**Launch**: April 2023
**Arrival**: July 2031

Key Features:
1. **Multi-Moon Focus**:
   - Will study Europa, Ganymede, and Callisto
   - First detailed comparison of these icy worlds

2. **Europa Flybys**:
   - Two dedicated Europa flybys
   - High-resolution imaging
   - Subsurface sounding

3. **Complementary Science**:
   - Will observe Europa in different wavelengths
   - Can compare with Clipper observations
   - Studies Jupiter's magnetosphere effects

### Future Mission Concepts

1. **Europa Lander** (Under Study):
   - Would directly sample surface material
   - Search for biosignatures
   - Analyze surface composition in detail

2. **Potential Network Missions**:
   - Multiple small landers or penetrators
   - Seismic studies of ice shell
   - Monitor geological activity

### What's Next

The skills you're learning here will be valuable for:
1. **Data Analysis**:
   - Similar techniques apply to new mission data
   - Understanding spectral features
   - Time series analysis

2. **Scientific Context**:
   - Historical perspective on Europa (satellites in general) studies
   - Understanding the basics of observation techniques
   - Importance of long-term monitoring. Mission data we could prep for in this project won't be available for years yet!

3. **Possible Next Steps**:
   - Scanning open problems with existing (public) data sets 
   - Exploring potential for comparative studies
   - Integration with other observations
   - Extended analysis of previous work
  
  ## **Task** 
  Continue to familiarize yourself with the MAST observation platform and work through `europa_uv_tutorial.py`