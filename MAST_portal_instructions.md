# MAST Portal Instructions: A Visual Guide

## Overview

The Mikulski Archive for Space Telescopes (MAST) is NASA's primary repository for astronomical data from space-based observatories. This guide will help you navigate the MAST Portal to find and download HST observations.

## Access Methods

### 1. Web Portal
The [MAST Portal](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html) provides a user-friendly interface for searching and downloading data. This guide focuses primarily on the web portal interface.

### 2. Programmatic Access
For programmatic access, MAST provides several options:

1. **Astroquery MAST Module**:
   - Documentation: [astroquery.mast documentation](https://astroquery.readthedocs.io/en/latest/mast/mast.html)
   - Python package for accessing astronomical data
   - Includes examples and tutorials
   - Part of the larger astroquery package

2. **MAST API**:
   - Documentation: [MAST API documentation](https://mast.stsci.edu/api/v0/)
   - RESTful web service
   - Supports multiple programming languages
   - Includes authentication for proprietary data

3. **Jupyter Notebooks**:
   - [MAST API Notebooks](https://github.com/spacetelescope/notebooks/tree/master/notebooks/MAST)
   - Example workflows and tutorials
   - Ready-to-run code samples

## Step-by-Step Guide

### 1. Basic Search
- Go to the MAST Portal homepage: [mast.stsci.edu/portal](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html)
- In the search box, enter "Europa"
- Under "Mission," select "HST" to filter for Hubble observations
- Click "Search" to see all HST observations of Europa

### 2. Finding Specific Observations
- Look for Proposal ID "8224" in the results
- Click on an observation to see its details
- The "Preview" tab shows quick-look images or spectra
- The "Files" tab lists all available data products

### 3. Filtering Data
Use the left sidebar to filter by:
- Instrument (e.g., "STIS")
- Data type (e.g., "SPECTROSCOPIC")
- Observation date
- Wavelength range

### 4. Downloading Data
- Select files by checking the boxes next to them
- Click "Download" to get selected files
- Choose between direct download or a download script
- For large datasets, use the download script option

### 5. Advanced Search Tips
- Use the "Advanced Search" form for complex queries
- Combine multiple filters (e.g., date range AND instrument)
- Search by coordinates using the "Resolve" button
- Save your searches for future use

### 6. Working with Data Products
Look for files ending in:
- _raw.fits: Raw telescope data
- _flt.fits: Calibrated, flat-fielded images
- _x1d.fits: Extracted 1D spectra
- _x2d.fits: 2D rectified spectra

Remember to download associated calibration files if needed.

## Visual Guide

Here are key screenshots to help you navigate the MAST Portal:

[Screenshot 1: MAST Portal Home]
- The main search interface at mast.stsci.edu
- Shows the search box where you enter "Europa"
- Highlights the "Mission" filter dropdown
- Note: The search box is in the center of the page

[Screenshot 2: Search Results]
- Results page showing HST observations of Europa
- Left sidebar shows available filters
- Each row represents one observation
- Note the columns for Proposal ID, Instrument, and Observation Date

[Screenshot 3: Observation Details]
- Detailed view of an observation from Proposal 8224
- Shows the "Preview" and "Files" tabs
- Displays observation metadata
- Note the available data products

[Screenshot 4: Data Product Selection]
- The file selection interface
- Shows different types of FITS files
- Checkboxes for selecting files
- Download options at the bottom

[Screenshot 5: Advanced Search Form]
- Shows the advanced search interface
- Multiple filter options
- Coordinate search tools
- Save search functionality

Note: Screenshots should be added to this guide to provide visual reference. The actual interface may change slightly over time, but the basic layout and functionality remain similar.

## Tips for Efficient Use

1. **Saving Time**:
   - Use filters early to narrow down results
   - Save common searches for quick access
   - Use the download script for multiple files

2. **Finding the Right Data**:
   - Check the observation details carefully
   - Look at previews before downloading
   - Note the data quality flags

3. **Managing Downloads**:
   - Create organized directories for different datasets
   - Keep track of calibration files
   - Document your search parameters

## Common Issues and Solutions

1. **Too Many Results**:
   - Add more specific filters
   - Use the Advanced Search
   - Filter by observation date

2. **Missing Files**:
   - Check if proprietary period has ended
   - Look for alternative datasets
   - Contact MAST Help if needed

3. **Download Problems**:
   - Try the download script instead of direct download
   - Check your internet connection
   - Clear browser cache if needed

## Getting Help

- MAST Help Desk: [archive@stsci.edu](mailto:archive@stsci.edu)
- Documentation: [MAST Documentation](https://archive.stsci.edu/docs/)
- FAQ: [MAST FAQ](https://archive.stsci.edu/docs/faqs/) 