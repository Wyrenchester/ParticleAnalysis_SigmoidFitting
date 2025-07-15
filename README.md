# Fluorescence Analysis and Sigmoid Fitting

#### This document describes a Python-based code for the analysis of live cell _de novo_ puncta data analysis. The code processes raw output from ICY in Excel file format, align and clean the data as needed, calculates normalized Förster resonance energy transfer (nFRET) and color balance (CB), applies sigmoid curve fitting to identify kinetic info, and export visual and Excel summaries of results. The script is intended for use with Jupyter Notebook and meant for data with 3 fluorescence channels (e.g. YFP, FRET, and mCherry).

### Overview of Purpose and Functionality
The primary purpose of the script is to support the analysis of fluorescence (and, optionally, XY) data from _de novo_ particles over time. 

#### Functions of this code:
1. Loads intensity and positional data from Excel files.
2. Corrects any stitching artifacts caused by combining multiple tracking segments in ICY.
3. Aligns and renames intensity data based on user-specified channel designations.
4. Calculates normalized FRET (nFRET) and color balance (CB) values using user-defined bleedthrough correction coefficients.
5. Applies sigmoid curve fitting to estimate transition metrics and calculate ec10, ec90, along with mean squared error as a measure of fitting effectiveness.
6. Generate and export XY trajectory, fluorescence intensity, and fitted nFRET and CB plots.
7. Exports a summary Excel workbook containing both raw and calculated/fitted data as well as kinetic info.

#### Input Files and Requirements
This code can load two Excel files:
1. The intensity file that includes fluorescence info. The script assumes that sheet indices 0, 6, and 12 are the channels of interest, which is mean intensity. Note that in python, numbering starts at 0 (e.g. the first tab is 0, and the 2nd tab is 1).
2. (OPTIONAL) The XY data file that contains x/y coordinates for each time point.

Requirements:
If you don't have the XY file, use the version of the code that ends in NOXY
You also do **not** need to modify the raw files in any way.
This code is assuming that there is only **1** particle being represented in each file.

Users must specify the file paths in the script under the variables *filepath1* and *filepath2*. Output files are saved to a designated directory defined by *filepath_out*.

#### Configuration Parameters
+ YFP, FRET, mCH: integers indicating which of the three channels corresponds to YFP, FRET, and mCherry. Numbering starts at 1, **not** 0. Should range from 1-3.
+ YFP_BT and mCH_BT: bleedthrough coefficients for YFP and mCherry.
+ Tempres: temporal resolution of imaging (in seconds between frames).
+ Dim: scaling factor to convert pixel units into microns.
+ fit_x_min and fit_x_max: Optional time bounds (in minutes) restricting the range of data used during sigmoid fitting. These allow users to focus the fit on a meaningful subset of the data. To use the full range, set to *none*.

#### Data Processing and Calculating Additional Metrics
After loading, the data is cleaned and aligned using the *correct_shift()* function. This corrects for issues from stitching together multiple tracks in ICY. The three channels are then merged and renamed according to user input.

From this cleaned dataset, the following derived metrics are computed:
  1. Time (sec and min): Based on the frame index and temporal resolution.
  2. nFRET: A normalized FRET signal calculated by:

          (FRET - YFP_BT * YFP - mCH_BT * mCH) / mCH
  4. Color Balance (CB): The ratio of YFP to total fluorescence, calculated by:

          YFP / (YFP + mCH)
  6. x (um), y (um): Positional data in microns. Based on the Dim scaling factor.

#### Curve Fitting and Kinetic Analysis
The script fits a sigmoid function with the formula:

          y = L / (1 + exp(-k(x - x0))) + b
      
to nFRET and CB, either over the entire time range or over a user-defined window.

The resulting fit is then used to calculate:
  1. ec10: The time at which the curve reaches 10% of its amplitude.
  2. ec90: The time at which the curve reaches 90% of its amplitude.
  3. PR duration: The time interval between ec10 and ec90, representing the dynamic response window.
  4. MSE: Mean squared error of the fit within the fitting window. A value of 0 indicates a perfect fit, so we want the value to be as low as possible.

*In the event of poor fit within the set fitting range, the function is designed to return NaN values instead of erroring.*

#### Output Files and Formats
Processed results are saved to the output directory under filenames automatically derived from the input base name. For instance, an input file named 2_3ROI1I.xls will generate:
  1. 2_3ROI1I_Summary.xlsx: An Excel workbook with 2 tabs:
      + All Values: The full dataset with all computed values and the sigmoid fits + residuals.
      + Kinetics: A summary table of fit parameters for both nFRET and CB.
  2. 2_3ROI1I_trajectory.png: A color-coded XY plot.
  3. 2_3ROI1I_intensity.png: A plot of raw fluorescence intensity (YFP and mCH) over time.
  4. 2_3ROI1I_nFRET_Fit.png: A comparison of measured vs. fitted nFRET values.
  5. 2_3ROI1I_CB_Fit.png: A comparison of measured vs. fitted color balance values.

All figures are saved in PNG format with 300 DPI resolution.

#### Dependencies and Environment
The script requires Python 3 and the following packages:
+ pandas
+ numpy
+ matplotlib
+ scipy
+ sklearn
+ xlrd
+ openpyxl
