# Fluorescence Analysis and Sigmoid Fitting

#### This document describes a Python-based code for the analysis of live cell puncta data analysis. The code processes raw output from ICY in Excel file format, align and clean the data as needed, calculates normalized Förster resonance energy transfer (nFRET) and color balance (CB), applies sigmoid curve fitting to identify kinetic info, and export visual and Excel summaries of results.  
There are 2 types of script included:
  1) .py - python script, usable in applications such as Spyder. 
  2) .ipynb  - python workbook, intended for use with Jupyter Notebook. More helpful if you want to know what each block of code is doing, or would like to easily output both the code + results in a pdf report. 
An alternate version of the code, which does ***not*** require the XY positional file is also included.

### Overview of Purpose and Functionality
The primary purpose of the script is to support the analysis of fluorescence (and, optionally, XY) data from _de novo_ particles over time. 
This code is meant for data with 3 fluorescence channels (e.g. YFP, FRET, and mCherry).

### Functions of this code:
1. Loads intensity and X-Y positional data from Excel files, and log files from Deltavision to automatically extract imaging parameters
2. Corrects any stitching or shifting in the raw excel caused by combining multiple tracking segments in ICY.
3. Pulls parameters (pixel to um conversion, temporal resolution, channel desinations) from Deltavision log files automatically to reduce user error.
4. Calculates normalized FRET (nFRET) and color balance (CB) values using user-defined bleedthrough correction coefficients.
5. Applies sigmoid curve fitting to estimate transition metrics and calculate ec10 (PR activation start time), ec90 (PR activation end time), ec90-ec10 (PR duration), along with mean squared error (MSE) as a measure of fitting effectiveness.
6. (OPTIONAL) Generate and export XY trajectory, fluorescence intensity, and fitted nFRET and CB plots.
7. Exports a summary Excel workbook summarizing kinetic results, calculated nFRET & CB values along with their sigmoid fits, and pngs of graphical outputs.

### Input Files and Requirements
This code can load 2 Excel files and 1 log file:
1. Intensity file that includes fluorescence info (.xls or .xlxs)
2. (OPTIONAL) The XY data file that contains x/y coordinates for each time point (.xls or .xlxs)
3. Deltavision log file (.dv)

Requirements:
+ If you don't have the XY file, use the version of the code that ends in NOXY
+ You also do **not** need to modify the raw files in any way.
+ This code is assuming that there is only **1** particle being represented in each file.
+ If you are running python for the first time, there may be some libraries that may need to be installed via the pip install input instead of just being loaded via the library command in-script. Whatever python application you use should tell you what is missing.

### Configuration Parameters
+ intensity_file: file path for the raw intensity file from ICY. Input in the form of r"Filepath"
+ XY_file: file path for the raw x-y position file from ICY. Input in the form of r"Filepath"
+ YFP_bleedthrough and mCH_bleedthrough: bleedthrough coefficients for YFP and mCherry, respectively.
      *The current bleedthrough values in the code is based on single fluorophore experiments done by Alex Z. in cells (values last updated 7/31/2025)
+ lower_bound_in_minutes and upper_bound_in_minutes: Optional time bounds (in minutes) restricting the lower and upper range of data used during sigmoid fitting, respectively. These allow users to focus
  the fit on a meaningful subset of the data. *To use the full range, set to **none**.*

### Data Processing and Calculating Additional Metrics
After loading, the data is cleaned and aligned using the *correct_shift()* function. This corrects for issues from stitching together multiple tracks in ICY. The three channels are then merged and renamed according to the log file.
*The code specifically looks for these channels: 488Laser/GFP, 488Laser/mCherry, 561Laser/mCherry, so data aquired using other filter sets  or fluororphores will **not** work without modification*

From this cleaned dataset, the following derived metrics are computed:
  1. Time (sec and min): Based on the frame index and temporal resolution.
  2. nFRET: A normalized FRET signal calculated by:
```math
$$nFRET = \frac{FRET - \text{YFP\_bleedthrough} * YFP - \text{mCH\_bleedthrough} * mCH}{mCH} $$
```
 \
  3. Color Balance (CB): The ratio of YFP to total fluorescence, calculated by:
```math
CB = \frac{YFP}{YFP + mCH}
```
 \
  4. x (um), y (um): Positional data in microns. Based on the Dim scaling factor.
  5. nFRET and CB sigmoid curve fits with residuals and MSE
     
### Curve Fitting and Kinetic Analysis
The script fits the following sigmoid function to nFRET and CB, either over the entire time range or over a user-defined window:
```math
f(x) = y_0+\frac{a}{1+e^{\frac{-(x-x_0)}{b}}}
```
 \
The resulting fit is then used to calculate:
  1. ec10: The time at which the curve reaches 10% of its amplitude.
  2. ec90: The time at which the curve reaches 90% of its amplitude.
  3. PR duration: The time interval between ec10 and ec90, representing the dynamic response window.
  4. MSE: Mean squared error of the fit within the fitting window. A value of 0 indicates a perfect fit, so we want the value to be as low as possible.
*In the event of poor fit within the set fitting range, the function is designed to return NaN values instead of erroring.*

### Output Files and Formats
Processed results are saved to the output directory under filenames automatically derived from the input base name. For instance, an input file named 2_3ROI1I.xls will generate:
  1. 2_3ROI1I_Summary.xlsx: An Excel workbook with 2 tabs:
      + All Values: The full dataset with all computed values and the sigmoid fits + residuals.
      + Kinetics: A summary table of fit parameters for both nFRET and CB.
  2. 2_3ROI1I_trajectory.png: A color-coded XY plot.
  3. 2_3ROI1I_intensity.png: A plot of raw fluorescence intensity (YFP and mCH) over time.
  4. 2_3ROI1I_nFRET_Fit.png: A comparison of measured vs. fitted nFRET values.
  5. 2_3ROI1I_CB_Fit.png: A comparison of measured vs. fitted color balance values.

All figures are saved in PNG format with 300 DPI resolution.

### Dependencies and Environment
The script requires Python 3 at a minimum.
