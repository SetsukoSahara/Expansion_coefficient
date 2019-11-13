# A common rule governing differentiation kinetics of mouse cortical progenitors

Setsuko Sahara, Takashi Kodama, and Charles F Stevens

Here we provide two R-scripts used in our paper to model the differentiation kinetics of mouse cortical progenitors.
"script_common_Tc.R" is to estimate the expansion coefficients by using the total cell cycle length of all progenitors, and "script_different_Tc.R" incorporates cell cycle lengths unique to each progenitor type.

For details please see our paper.
--here we need to cite the paper after the acceptane--


## Requirement
R 3.3.3 or above

## Usage
1. Make sure two csv files (cell_fraction.csv, Tc.csv) and these script files are all in the current working directory,
   In console, run "source("script_common_Tc.R")"
   In console, run "expansion.coefficients()"
2. In console, run “source("script_common_Tc.R")" and “source("script_different_Tc.R")”
3. To replicate the modeled cell fraction (Fig. 1U), the cell cycle length of all progenitors (Fig. 2L), and the expansion coefficients estimated based on the cell cycle length of all progenitors, run “expansion.coefficients()” in console.
4. To replicate the expansion coefficients estimated based on the cell cycle lengths unique to each progenitor type, run “expansion.coefficients2()” in console.


