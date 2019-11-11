## A common rule governing differentiation kinetics of mouse cortical progenitors

Setsuko Sahara, Takashi Kodama, and Charles F Stevens

Here We provide two R-scripts used in our paper to model the differentiation kinetics of mouse cortical progenitors.
"script_common_Tc.R" is to etimate the expansion cofficients by using the total cell cycle length of all progenitors, and "script_different_Tc.R" is a modified script by using cell cycle length unique to an each progenitor type.

For details please see our paper.
--here we need to cite the paper after the acceptane--


# Requirement
R 3.3.3 or above

# Usage
1. Make sure two csv files (cell_fraction.csv, Tc.csv) and this script file are all in the current working directory,
In console, run "source("script_common_Tc.R")"
In console, run "expansion.coefficients()"


