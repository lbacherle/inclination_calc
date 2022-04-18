&emsp;inclination_calc.jl

Simple Julia program to calculate launchsite parameters.

Modes:\
&emsp;1) Calculate inertial Azimuth (disregarding Earths rotation) for a launchsite or latitude and a desired inclination\
&emsp;2) Calculate rotational Azimuth (taking into account Earth rotation) for a launchsite or latitude and a desired inclination\
&emsp;3) Calculate resulting Inclination for a launchsite or latitude and an Azimuth angle\
&emsp;4) Compare propellant mass needed for different latitudes (assuming an eastwards launch)

Conventions:\
&emsp;1) All angles must be given in degrees\
&emsp;2) Inclination and Latitude are defined from -90° to +90°, Azimuth from 0° to 360°

Input:\
&emsp;If file called 'values.csv' in the correct format (see example values.csv) exists, reads the values from this csv and prints all output to files.\
&emsp;If this file does not exist, starts an interactive mode and allows user to freely input variables.

Dependencies:\
&emsp;1) The Plots package is needed to show the resulting plots. Installation can be done by typing julia into the command line to start the REPL, then typing the following:\
&emsp;&emsp;julia> using package\
&emsp;&emsp;julia> Pkg.add("Plots")

Sources:\
&emsp;1) Formulas: The derivation of formulas for inertial and rotational calculations can be found under https://www.orbiterwiki.org/index.php?title=Launch_Azimuth&oldid=17141
                 or https://ntrs.nasa.gov/api/citations/19980227091/downloads/19980227091.pdf ('Technical Note D-233, Determination of Azimuth Angle at Burnout for placing Satellite over a selected Earth position', T.H. Skopinski and K.G Johnson, NASA, 1960)\
&emsp;2) Launchsite Latitudes: 'Space Mission Engineering - The New SMAD' by Wertz, Everett and Puschell, published 2011 by Microcosm Press
