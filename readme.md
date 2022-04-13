    inclination_calc.jl

Simple Julia program to calculate launchsite parameters.

Modes:
    1) Calculate inertial Azimuth (disregarding Earths rotation) for a launchsite or latitude and a desired inclination
    2) Calculate rotational Azimuth (taking into account Earth rotation) for a launchsite or latitude and a desired inclination
    3) Calculate resulting Inclination for a launchsite or latitude and an Azimuth angle

Conventions:
\t1) All angles must be given in degrees
    2) Inclination and Latitude are defined from -90° to +90°, Azimuth from 0° to 360°

Sources:
    1) Formulas: The derivation of formulas for inertial and rotational calculations can be found under https://www.orbiterwiki.org/index.php?title=Launch_Azimuth&oldid=17141
                 or https://ntrs.nasa.gov/api/citations/19980227091/downloads/19980227091.pdf ('Technical Note D-233, Determination of Azimuth Angle at Burnout for placing Satellite over a selected Earth position', T.H. Skopinski and K.G Johnson, NASA, 1960)
    2) Launchsite Latitudes: 'Space Mission Engineering - The New SMAD' by Wertz, Everett and Puschell, published 2011 by Microcosm Press
