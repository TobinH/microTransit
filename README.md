# microtransit
qPLM tools in R

microtransit provides a set R language tools for importing and analyzing quantitative polarized light microscopy (qPLM) data. Current task sets include:

1. read intensity, retardance, and azimuth bitmap images from Werner Kaminsky's Rotopol software into R, with calibration options for specimen thickness, illumination wavelength, and expected birefringence.
2. combine retardance and azimuth data into a single false-color output image using the CIELUV perceptual color space.
3. compute summary spherical statistics of slow-axis orientations per sample.
4. spatial eigenfunction analysis of slow-axis orientation.

task sets in development include:

5. Jones matrix modeling and qPLM simulation tools for simulating expected qPLM data distributions in idealized tissue arrangements.
6. manual qPLM input from multiple images at known analyzer orientations.

