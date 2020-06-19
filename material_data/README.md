# Fitting Data
The process for fitting experimental index or permittivity data is (unfortunately)
not well streamlined as of yet.

## Generate a Data File
To fit your own data using the provided scripts, first make a space-delimited
file in your local copy of this directory with the data in the form:
(wavelength, um) (n) (k) (like the provided examples).

## Using Levy Estimation for Initial Fitting (optional)
Open up Matlab, if you intend to use Levy estimation to generate initial guesses (which can sometimes
improve the fit results, but not always), and run `initial_fits_v1.m`. Change the
input file to your newly-created data file, and potentially adjust `num_poles`.
`num_poles` should be `2n + 1` where `n` is the number of poles you think you
need to accurately fit your data.


When you run `initial_fits_v1`, Matlab
should print a list of poles and residues (`poles` and `res`). The poles are the
initial guesses for our `a` values while the residues are for `c`. These values
are then used in the next step.

## Fitting using Mathematica
Now that you have a data file (and optionally initial estimates for the parameters),
open `material_fitting.nb`. Analogously to the demo data already there, add lines
to read in your data (`Import[...]`), trim it to the wavelength range of interest
(`Select[...]`), and convert it to epsilon as a function of eV (`epsFunc[...]`).
Run the code sections that define the fitting, plotting, and output functions, as
well as the shared constraint definition.
Finally, copy one of the other material blocks and then update `n`, `data`, and
the estimate lists (`aReEst, aImEst,...`) to use your data. When you run that
section, it should produce a (hopefully good!) fit for the dielectric function of your data.

The produced values can then either be added to the material header file (`include/fdtd/material_data.h`)
or be used with the custom material entry on the website. Note that you might need
to keep an eye on the simulation the first time you generate a new fit as it is
possible that it will make the simulation unstable (this will be very obvious).
