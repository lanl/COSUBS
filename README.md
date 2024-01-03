# COSUBS (COrrelated Sampling Using Batch Statistics)

COSUBS reads MCNP6 output files and uses tally fluctiation charts
(TFCs) to apply batch statistics to combinations of tallies.
The full descrition and user manual is published in 
Nuclear Science and Technology Open Research (NSTOR):

> Jeffrey A. Favorite, "A General-Purpose Code for Correlated Sampling Using
> Batch Statistics with MCNP6 for Fixed-Source Problems,"
> Nuclear Science and Technology Open Research, submitted (2023).

COSUBS has been built and tested only on a Linux platform with an Intel
FORTRAN compiler.

# Directory Contents

src--will contain source code when approved for open source by LANL

test--complete set of test inputs and outputs used in the NSTOR paper

NSTOR.extended_data.pdf--Appendices for the NSTOR paper
- Appendix A. Theoretical Variances for Analog Monte Carlo Calculations
- Appendix B. MCNP Base-Case Input File Listing
- Appendix C. Modifications to MCNP6.3 to Make the Maximum Number of TFC Entries a User Input

