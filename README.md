# COSUBS (COrrelated Sampling Using Batch Statistics)

COSUBS reads MCNP6 output files and uses tally fluctuation charts
(TFCs) to apply batch statistics to combinations of tallies.
The full description and user manual is published in 
Nuclear Science and Technology Open Research (NSTOR):

> Jeffrey A. Favorite, "A General-Purpose Code for Correlated Sampling Using
> Batch Statistics with MCNP6 for Fixed-Source Problems,"
> Nuclear Science and Technology Open Research, awaiting open peer review (2024);
> https://doi.org/10.12688/nuclscitechnolopenres.17496.1.

COSUBS has been built and tested only on a Linux platform with an Intel
FORTRAN compiler.

# Directory Contents

src--source code

test--complete set of test inputs and outputs used in the NSTOR paper

NSTOR.extended_data.pdf--Appendices for the NSTOR paper
- Appendix A. Theoretical Variances for Analog Monte Carlo Calculations
- Appendix B. MCNP Base-Case Input File Listing
- Appendix C. Modifications to MCNP6.3 to Make the Maximum Number of TFC Entries a User Input

# Copyright

© 2023. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so. 

# License

This program is Open-Source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

