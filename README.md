# ptb-beta-paper
Code to implement simulations and data analyses from "A perturbation method for inference on regularized regression estimates"

The functions in this code can be used to obtain confidence intervals for regression coefficients from certain regularized regression methods. Please see the article for reference.

Minnier, J., Tian, L. & Cai, T. (2011), "A perturbation method for inference on regularized regression estimates"", Journal of the American Statistical Association 106(496), 1371–1382.

# Instructions

To replicate the simulations, run the files:
1-run-ResampleALASSO-simulations.R
2-run-CompressResample-simulations.R
3-run-Summary-simulations.R

for various values of `paramnum` to obtain different settings.

To replicate the data analysis, run the file:
ALASSO-HIVexample-May2012.R

# License

    Code to implement simulations and data analyses from 
    "A perturbation method for inference on regularized regression estimates"
    Copyright (C) 2010-2016  Jessica Minnier <minnier-at-ohsu.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
