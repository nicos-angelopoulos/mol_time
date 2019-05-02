
R function for calculating the timing of DNA aberration events from 
NGS data.

Written by Peter Campbell and Francesco Maura (May 2017)
Tidy up/restructuring by Nicos  Angelopoulos, (2017/06/07)

This software is distributed under the MIT license:
https://opensource.org/licenses/MIT

To load simply do 
% R
> source("mol_time.R")

To run call the function with something like,
> mol_time(data.dir="data", res.dir="res-defaults")

To see if you getting the same results to ours, compare the 
outputs of the above call to the files in directory "test-defaults".

See documentation in man/mol_time.pdf

Authors: Peter Campbell, Francesco Maura & Nicos Angelopoulos
Maintainer & docs: Nicos Angelopoulos

CASM: Cancer Ageing and Somatic Mutations
Wellcome Sanger Institute, 
June, 2017 - Jan, 2018
---
http://www.sanger.ac.uk/science/programmes/cancer-genetics-and-genomics
