A Multiscale Quantitative Systems Pharmacology Model for the Development and Optimization
of mRNA Vaccines

Lorenzo Dasti, Stefano Giampiccolo, Elisa Pettinà, Giada Fiandaca, Natascia Zangani, 
Lorena Leonardelli, Fabio De Lima Hedayioglu, Elio Campanile, Luca Marchetti

SUPPLEMENTARY FILE S2

This folder contains MATLAB scripts for simulating the tissue layer, as described in Step 
3 of Supplementary Section S5.

Contents

.mat Files:
- BNT162b2_fit.mat: Parameter set based on antibody data for the BNT162b2 vaccine 
  in the general population.
- BNT162b2_fit_over60.mat: Parameter set based on antibody data for the BNT162b2 
  vaccine in individuals over 60 years of age.
- Liang_data.mat: Dataset sourced from Liang et al.
- Liang_fit.mat: Parameter set obtained using data from Liang et al.
- mRNA-1273_fit.mat: Parameter set based on antibody data for the mRNA-1273 vaccine.

Scripts:
- main.m: Simulates tissue layer dynamics, using antigen-presenting cell (APC) parameters 
  and either the BNT162b2 vaccine (general or over 60 y.o. populations) or the mRNA-1273 
  vaccine parameters.

Functions:
- model_equations.m: Defines all equations governing the tissue layer.
- parameters.m: Contains parameter values for the tissue layer.
- simulation.m: Integrates the system of ODEs in the tissue layer model.

Text Files:
- README.txt: This file, which provides an overview of folder contents.