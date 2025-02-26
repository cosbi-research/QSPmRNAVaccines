SUPPLEMENTARY FILES S2

"A Multiscale Quantitative Systems Pharmacology model for the development and optimization of mRNA vaccines"

Lorenzo Dasti, Stefano Giampiccolo, Elisa Pettina’, Giada Fiandaca, Natascia Zangani, 
Lorena Leonardelli, Fabio De Lima Hedayioglu, Elio Campanile, Luca Marchetti
__________________________________________________________________________

The present folder contains the script for the reproduction of Figure7 in 
    the main text (mRNA-1273 antibodies data fit and validation).


List of the files contained -----------------------------------------------


.mat FILES:

- JochumData.mat: contains the data from Jochum et al. (validation);

- KeshavarzData.mat: contains the data from Keshavarz et al. (fit);

- KirsteData.mat: contains the data from Kirste et al. (validation);

- Liang_data.mat: contains the data retrieved from Liang et al.;

- Liang_fit.mat: contains the parametrization obtained with Liang data;

- mRNA-1273_fit.mat: contains the parametrization obtained with mRNA-1273
    antibodies data;


SCRIPTS:

- Figure_7.m: simulates the dynamics of the tissue layer at different 
    doses and reproduces all the panels of Figure7;


FUNCTIONS:

- logmeanandci.m: script to obtain the geometric mean and 95% confidence 
    interval from a vector of values that are supposed to have a lognormal 
    distribution;

- model_equations.m: contains all the equations of the tissue layer of the 
    model;

- parameters.m: contains all the parameters of the tissue layer of the 
    model;

- simulation.m: the file that integrates the system of ODEs above;


.txt FILES:

- readme.txt: this file.


.png FILES:

- Antibodies_in_blood.png: figure produced by "Figure_7.m" script saved in png format