SUPPLEMENTARY FILES S2

"A Multiscale Quantitative Systems Pharmacology model for the development and optimization of mRNA vaccines"

Lorenzo Dasti, Stefano Giampiccolo, Elisa Pettina’, Giada Fiandaca, Natascia Zangani, 
Lorena Leonardelli, Fabio De Lima Hedayioglu, Elio Campanile, Luca Marchetti
__________________________________________________________________________

The present folder contains the script for the reproduction of the panel b) 
    of Figure3 in the main text (BNT162b2 antibodies data fit).


List of the files contained -----------------------------------------------


.mat FILES:

- BNT162b2_fit.mat: contains the parametrization obtained with BNT162b2
    antibodies data for the general population;

- GoelData.mat: contains the data from Goel et al.;

- KeshavarzData.mat: contains the data from Keshavarz et al.;

- Liang_data.mat: contains the data retrieved from Liang et al.;

- Liang_fit.mat: contains the parametrization obtained with Liang data;

- SahinData.mat: contains the data from Sahin et al.;


SCRIPTS:

- Figure_3b.m: simulates the dynamics of the tissue layer at different 
    doses and reproduces panel b) of Figure3;


FUNCTIONS:

- model_equations.m: contains all the equations of the tissue layer of the 
    model;

- parameters.m: contains all the parameters of the tissue layer of the 
    model;

- simulation.m: the file that integrates the system of ODEs above;


.txt FILES:

- readme.txt: this file.

.png FILES:

- Antibodies_in_blood.png: figure produced by "Figure_3b.m" script saved in png format
