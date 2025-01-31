SUPPLEMENTARY FILES S2

"A Multiscale Quantitative Systems Pharmacology model for the development and optimization of mRNA vaccines"

Lorenzo Dasti, Stefano Giampiccolo, Elisa Pettina’, Giada Fiandaca, Natascia Zangani, 
Lorena Leonardelli, Fabio De Lima Hedayioglu, Elio Campanile, Luca Marchetti
__________________________________________________________________________

The present folder contains the script for the reproduction of Figure4 in 
    the main text (BNT162b2 antibodies data validation).


List of the files contained -----------------------------------------------


.mat FILES:

- BNT162b2_fit.mat: contains the parametrization obtained with BNT162b2
    antibodies data for the general population;

- Liang_data.mat: contains the data retrieved from Liang et al.;

- Liang_fit.mat: contains the parametrization obtained with Liang data;

- NaaberData.mat: contains the data from Naaber et al.;

- SahinData.mat: contains the data from Sahin et al.;

- SahinData.mat: contains the data from Takeuchi et al.;


SCRIPTS:

- Figure_4.m: simulates the dynamics of the tissue layer at different 
    doses and reproduces Figure4;


FUNCTIONS:

- model_equations.m: contains all the equations of the tissue layer of the 
    model;

- parameters.m: contains all the parameters of the tissue layer of the 
    model;

- simulation.m: the file that integrates the system of ODEs above;

- simulation_tris.m: the file that integrates the system of ODEs above 
    when a third dose is administered;


.txt FILES:

- readme.txt: this file.

.png FILES:

- Antibodies_in_blood.png: figure produced by "Figure_4.m" script saved in png format