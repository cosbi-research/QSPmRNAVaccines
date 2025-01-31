SUPPLEMENTARY FILES S2

"A Multiscale Quantitative Systems Pharmacology model for the development and optimization of mRNA vaccines"

Lorenzo Dasti, Stefano Giampiccolo, Elisa Pettina’, Giada Fiandaca, Natascia Zangani, 
Lorena Leonardelli, Fabio De Lima Hedayioglu, Elio Campanile, Luca Marchetti
__________________________________________________________________________

The present folder contains the script for the reproduction of the panel a) 
    of Figure3 in the main text (antigen presenting cells data fit).


List of the files contained -----------------------------------------------


.mat FILES:

- Liang_data.mat: contains the data retrieved from Liang et al.;

- Liang_fit.mat: contains the parametrization obtained with Liang data;


SCRIPTS:

- Figure_3a.m: simulates the dynamics of the antigen presenting cells part 
    of tissue layer and reproduces panel a) of Figure3;


FUNCTIONS:

- model_equations.m: contains all the equations of the APCs part of the 
    tissue layer (both in the injection site and in lymph node);

- parameters.m: contains all the parameters of the APCs part of the tissue 
    layer;

- simulation.m: the file that integrates the system of ODEs above;


.txt FILES:

- readme.txt: this file.


.png FILES:

- Liang_Data.png: figure produced by "Figure_3a.m" script saved in png format