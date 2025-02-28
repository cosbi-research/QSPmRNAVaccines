# SUPPLEMENTARY FILES S2

## "A Multiscale Quantitative Systems Pharmacology model for the development and optimization of mRNA vaccines"

### Authors
Lorenzo Dasti, Stefano Giampiccolo, Elisa Pettinà, Giada Fiandaca, Natascia Zangani,  
Lorena Leonardelli, Fabio De Lima Hedayioglu, Elio Campanile, Luca Marchetti

---

The present folder contains the script for the reproduction of panel a) of Figure 3 in the main text (antigen-presenting cells data fit).

## List of the files contained

### .mat FILES
- **Liang_data.mat**: Contains the data retrieved from Liang et al.
- **Liang_fit.mat**: Contains the parametrization obtained with Liang data.

### SCRIPTS
- **Figure_3a.m**: Simulates the dynamics of the antigen-presenting cells part of the tissue layer and reproduces panel a) of Figure 3.

### FUNCTIONS
- **model_equations.m**: Contains all the equations of the APCs part of the tissue layer (both in the injection site and in the lymph node).
- **parameters.m**: Contains all the parameters of the APCs part of the tissue layer.
- **simulation.m**: The file that integrates the system of ODEs above.

### .txt FILES
- **readme.txt**: This file.

### .png FILES
- **Liang_Data.png**: Figure produced by "Figure_3a.m" script saved in PNG format.
