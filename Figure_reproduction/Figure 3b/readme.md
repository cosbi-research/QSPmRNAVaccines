The present folder contains the script for the reproduction of panel b) of Figure 3 in the main text (BNT162b2 antibodies data fit).

## List of the files contained

### .mat FILES
- **BNT162b2_fit.mat**: Contains the parametrization obtained with BNT162b2 antibodies data for the general population.
- **GoelData.mat**: Contains the data from Goel et al.
- **KeshavarzData.mat**: Contains the data from Keshavarz et al.
- **Liang_data.mat**: Contains the data retrieved from Liang et al.
- **Liang_fit.mat**: Contains the parametrization obtained with Liang data.
- **SahinData.mat**: Contains the data from Sahin et al.

### SCRIPTS
- **Figure_3b.m**: Simulates the dynamics of the tissue layer at different doses and reproduces panel b) of Figure 3.

### FUNCTIONS
- **model_equations.m**: Contains all the equations of the tissue layer of the model.
- **parameters.m**: Contains all the parameters of the tissue layer of the model.
- **simulation.m**: The file that integrates the system of ODEs above.

### .txt FILES
- **readme.txt**: This file.

### .png FILES
- **Antibodies_in_blood.png**: Figure produced by "Figure_3b.m" script saved in PNG format.
