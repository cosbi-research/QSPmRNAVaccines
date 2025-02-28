The present folder contains the script for the reproduction of Figure 4 in the main text (BNT162b2 antibodies data validation).

## List of the files contained

### .mat FILES
- **BNT162b2_fit.mat**: Contains the parametrization obtained with BNT162b2 antibodies data for the general population.
- **Liang_data.mat**: Contains the data retrieved from Liang et al.
- **Liang_fit.mat**: Contains the parametrization obtained with Liang data.
- **NaaberData.mat**: Contains the data from Naaber et al.
- **PayneData.mat**: Contains the data from Payne et al.
- **SahinData.mat**: Contains the data from Sahin et al.
- **TakeuchiData.mat**: Contains the data from Takeuchi et al.

### SCRIPTS
- **Figure_4.m**: Simulates the dynamics of the tissue layer at different doses and reproduces Figure 4.

### FUNCTIONS
- **model_equations.m**: Contains all the equations of the tissue layer of the model.
- **parameters.m**: Contains all the parameters of the tissue layer of the model.
- **simulation.m**: The file that integrates the system of ODEs above.
- **simulation_tris.m**: The file that integrates the system of ODEs above when a third dose is administered.

### .txt FILES
- **readme.txt**: This file.

### .png FILES
- **Antibodies_in_blood.png**: Figure produced by "Figure_4.m" script saved in PNG format.
