The present folder contains the script for the reproduction of Figure 7 in the main text (mRNA-1273 antibodies data fit and validation).

## List of the files contained

### .mat FILES
- **JochumData.mat**: Contains the data from Jochum et al. (validation).
- **KeshavarzData.mat**: Contains the data from Keshavarz et al. (fit).
- **KirsteData.mat**: Contains the data from Kirste et al. (validation).
- **Liang_data.mat**: Contains the data retrieved from Liang et al.
- **Liang_fit.mat**: Contains the parametrization obtained with Liang data.
- **mRNA-1273_fit.mat**: Contains the parametrization obtained with mRNA-1273 antibodies data.

### SCRIPTS
- **Figure_7.m**: Simulates the dynamics of the tissue layer at different doses and reproduces all the panels of Figure 7.

### FUNCTIONS
- **logmeanandci.m**: Script to obtain the geometric mean and 95% confidence interval from a vector of values that are supposed to have a lognormal distribution.
- **model_equations.m**: Contains all the equations of the tissue layer of the model.
- **parameters.m**: Contains all the parameters of the tissue layer of the model.
- **simulation.m**: The file that integrates the system of ODEs above.

### .txt FILES
- **readme.txt**: This file.

### .png FILES
- **Antibodies_in_blood.png**: Figure produced by "Figure_7.m" script saved in PNG format.
