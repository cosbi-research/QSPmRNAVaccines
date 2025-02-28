The present folder contains the script for the reproduction of Figure 6 in the main text (BNT162b2 antibodies data fit for over 60 y.o. population and comparison with BNT162b2 data fit for the general population).

## List of the files contained

### .mat FILES
- **BNT162b2_fit.mat**: Contains the parametrization obtained with BNT162b2 antibodies data for the general population.
- **BNT162b2_fit_over60.mat**: Contains the parametrization obtained with BNT162b2 antibodies data for over 60 years old population.
- **GoelData_over60.mat**: Contains the data from Goel et al. for over 60 y.o. population.
- **Liang_data.mat**: Contains the data retrieved from Liang et al.
- **Liang_fit.mat**: Contains the parametrization obtained with Liang data.
- **RoltgenData_over60.mat**: Contains the data from Roltgen et al. for over 60 y.o. population.

### SCRIPTS
- **Figure_6.m**: Simulates the dynamics of the tissue layer at different doses and reproduces Figure 6.

### FUNCTIONS
- **model_equations.m**: Contains all the equations of the tissue layer of the model.
- **parameters.m**: Contains all the parameters of the tissue layer of the model.
- **simulation.m**: The file that integrates the system of ODEs above.

### .txt FILES
- **readme.txt**: This file.

### .png FILES
- **Antibodies_in_blood_-_over_60_years_old.png**: Figure produced by "Figure_6.m" script saved in PNG format.
