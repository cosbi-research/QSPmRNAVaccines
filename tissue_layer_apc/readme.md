This folder contains MATLAB scripts for simulating the antigen-presenting cell (APC) component of the tissue layer, as described in Step 1 of Supplementary Section S5.

## Contents

### .mat Files
- **Liang_data.mat**: Dataset sourced from Liang et al.
- **Liang_fit.mat**: Parameter set derived from Liang et al.'s data.

### Scripts
- **main.m**: Simulates the dynamics of the APC component within the tissue layer, using the APC parameterization from Liang et al.'s data.

### Functions
- **model_equations.m**: Defines all equations for the APC component of the tissue layer, covering both the injection site and the lymph node.
- **parameters.m**: Contains parameter values specific to the APC component of the tissue layer.
- **simulation.m**: Integrates the system of ODEs for the APC dynamics.

### Text Files
- **README.txt**: This file, providing an overview of folder contents.
