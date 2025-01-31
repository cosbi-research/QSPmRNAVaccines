%% retrieve experimental data

load Liang_data.mat


%% Setting the product-specific parameters

dose_g = 50e-6;% g of mRNA in Liang dose protocol (50 micrograms)
load Liang_fit.mat % loading the parametrization optimized by means of antigen prenenting cells data (Liang et al.)
Ag_MW= 232340; % single strand molecular weight of MCitrin molecule, g/mol


%% general parameters for the simulation

t_end = 100; % simulation end-point (days)
Tolerance = 1e-10; % Tolerance for the integration (MATLAB default is 1e-6)


%% simulation of the dynamics
[t_record, y_record] = simulation(Liang_fit, Liang_data, dose_g, Ag_MW, t_end, Tolerance);