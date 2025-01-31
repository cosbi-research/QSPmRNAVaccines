%% Setting the vaccine-specific parameters (Pfizer BNT162b2 or Moderna mRNA-1273)

%Pfizer BNT162b2 vaccine
dose_g = 30e-6; % g of mRNA in Pfizer vaccine (30 micrograms)
t_2nddose = 21; % day at which the second dose is administered
load BNT162b2_fit.mat % loading the parametrization optimized by means of Pfizer antibodies data
Ag_MW = 1377479.8; % (g/mol) molecular weight of BNT162b2 mRNA sequence 

% %Pfizer BNT162b2 vaccine (over 60 years old population)
% dose_g = 30e-6; % g of mRNA in Pfizer vaccine (30 micrograms)
% t_2nddose = 21; % day at which the second dose is administered
% load BNT162b2_fit_over60.mat % loading the parametrization optimized by means of Pfizer antibodies data
% Ag_MW = 1377479.8; % (g/mol) molecular weight of BNT162b2 mRNA sequence 

% %Moderna mRNA-1273 vaccine
% dose_g = 100e-6; % g of mRNA in Moderna vaccine (100 micrograms)
% t_2nddose = 28; % day at which the second dose is administered
% load mRNA-1273_fit.mat % loading the parametrization optimized by means of Pfizer antibodies data
% Ag_MW = 1377479.8; % (g/mol) molecular weight of BNT162b2 mRNA sequence 


%% Loading antigen prenenting cells data (Liang et al.) and fitted parameters

load Liang_data.mat % Liang data for initial conditions
load Liang_fit.mat % fitted parameters from Liang data

%% general parameters for the simulation

t_end = t_2nddose+38*7; % simulation end-point (38 weeks after the second dose)
Tolerance = 1e-10; % Tolerance for the integration (MATLAB default is 1e-6)

%% simulation of the dynamics
[t_record, y_record, pars] = simulation(Liang_fit, opt_param, Liang_data, dose_g, Ag_MW, t_2nddose, t_end, Tolerance);