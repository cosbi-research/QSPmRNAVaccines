%% Setting the product-specific parameters

dose_g = 30e-6; % g of mRNA in Pfizer vaccine (30 micrograms)
Ag_MW = 1377479.8; % (g/mol) molecular weight of BNT162b2 mRNA sequence 

% general parameters for the simulation
t_end = 5; %[day]


%% simulation of the (molecular) dynamics of DC

[t_record, y_record] = simulation_DC(dose_g, Ag_MW, t_end); %(in minutes)