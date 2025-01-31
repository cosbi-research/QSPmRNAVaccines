function [t_record, y_record] = simulation_DC(mRNA0_g, Ag_MW, t_end)


% computation of intracellular model parameters

pars_int = parameters_subcellular(mRNA0_g, Ag_MW);


%% simulation of the dynamics of the intracellular model (in minutes)

% Initial condition vector
IC_0 = [pars_int.mRNA0;0;0;0;pars_int.MHC_PM_un_0;pars_int.MHC_INT_un_0;0;0];

options = odeset('RelTol',1e-6, 'AbsTol',1e-6,'NonNegative',1:length(IC_0));

% simulation of the dynamics (in minutes)
[t_record, y_record] = ode15s(@(t,Y) model_equations_subcellular(t,Y,pars_int), [0 60*24*t_end], IC_0, options);


end