function [t_record, y_record, pars]=simulation(Liang_fit, par_fit, Liang_data, mRNA0_g, Ag_MW, t_2nddose, t_end, Tolerance)

% load parameters collection for the simulation
pars = Parameters(Liang_fit, par_fit, Liang_data, mRNA0_g, Ag_MW);

% Initial condition vector
IC_0 = [pars.mRNA0_pmol; ...  %-----------------------------% mRNA @IS
    pars.NP0_IS;0;0; ...  %---------------------------------% NP;NPL;NPAg  @IS
    pars.MN0_IS;0;0; ...  %---------------------------------% MN;MNL;MNAg  @IS
    pars.mIDC0_IS;0;0;0;0;0;0;0; ...  %---------------------% mDCs @IS
    pars.pIDC0_IS;0;0;0;0;0;0;0; ...  %---------------------% pDCs @IS
    0;0; ... %----------------------------------------------% NPL;NPAg  @LN
    0;0; ... %----------------------------------------------% MNL;MNAg  @LN
    0;0;0;0;0;0;0; ... %------------------------------------% mDCs  @LN
    0;0;0;0;0;0;0; ... %------------------------------------% pDCs  @LN
    pars.NT0;pars.AT_N0;pars.MT0;pars.AT_M0;pars.FT0; ...   % T-cells @LN
    pars.NB0;pars.AB_N0;pars.MB0;pars.AB_M0; ... %----------% B-cells  @LN
    pars.SP0;pars.LP0; ... %--------------------------------% Plasma cell  @LN
    pars.NP0_BL;... %---------------------------------------% NP @BL
    pars.MN0_BL; ... %--------------------------------------% MN @BL
    pars.mIDC0_BL;... %-------------------------------------% mDC @BL
    pars.pIDC0_BL;... %-------------------------------------% pDC @BL
    zeros(17,1);...  %--------------------------------------% Ab @BL
    zeros(17,1);...  %--------------------------------------% SP @BL
    zeros(17,1);... %---------------------------------------% LP @BL
    zeros(17,1)]; %-----------------------------------------% GCB @LN
    

% ODEs integrator parameters
options = odeset('RelTol', Tolerance, 'AbsTol', Tolerance,'NonNegative',1:length(IC_0));

% 1st dose
[t_1,y_1]=ode15s(@(t,Y) model_equations(t,Y,pars), [0 t_2nddose], IC_0, options);

% 2nd dose
IC_1 = y_1(end,:);
IC_1(1) = IC_1(1)+pars.mRNA0_pmol;   % 2nd dose
[t_2,y_2]=ode15s(@(t,Y) model_equations(t,Y,pars), [0 t_end-t_2nddose], IC_1, options);

% merging the simulations
t_record = [t_1(1:(end-1)); t_2(1:end)+t_1(end)]; % simulation time vector
y_record = [y_1(1:(end-1),:); y_2(1:end,:)]; % simulation dynamics vector

end