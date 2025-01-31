function [T,SV] = simulation(Liang_fit, Liang_data, mRNA0_g, Ag_MW,t_end,Tolerance)


%% retrieve experimental data______________________________________________

%NP, neutrophils
NPLNPisexp = Liang_data.NPISLNP; 
NPpisexp = Liang_data.NPISp;  
NPLNPlnexp = Liang_data.NPLNLNP;
NPplnexp = Liang_data.NPLNp;
NPtotexp = Liang_data.NPT;

%MN, monocytes
MNLNPisexp = Liang_data.MNISLNP;
MNpisexp = Liang_data.MNISp;
MNLNPlnexp = Liang_data.MNLNLNP;
MNplnexp = Liang_data.MNLNp;
MNtotexp = Liang_data.MNT;

%mDC, myeloid dendritic cells
mDCLNPisexp = Liang_data.mDCISLNP;
mDCLNPlnexp = Liang_data.mDCLNLNP;
mDCtotexp = Liang_data.mDCT;

%pDC, plasmacytoid dendritic cells
pDCLNPisexp = Liang_data.pDCISLNP;
pDCpisexp = Liang_data.pDCISp;
pDCLNPlnexp = Liang_data.pDCLNLNP;
pDCtotexp = Liang_data.pDCT;


%% initial conditions------------------------------------------------------

%NP, neutrophils
NP0 = NPtotexp(1);
NPLNPis0 = NPLNPisexp(1);
NPpis0 = NPpisexp(1);
NPLNPln0 = NPLNPlnexp(1);
NPpln0 = NPplnexp(1);

%MN, monocytes
MN0 = MNtotexp(1);
MNLNPis0 = MNLNPisexp(1);
MNpis0 = MNpisexp(1);
MNLNPln0 = MNLNPlnexp(1);
MNpln0 = MNplnexp(1);

%mDC, myeloid dendritic cells
mDC0 = mDCtotexp(1);
mDCLNPis0 = mDCLNPisexp(1);
mDCLNPln0 = mDCLNPlnexp(1);

%pDC, plasmacytoid dendritic cells
pDC0 = pDCtotexp(1);
pDCLNPis0 = pDCLNPisexp(1);
pDCpis0 = pDCpisexp(1);
pDCLNPln0 = pDCLNPlnexp(1);


C0 = [NP0;MN0;mDC0;pDC0];



%% Model parameters
pars = Parameters(Liang_fit, C0, mRNA0_g, Ag_MW);



%% Initial condition vector
X0 = [pars.mRNA0_pmol;...                                 % mRNA0 at IS
    NP0;NPLNPis0;NPpis0;...                               % NP0 at IS
    MN0;MNLNPis0;MNpis0;...                               % MN0 at IS
    mDC0;mDCLNPis0;0;...                                  % mDC at IS
    pDC0;pDCLNPis0;pDCpis0;...                            % pDC0 at IS
    NPLNPln0;NPpln0;...                                   % NP0 at LN
    MNLNPln0;MNpln0;...                                   % MN0 at LN
    mDCLNPln0;0;...                                       % mDC at LN
    pDCLNPln0;0;...                                       % pDC0 at LN
    pars.NP0_BL;pars.MN0_BL;pars.mIDC0_BL;pars.pIDC0_BL]; % NP0, MN0, mDC0, pDC0 in BL


%% Integration-------------------------------------------------------------

%options
non_negative = 1:length(X0); %All positive
opt = odeset('RelTol', Tolerance, 'AbsTol', Tolerance,'NonNegative', non_negative);    

[T,SV] = ode15s(@(t,X) model_equations(t, X, pars), [0 t_end], X0, opt);


end