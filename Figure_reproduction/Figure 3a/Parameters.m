function pars = Parameters(Liang_fit, C0, mRNA0_g, Ag_MW)


%% General parameters------------------------------------------------------

% Conversion hr-1 to day-1
pars.hr2day = 24;
% Lymph node volume
pars.V_LN = 5e-4; %L
% Number of Lymph nodes involved in the immune response
pars.N_LN = 5;
% Blood volume in rhesus macaques, L
pars.V_BL = (mean([68.75, 55.8])*(9.7))/1000;

% mRNA degradation rate: half-life 10h converted into days
pars.kdmrna = log(2)/(10/pars.hr2day);

pars.mRNA0_g  = mRNA0_g;
pars.Ag_MW  = Ag_MW;
% initial amount of mRNA encapsulated in the LNPs, pmol
pars.mRNA0_pmol = (10^12)*pars.mRNA0_g/(pars.Ag_MW); 



%% General parameters
% Muscle uptake rate and degradation rate of LNP degradation rate as muscles consumption
pars.kdeg = Liang_fit(1);
% Recruitment in dependence of the dose
pars.K_mRNA = Liang_fit(2);
% scaling factor for the uptake rate
pars.SF_Gamma = Liang_fit(3);



%% Neutropihls parameters--------------------------------------------------
% NP0_IS: the initial number of Neutrophils cells in injection site
pars.NP0_IS = C0(1); % cells
% OmegaNP: migration rate from IS to BL of NP
pars.OmegaNP = Liang_fit(4);
% BetaNP: death rate of Neutrophils
pars.BetaNP = log(2)/(7/24); %d-1 
% EtaNP: recruitment rate of NP function of LNP
pars.EtaNP = Liang_fit(5);
% GammaNP: LNP uptake of NP
pars.GammaNP = Liang_fit(6);
% DeltaNP: antigen expression rate of NP
pars.DeltaNP = Liang_fit(7);
% MuNP: migration rate to from IS to  LN of NP w/LNP
pars.MuNP = Liang_fit(8);
% XiNP: migration rate from IS to LN of NP w/Ag
pars.XiNP = Liang_fit(9);
% AlphaNP_IS: birth rate of Neutrophils in injection site
pars.AlphaNP_IS = pars.NP0_IS*(pars.BetaNP+pars.OmegaNP);
% NP0_BL: the initial number of Neutrophils cells in blood
pars.NP0_BL = 3.19*10^9*pars.V_BL; % cells
% AlphaNP_BL: birth rate of Neutrophils in blood
pars.AlphaNP_BL = - pars.OmegaNP*pars.NP0_IS + pars.BetaNP*pars.NP0_BL;



%% Monocytes parameters----------------------------------------------------
% NP0_IS: the initial number of Neutrophils cells in injection site
pars.MN0_IS = C0(2); % cells
% OmegaMN: migration rate from IS to BL of MN
pars.OmegaMN = Liang_fit(10);
% BetaMN: death rate of Monocytes
pars.BetaMN = pars.hr2day*log(2)/24; %d-1 
% EtaMN: recruitment rate of MN function of LNP
pars.EtaMN = Liang_fit(11);
% GammaMN: LNP uptake of MN
pars.GammaMN = Liang_fit(12);
% DeltaMN: antigen expression rate of MN
pars.DeltaMN = Liang_fit(13);
% MuMN: migration rate to from IS to  LN of MN w/LNP
pars.MuMN = Liang_fit(14);
% XiMN: migration rate from IS to LN of MN w/Ag
pars.XiMN = Liang_fit(15);
% AlphaNP_IS: birth rate of Neutrophils in injection site
pars.AlphaMN_IS = pars.MN0_IS*(pars.BetaMN+pars.OmegaMN);
% MN0_BL: the initial number of Neutrophils cells in blood
pars.MN0_BL = ((461+291)/2)*10^6*pars.V_BL; % cells
% AlphaNP_BL: birth rate of Neutrophils in blood
pars.AlphaMN_BL = - pars.OmegaMN*pars.MN0_IS + pars.BetaMN*pars.MN0_BL;



%% mDC parameters...........................................................

% mDC0_IS: the initial number of mDC cells in the injection site
pars.mIDC0_IS = C0(3); % cells
% BetamIDC: death rate of immature mDC
pars.BetamIDC = 0.0924; %d-1 
% BetamDC: death rate of mature mDC
pars.BetamDC = 0.2310; %d-1
% OmegamDC: migration rate from IS to BL of mDC 
pars.OmegamDC = Liang_fit(16);
% EtamDC: recruitment rate of mDC function of LNP
pars.EtamDC = Liang_fit(17);
% GammamDC: LNP uptake of mDC
pars.GammamDC = Liang_fit(18);
% DeltamDC: Ag expression rate for mDC
pars.DeltamDC = Liang_fit(19);
% MumDC: migration rate to from IS to  LN of mDC w/LNP
pars.MumDC = Liang_fit(20);
% XimDC: migration rate from IS to LN of mDC w/Ag
pars.XimDC = Liang_fit(21);
% TaumDC: death rate of mature mDC (at IS and LN)
pars.TaumDC = Liang_fit(28);
% AlphapIDC_IS: birth rate of immature pDC in IS
pars.AlphamIDC_IS = pars.mIDC0_IS*(pars.BetamIDC+pars.OmegamDC);
% mDC0_BL: the initial number of pDC cells in blood
pars.mIDC0_BL = ((9.97+14.80)/2)*10^6*pars.V_BL; % cells 
% AlphapIDC_BL: birth rate of immature pDC in BL
pars.AlphamIDC_BL = - pars.OmegamDC*pars.mIDC0_IS + pars.BetamIDC*pars.mIDC0_BL;



%% pDC parameters...........................................................

% pDC0_IS: the initial number of pDC cells in the injection site
pars.pIDC0_IS = C0(4); % cells
% BetapIDC: death rate of immature pDC
pars.BetapIDC = 0.0924; %d-1 
% BetapDC: death rate of mature pDC
pars.BetapDC = 0.2310; %d-1
% OmegapDC: migration rate from IS to BL of pDC 
pars.OmegapDC = Liang_fit(22);
% EtapDC: recruitment rate of pDC function of LNP
pars.EtapDC = Liang_fit(23);
% GammapDC: LNP uptake of pDC
pars.GammapDC = Liang_fit(24);
% DeltapDC: Ag expression rate for pDC
pars.DeltapDC = Liang_fit(25);
% MupDC: migration rate to from IS to  LN of pDC w/LNP
pars.MupDC = Liang_fit(26);
% XipDC: migration rate from IS to LN of pDC w/Ag
pars.XipDC = Liang_fit(27);
% TaupDC: death rate of mature pDC (at IS and LN)
pars.TaupDC = Liang_fit(29);
% AlphapIDC_IS: birth rate of immature pDC in IS
pars.AlphapIDC_IS = pars.pIDC0_IS*(pars.BetapIDC+pars.OmegapDC);
% pDC0_BL: the initial number of pDC cells in blood
pars.pIDC0_BL = 7*10^6*pars.V_BL; % cells 
% AlphapIDC_BL: birth rate of immature pDC in BL
pars.AlphapIDC_BL = - pars.OmegapDC*pars.pIDC0_IS + pars.BetapIDC*pars.pIDC0_BL;



end