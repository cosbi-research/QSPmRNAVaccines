function pars = Parameters(Liang_fit, Ab_fit, Liang_data, mRNA0_g, Ag_MW)


%% General parameters---------------------------------------------------

% Conversion hr-1 to day-1
pars.hr2day = 24;
% Lymph node volume
pars.V_LN = 5e-4; %L
% Blood volume
pars.V_BL = 5; %L
% NA: Avogadro constant
pars.NA = 6.022e23;% 
% mRNA half-life 10 h, converted into minutes
pars.kdmrna = log(2)/(10/24); 



%% mRNA and injection site (IS) parameters

% g of mRNA
pars.mRNA0_g = mRNA0_g;
% molecular weight of encoded antigen, g/mol
pars.Ag_MW = Ag_MW;
% initial amount of mRNA encapsulated in the LNPs, pmol
pars.mRNA0_pmol = (10^12)*pars.mRNA0_g/(pars.Ag_MW);

% Muscle uptake rate and degradation rate of LNP by muscle consumption
pars.kdeg = Liang_fit(1);
% Recruitment in dependence of the dose
pars.K_mRNA = Liang_fit(2); 
% scaling factor for the uptake rate Gamma, pmol^-1
pars.SF_Gamma = Liang_fit(3);
% pars.mRNA_max: quantity of mRNA that triggers a regulating behavior at the IS
pars.mRNA_max = Ab_fit(13);
% ksat: mRNA degradation rate due to the regulative behavior
pars.ksat = Ab_fit(14);
% k_slope: slope of the regulating function, pmol
pars.k_slope = 1;



%% Neutropihls (NP) parameters------------------------------------------

% OmegaNP: migration rate from IS to BL of NP
pars.OmegaNP = Liang_fit(4);
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

% NP0_IS: initial number of NP in IS
pars.NP0_IS = Liang_data.NPT(1); % cells
% NP0_BL: initial number of NP cells in BL
pars.NP0_BL = 3.19*10^9*pars.V_BL; % cells
% BetaNP: death rate of NP
pars.BetaNP = pars.hr2day*log(2)/7; %d^(-1)
% AlphaNP_IS: birth rate of NP in IS
pars.AlphaNP_IS = pars.NP0_IS*(pars.BetaNP+pars.OmegaNP);
% AlphaNP_BL: birth rate of NP in BL
pars.AlphaNP_BL = - pars.OmegaNP*pars.NP0_IS + pars.BetaNP*pars.NP0_BL;



%% Monocytes (MN) parameters--------------------------------------------

% OmegaMN: migration rate from IS to BL of MN
pars.OmegaMN = Liang_fit(10);
% EtaMN: recruitment rate of MN function of LNP
pars.EtaMN = Liang_fit(11);
% GammaMN: LNP uptake of MN, day^-1
pars.GammaMN = Liang_fit(12);
% DeltaMN: antigen expression rate of MN
pars.DeltaMN = Liang_fit(13);
% MuMN: migration rate to from IS to LN of MN w/LNP
pars.MuMN = Liang_fit(14);
% XiMN: migration rate from IS to LN of MN w/Ag
pars.XiMN = Liang_fit(15);

% MN0_IS: initial number of MN cells in IS
pars.MN0_IS = Liang_data.MNT(1); % cells
% MN0_BL: initial number of MN cells in BL
pars.MN0_BL = ((461+291)/2)*10^6*pars.V_BL; % cells
% BetaMN: death rate of MN
pars.BetaMN = pars.hr2day*log(2)/24; %d^(-1)
% AlphaMn_IS: birth rate of MN in IS
pars.AlphaMN_IS = pars.MN0_IS*(pars.BetaMN+pars.OmegaMN);
% AlphaMN_BL: birth rate of MN in BL
pars.AlphaMN_BL = - pars.OmegaMN*pars.MN0_IS + pars.BetaMN*pars.MN0_BL;



%% mDC parameters.......................................................

% OmegamDC: migration rate from IS to BL of mDC 
pars.OmegamDC = Liang_fit(16);
% EtamDC: recruitment rate of mDC function of LNP
pars.EtamDC = Liang_fit(17);
% GammamDC: LNP uptake of mDC, day^-1
pars.GammamDC = Liang_fit(18);
% DeltamDC: antigen expression rate of mDC
pars.DeltamDC = Liang_fit(19); 
% MumDC: migration rate to from IS to  LN of mDC w/LNP
pars.MumDC = Liang_fit(20);
% XimDC: migration rate from IS to LN of mDC w/Ag
pars.XimDC = Liang_fit(21);
% BetamDC: death rate of mature mDC
pars.BetamDC = Liang_fit(28); %d-1

% mDC0_IS: initial number of mDC cells in IS
pars.mIDC0_IS = Liang_data.mDCT(1); % cells
% mDC0_BL: the initial number of pDC cells in blood
pars.mIDC0_BL = ((9.97+14.80)/2)*10^6*pars.V_BL; % cells
% BetamIDC: death rate of immature mDC
pars.BetamIDC = 0.0924; %day^(-1) 
% AlphamIDC_IS: birth rate of immature pDC in IS
pars.AlphamIDC_IS = pars.mIDC0_IS*(pars.BetamIDC+pars.OmegamDC);
% AlphapIDC_BL: birth rate of immature pDC in BL
pars.AlphamIDC_BL = - pars.OmegamDC*pars.mIDC0_IS + pars.BetamIDC*pars.mIDC0_BL;



%% pDC parameters.......................................................

% OmegapDC: migration rate from IS to BL of pDC 
pars.OmegapDC = Liang_fit(22);
% EtapDC: recruitment rate of pDC function of LNP
pars.EtapDC = Liang_fit(23);
% GammapDC: LNP uptake of pDC, day^-1
pars.GammapDC = Liang_fit(24);
% DeltapDC: antigen expression rate of pDC
pars.DeltapDC = Liang_fit(25);
% MupDC: migration rate to from IS to  LN of pDC w/LNP
pars.MupDC = Liang_fit(26);
% XipDC: migration rate from IS to LN of pDC w/Ag
pars.XipDC = Liang_fit(27);
% BetapDC: death rate of mature pDC
pars.BetapDC = Liang_fit(29); %d-1

% pDC0_IS: the initial number of pDC cells in IS
pars.pIDC0_IS = Liang_data.pDCT(1); % cells
% pDC0_BL: the initial number of pDC cells in BL
pars.pIDC0_BL = 7*10^6*pars.V_BL; % cells
% BetapIDC: death rate of immature pDC
pars.BetapIDC = 0.0924; %day^(-1) 
% AlphapIDC_IS: birth rate of immature pDC in IS
pars.AlphapIDC_IS = pars.pIDC0_IS*(pars.BetapIDC+pars.OmegapDC);
% AlphapIDC_BL: birth rate of immature pDC in BL
pars.AlphapIDC_BL = - pars.OmegapDC*pars.pIDC0_IS + pars.BetapIDC*pars.pIDC0_BL;



%% Dendritic cells------------------------------------------------------

%DC antigen expression/maturation weights...............................
pars.p_L = 0.2;  % Low antigen express weight
pars.p_M = 0.65;  % Medium antigen express weight
pars.p_H = 1;  % High antigen express weight

% Optimized DC transition rates with molecular layer
pars.k_trM_mDC = 5.984386009847674;
pars.k_trH_mDC = 4.023456291503178;
pars.k_atrH_mDC = 2.001144541718650;
pars.k_atrM_mDC = 1.182525483399790;
pars.k_atrL_mDC = 0.393936426260313;
pars.k_trM_pDC = 5.985073797735384;
pars.k_trH_pDC = 4.025034796217001;
pars.k_atrH_pDC = 2.000405122440113;
pars.k_atrM_pDC = 1.182561965521271;
pars.k_atrL_pDC = 0.393949613591618;
pars.NmaxMHC_mDC = 1.226433347931808e+04;
pars.max_perc_mDC = 0.061321667396590;
pars.NmaxMHC_pDC = 1.226433347931808e+04;
pars.max_perc_pDC = 0.061321667396590;

% maximum of bounded MHC II exposed on the PC, pmol
pars.maxMHC_mDC=pars.NmaxMHC_mDC/pars.NA*1E12;
pars.maxMHC_pDC=pars.NmaxMHC_pDC/pars.NA*1E12;



%% scaling factor for antigen's quantity in LN
pars.expAg = Ab_fit(7);



%% T cells parameters---------------------------------------------------

% NT0: the initial number of naive T cells
pars.NT0 = 1.445e3; % cells
% AT_N0: initial number of activated T cell derived from naive T cells
pars.AT_N0 = 0.0; % cells
% AT_M0: initial number of activated T cell derived from memory T cells
pars.AT_M0 = 0.0; % cells
% MT0: initial number of memory T cell
pars.MT0 = 0.0; % cells
% FT0: initial number of functional T cell
pars.FT0 = 0.0; % cells
% n of T-epitope-MHCII-peptide complexes on DC membrane to achieve 50% activation rate of NT cell
pars.KNT= 400; 
% n of T-epitope-MHCII-peptide complexes on DC membrane to achieve 50% activation rate of MT cell
pars.KMT= 40;

% BetaNT: death rate of naive helper T cells
pars.BetaNT = 0.0029 ; % day-1
%BetaAT: death rate of activated helper T cells
pars.BetaAT = 0.18; % day-1
% BetaMT: death rate of memory helper T cells
pars.BetaMT = 2.7397e-004; % day-1
% BetaFT:	death rate of functional helper T cells
pars.BetaFT = 0.18; % day-1
% f1: differentiation percentage for activated helper T cells to become memory T cells
pars.f1 = 0.5;

% SigmaNT: maximum activation rate from naive T cells
pars.SigmaNT = Ab_fit(8);
% SigmaMT: maximum activation rate from memory T cells
pars.SigmaMT = Ab_fit(9);
% RhoAT: maximum proliferation rate for activated helper T cells
pars.RhoAT = Ab_fit(10);



%% B cells parameters---------------------------------------------------

% NB0: initial number of naive B cells
pars.NB0 = (normpdf(-3:6/16:3)./sum(normpdf(-3:6/16:3))*1.3e9*0.000004)'; % cells
% BRN: number of B cell receptor on each B cell
pars.BRN = 75000;  % number of BCR/cell
% Occupied BRC number to achieve 50% activation rate of naive B cells
pars.K_R = 1; %no unit
% J: the number of B cell subclones
pars.J = 17;
pars.Ka = zeros(pars.J,1);
for j = 1:17
    pars.Ka(j,1) = 1e-6*2^(j-(pars.J+1)/2); % Ka is in the unit of pM-1(L/pmole)
end
% AB_N0: initial number of activated B cells derived from naive B cells
pars.AB_N0 = ones(pars.J,1)*0.0; % cells
% AB_M0: initial number of activated B cells derived from memory B cells
pars.AB_M0 = ones(pars.J,1)*0.0; % cells
% SP0: initial number of short-lived plasma (antibody secreting) B cells
pars.SP0=ones(pars.J,1)*0.0; % cells
%LP0: initial number of long-lived (antibody secreting) B cells
pars.LP0=ones(pars.J,1)*0.0; % cells
% BM0: initial number of memory B cells
pars.MB0=ones(pars.J,1)*0.0; % cells

% SigmaNB: maximum activation rate from naive B cells
pars.SigmaNB = 2.48; 
% BetaAB: death rate of activated B cells
pars.BetaAB = 0.2518; 
% BetaMB: death rate of memory B cells
pars.BetaMB = 7.8278e-005; 
% BetaNB: death rate of Naive B cells
pars.BetaNB = 0.029; 
% g1: percentage for activated B cells to differentiate to memory B cells
pars.g1 = 0.5; % no unit
% g2: percentage for activated B cells to differentiate to short-lived plasma cells
pars.g2 = 0.4; % no unit

% CC_N: the carrying capacity for 1 FT cell to stimulate the activation and proliferation of target naive B cells.
pars.CC_N = Ab_fit(1); % no unit
% CC_M: the carrying capacity for 1 FT cell to stimulate the activation and proliferation of target memory B cells.
pars.CC_M = pars.CC_N*10; % no unit
% SigmaMB: maximum activation rate from memory B cells
pars.SigmaMB = Ab_fit(2); % day-1
% RhoAB_N: maximum proliferation rate for activated B cells derived from NB cells.
pars.RhoAB_N = Ab_fit(3);
% RhoAB_M: maximum proliferation rate for activated B cells derived from MB cells.
pars.RhoAB_M = Ab_fit(4);
% LambdaSP: migration rate LN2BL of SP
pars.LambdaSP = Ab_fit(5); %day-1
% LambdaLP: migration rate LN2BL of LP
pars.LambdaLP = pars.LambdaSP; %day-1 
% BetaSP: death rate of  short-lived plasma cells
pars.BetaSP = Ab_fit(12);
% BetaLP: death rate of  long-lived plasma cells
pars.BetaLP = Ab_fit(11);
% DelayB: Migration rate for the activated B cells to germinal center
pars.delayB = Ab_fit(15);



%% Antibodies-----------------------------------------------------------

pars.Ab0 = zeros(pars.J,1);
% A0: initial total amount of antibody
pars.A0=ones(pars.J,1)*0.0; % pmole

% AlphaA: secretion rate of antibody by plasma cells
pars.AlphaA = 8.64*10^8; % day-1

% BetaA: elimination rate of antibody
pars.BetaA = Ab_fit(6);



end