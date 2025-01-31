function dydt=model_equations(t,y,pars)


%% State variables for the differential equations__________________________
% mRNA, -----------------------------------------
mRNA = y(1); % mRNA encapsulated in the LNPs, pmol

%% IS: injection site------------------------------------------------------

% Neutrophils..............................................................
% N, Neutrophils, cells
NP_IS = y(2);
% NP_LNP, Neutrophils + LNP, cells
NPL_IS = y(3);
% NP_Ag, Neutrophils + LNP, cells
NPAg_IS = y(4);

% Monocytes................................................................
% MN, monocytes, cells
MN_IS = y(5); 
% M_LNP, monocytes + LNP, cells
MNL_IS = y(6);
% M_Ag, monocytes + antigen , cells
MNAg_IS = y(7);

% Myeloid Dendritic Cells..................................................
% mDC, myeloid DC, naive cells
mDC_IS = y(8);     %mDCIS_naive
% mDCLNP, myeloid DC + LNP, cells
mDCL_IS = y(9);
% mDCAgLon, myeloid dendritic cells + Low antigen (binding phase), cells
mDCAgLon_IS = y(10);
% mDCAgMon, myeloid dendritic cells + Medium antigen (binding phase) , cells
mDCAgMon_IS = y(11);
% mDCAgHon, myeloid dendritic cells + High antigen , cells
mDCAgH_IS = y(12);
% mDCAgMoff, myeloid dendritic cells + Medium antigen (unbinding phase) , cells
mDCAgMoff_IS = y(13);
% mDCAgLoff, myeloid dendritic cells + Low antigen (unbinding phase), cells
mDCAgLoff_IS = y(14);
%mDCoff_IS, mature myeloid dendritic cells without Ag exposition (basic off)
mDCoff_IS=y(15);

% Plasmacytoid Dendritic Cells.............................................
% pDC, plasmacitoid DC, cells
pDC_IS = y(16);     %pDCIS_naive                        
% pDCLNP, plasmacitoid DC + LNP, cells
pDCL_IS = y(17);
% pDCAgLon, plasmacitoid dendritic cells + Low antigen (binding phase), cells
pDCAgLon_IS = y(18);
% pDCAgMon, plasmacitoid dendritic cells + Medium antigen (unbinding phase), cells
pDCAgMon_IS = y(19);
% pDCAgHon, plasmacitoid dendritic cells + High antigen, cells
pDCAgH_IS = y(20);
% pDCAgMoff, plasmacitoid dendritic cells + Medium antigen (unbinding phase), cells
pDCAgMoff_IS = y(21);
% pDCAgLoff, plasmacitoid dendritic cells + Low antigen (unbinding phase), cells
pDCAgLoff_IS = y(22);
% pDCoff, mature plasmacitoid dendritic cells without antigen exposition, cells
pDCoff_IS = y(23);

%% LYMPH NODE--------------------------------------------------------------
% Neutrophils..............................................................
% NLNP, Neutrophils + LNP, cells
NPL_LN = y(24);
% NLAg, Neutrophils + LNP, cells
NPAg_LN = y(25);

%Monocytes.................................................................
% MN_LNP, monocytes + LNP, cells
MNL_LN = y(26);
% MN_Ag, monocytes + LNP, cells
MNAg_LN = y(27);

%Myeloid Dendritic Cells...................................................
% mDCLNP, myeloid DC + LNP, cells
mDCL_LN = y(28); 
% mDCAgLon, myeloid dendritic cells + Low antigen (binding phase) , cells
mDCAgLon_LN = y(29);
% mDCAgMon, myeloid dendritic cells + Medium antigen (binding phase), cells
mDCAgMon_LN = y(30);
% mDCAgH, myeloid dendritic cells + High antigen , cells
mDCAgH_LN = y(31);
% mDCAgMoff, myeloid dendritic cells + Medium antigen (unbinding phase), cells
mDCAgMoff_LN = y(32);
% mDCAgLoff, myeloid dendritic cells + Low antigen (unbinding phase), cells
mDCAgLoff_LN = y(33);
% mDC, mature myeloid DC without antigen exposition, cells
mDC_LN = y(34);  %mDC_LN basic off

%Plasmacytoid Dendritic Cells..............................................
% pDC_LNP, plasmacitoid DC + LNP, cells
pDCL_LN = y(35);
% pDCAgLon, plasmacitoid dendritic cells + Low antigen (binding phase), cells
pDCAgLon_LN = y(36);
% pDCAgMon, plasmacitoid dendritic cells + Medium antigen (binding phase), cells
pDCAgMon_LN = y(37);
% pDCAgHon, plasmacitoid dendritic cells + High antigen , cells
pDCAgH_LN = y(38);
% pDCAgMoff, plasmacitoid dendritic cells + Medium antigen (unbinding phase), cells
pDCAgMoff_LN = y(39);
% pDCAgLoff, plasmacitoid dendritic cells + Low antigen (unbinding phase), cells
pDCAgLoff_LN = y(40);
% pDC, mature plasmacitoid DC cells without antigen exposition, cells
pDC_LN = y(41); %pDC_LN basic off

% T-cells..................................................................
% NT, naive helper T cells,  cells
NT = y(42);
% AT_N, activated helper T cells derived from NT,  cells
AT_N = y(43);
% MT, memory T helper cell, cells
MT = y(44);
% AT_M,	activated helper T cells derived from MT, cells
AT_M = y(45);
% FT, functional helper T cell, cells
FT = y(46);

% B-cells..................................................................
% NB, naive B cells,   cells
NB = y(47:46+pars.J);
% AB_N, activated B cells derived from naive B cells,  cells
AB_N = y(47+pars.J:46+2*pars.J);
% AB_N2, activated B cells at a "mature" state,  cells
GCB = y(51+9*pars.J:50+10*pars.J);
% MB, memory B cells, cells
MB = y(47+2*pars.J:46+3*pars.J);
% AB_M, activated B cells derived from memory B cells, cells
AB_M = y(47+3*pars.J:46+4*pars.J);
% SP, short-lived plasma (antibody secreting) B cells,  cells
SP = y(47+4*pars.J:46+5*pars.J);
% LP, long-lived (antibody secreting) B cells , cells
LP = y(47+5*pars.J:46+6*pars.J);

%% Blood-------------------------------------------------------------------

% Innate immune system.....................................................
% NP, Neutrophils, cells in blood
NP_BL = y(47+6*pars.J);
% MN, monocytes, cells in blood
MN_BL = y(48+6*pars.J);
% Naive myeloid Dendritic Cells in blood
mDC_BL = y(49+6*pars.J);
% Naive plasmacytoid Dendritic Cells in blood
pDC_BL = y(50+6*pars.J);

% B-cells..................................................................
% SP, short-lived plasma (antibody secreting) B cells,  cells in blood
SP_BL = y(51+6*pars.J:50+7*pars.J);
% LP, long-lived (antibody secreting) B cells , cells in blood
LP_BL = y(51+7*pars.J:50+8*pars.J);

% Antibodies...............................................................
%A, total antibody, pmole in blood
Ab_BL = y(51+8*pars.J:50+9*pars.J);



%% Total cells in the three different maturation levels

%LN------------------------------------------------------------------------
%mDC_LN
mDCAgL_LN=mDCAgLon_LN+mDCAgLoff_LN; % mDC sum of Low antigen exposition cells at LN
mDCAgM_LN=mDCAgMon_LN+mDCAgMoff_LN; % mDC sum of Medium antigen exposition cells at LN
%pDC_LN
pDCAgL_LN=pDCAgLon_LN+pDCAgLoff_LN; % pDC sum of Low antigen exposition cells at LN
pDCAgM_LN=pDCAgMon_LN+pDCAgMoff_LN; % pDC sum of Low antigen exposition cells at LN

%% Cell counts
% number of mature mDC
mDC_mat = mDCAgL_LN*pars.p_L + mDCAgM_LN*pars.p_M + mDCAgH_LN*pars.p_H; 
% number of mature pDC
pDC_mat = pDCAgL_LN*pars.p_L + pDCAgM_LN*pars.p_M + pDCAgH_LN*pars.p_H;
% number of mDC
mDC_tot = mDCAgL_LN + mDCAgM_LN + mDCAgH_LN; 
% number of pDC
pDC_tot = pDCAgL_LN + pDCAgM_LN + pDCAgH_LN; 

if (mDC_mat<pars.p_L)
    mDC_mat = 0;
end
if (pDC_mat<pars.p_L)
    pDC_mat = 0;
end

%% Naive T cell activation parameters
%% mdC contribution to NT activation
% medium number of MHCII-antigen complexes presented on the mDCs' PM

if (mDC_mat < pars.p_L)
    Nmedio_mDC = 0;
else
    Nmedio_mDC = ((mDC_mat).*pars.NmaxMHC_mDC)./(mDC_tot); 
end

% probability that a mDC meets a T-cell
sum_mDC_T = (mDC_tot./...
          (mDC_tot+NT+AT_N+AT_M+MT));

% mDC-Tcell activation parameters
D_mDC_N = sum_mDC_T .* Nmedio_mDC ./ (Nmedio_mDC+pars.KNT);

%% pdC contribution to NT activation
% medium number of MHCII-antigen complexes presented on the pDCs' PM

if (pDC_mat < pars.p_L)
    Nmedio_pDC= 0;
else
    Nmedio_pDC = ((pDC_mat).*pars.NmaxMHC_pDC)./(pDC_tot); 
end

% probability that a pDC meets a T-cell

sum_pDC_T = (pDC_tot./...
          (pDC_tot+NT+AT_N+AT_M+MT));

% pDC-Tcell activation parameters
D_pDC_N = sum_pDC_T .* Nmedio_pDC ./ (Nmedio_pDC+pars.KNT);


%% weighted contribution of DC cells to NT-cells activation
% mDC contribution weight to NT-cell activation
if (mDC_tot+pDC_tot == 0)
    weight_mDC = 0;
else
    weight_mDC = (mDC_tot./(mDC_tot+pDC_tot));
end
% pDC contribution weight to NT-cell activation
if (mDC_tot+pDC_tot == 0)
    weight_pDC = 0;
else
    weight_pDC = (pDC_tot./(mDC_tot+pDC_tot));
end

% Weighted activation parameters for Naive T-cell following DC contact
D_N = D_mDC_N .* weight_mDC + D_pDC_N .* weight_pDC;
      
      
% %%  Memory T-cells activation parameters
% mDC's contribution to MT-cell activation
D_mDC_M = sum_mDC_T .* Nmedio_mDC ./ (Nmedio_mDC+pars.KMT);
% mDC's contribution to MT-cell activation
D_pDC_M = sum_pDC_T .* Nmedio_pDC ./ (Nmedio_pDC+pars.KMT);

% Weighted contribution to Memory T-cell activation following DC contact
D_M = D_mDC_M .* weight_mDC + D_pDC_M .* weight_pDC;


%% The proliferation/differentiation function E for helper T cells

% mDC contribution for naive T-cell proliferation/differentiation
E_mDC_N = sum_mDC_T .* (Nmedio_mDC -pars.KNT) ./ (Nmedio_mDC+pars.KNT);

% mDC contribution for memory T-cell proliferation/differentiation
E_mDC_M = sum_mDC_T .* (Nmedio_mDC -pars.KMT) ./ (Nmedio_mDC+pars.KMT);

% check for NaN values
if isnan(E_mDC_N)
    E_mDC_N = 0;
end

if isnan(E_mDC_M)
    E_mDC_M = 0;
end

%% pDC contribution to proliferation/differentiation function E for helper T cells
% pDC contribution for naive T-cell proliferation/differentiation
E_pDC_N = sum_pDC_T .* (Nmedio_pDC -pars.KNT) ./ (Nmedio_pDC+pars.KNT);

% pDC contribution for memory T-cell proliferation/differentiation
E_pDC_M = sum_pDC_T .* (Nmedio_pDC -pars.KMT) ./ (Nmedio_pDC+pars.KMT);

% check for NaN values
if isnan(E_pDC_N)
    E_pDC_N = 0;
end

if isnan(E_pDC_M)
    E_pDC_M = 0;
end

%% weighted contribution of DC cells to MT-cells activation
% DC contribution for naive T-cell proliferation/differentiation
E_N = E_mDC_N .* weight_mDC + E_pDC_N .* weight_pDC;

% DC contribution for memory T-cell proliferation/differentiation
E_M = E_mDC_M .* weight_mDC + E_pDC_M .* weight_pDC ;


if t<0.1
E_N = max(E_N, 0);
E_M = max(E_M, 0);
end

%% B-cell activation
% BCR occupation
BCR=(NB+AB_N+GCB+AB_M+MB)*pars.BRN/pars.NA*1E12; % BCR is the amount of BCR, pmoleAg_i=Ag(i);

% Ag, total antigen production estimate, pmole
Ag = pars.expAg*(pars.maxMHC_mDC .* (...
    + (abs(pars.p_L-pars.p_M)/2+pars.p_L) .* (mDCAgL_LN)...
    + (abs(pars.p_M-pars.p_H)/2+pars.p_M) .* (mDCAgM_LN)...
    + (abs(pars.p_H-1)/2+pars.p_H).*(mDCAgH_LN))...
    + pars.maxMHC_pDC .* (...
    + (abs(pars.p_L-pars.p_M)/2+pars.p_L) .* (pDCAgL_LN)...
    + (abs(pars.p_M-pars.p_H)/2+pars.p_M) .* (pDCAgM_LN)...
    + (abs(pars.p_H-1)/2+pars.p_H).*(pDCAgH_LN))); %pmol

%free antigen computation (availabale for B-cell activation)
Eq1 = @(x) Ag/pars.V_LN-x*(1+ sum(pars.Ka.*BCR/pars.V_LN./(1+pars.Ka*x)));
Options = optimset('TolX',5e-10);
% Options = optimset('TolX',5e-16);
Agfree_con = fzero(Eq1,Ag/pars.V_LN,Options); 
% Agfree_con = fsolve(Eq1,Ag/pars.V_LN,Options); 

%receptor occupancy
ro = pars.Ka*Agfree_con./ (1 + pars.Ka*Agfree_con);
% Calculate occupied BCR number "R"
R = ro.*pars.BRN;
% Calculate the functions based on BCR occupancy for B cell activation, proliferation or differentiation
F = R./(pars.K_R+R);
G = (1-ro).*F;
H = (R-pars.K_R)./(R+pars.K_R);

% Functions based on help from functional T cells, for B cell activation, proliferation or differentiation
% function for naive B cells
P_N = pars.CC_N*sum(FT)./(pars.CC_N*sum(FT)+sum(NB)+sum(AB_N)+sum(GCB)+sum(AB_M)+sum(MB));
% P_N = pars.CC_N*sum(FT)./(pars.CC_N*sum(FT)+sum(NB));
% function for memory B cells
P_M = pars.CC_M*sum(FT)./(pars.CC_M*sum(FT)+sum(NB)+sum(AB_N)+sum(GCB)+sum(AB_M)+sum(MB));
% P_M = pars.CC_M*sum(FT)./(pars.CC_M*sum(FT)+sum(NB));


%% ODEs____________________________________________________________________

%% INJECTION SITE
% % mRNA encapsulated in the LNPs, pmol____________________________________


dydt(1,1) = - pars.GammapDC*mRNA*pDC_IS...
            - pars.GammamDC*mRNA*mDC_IS...
            - pars.GammaMN*mRNA*MN_IS...
            - pars.GammaNP*mRNA*NP_IS...
            - pars.kdmrna*mRNA - pars.kdeg*mRNA - pars.ksat*mRNA*(1/2)*(1+tanh((mRNA - pars.mRNA_max)/pars.k_slope));

% NP, Neutrophils__________________________________________________________
% NP, Neutrophils, cells : birth (omeostasis) + recruitment - uptake - death - migration IS2BL 
dydt(2,1) = pars.AlphaNP_IS + pars.EtaNP*mRNA/(pars.K_mRNA+mRNA)*NP_BL - pars.GammaNP*pars.SF_Gamma*mRNA*NP_IS - pars.BetaNP*NP_IS   - pars.OmegaNP*NP_IS;
% NPL, Neutrophils + LNP, cells
dydt(3,1) =  + pars.GammaNP*pars.SF_Gamma*mRNA*NP_IS - pars.DeltaNP*NPL_IS - pars.BetaNP*NPL_IS  - pars.MuNP*NPL_IS;
% NPAg, Neutrophils + Ag, cells
dydt(4,1) = + pars.DeltaNP*NPL_IS - pars.BetaNP*NPAg_IS - pars.XiNP*NPAg_IS;

% MN, monocytes____________________________________________________________
% MN, monocytes, cells : birth (omeostasis) + recruitment - uptake - death - migration IS2BL
dydt(5,1) = pars.AlphaMN_IS + pars.EtaMN*mRNA/(pars.K_mRNA+mRNA)*MN_BL - pars.GammaMN*pars.SF_Gamma*mRNA*MN_IS  - pars.BetaMN*MN_IS   - pars.OmegaMN*MN_IS;
% MNL, monocytes + LNP, cells
dydt(6,1) = + pars.GammaMN*pars.SF_Gamma*mRNA*MN_IS - pars.DeltaMN*MNL_IS - pars.BetaMN*MNL_IS  - pars.MuMN*MNL_IS;
% MNAg, monocytes + Ag, cells
dydt(7,1) =  + pars.DeltaMN*MNL_IS - pars.BetaMN*MNAg_IS - pars.XiMN*MNAg_IS;

% IS_________________DC:reversible maturation process______________________

% mDC, myeloid dendritic cells_____________________________________________
% mDC_IS basic on : birth (omeostasis) + recruitment - uptake - death - migration IS2BL
dydt(8,1) = pars.AlphamIDC_IS + pars.EtamDC*mRNA/(pars.K_mRNA+mRNA)*mDC_BL - pars.GammamDC*pars.SF_Gamma*mRNA*mDC_IS - pars.BetamIDC*mDC_IS  - pars.OmegamDC*mDC_IS;
% mDCL, myeloid DC + LNP, cells
dydt(9,1) = + pars.GammamDC*pars.SF_Gamma*mRNA*mDC_IS - pars.DeltamDC*mDCL_IS - pars.BetamDC*mDCL_IS  - pars.MumDC*mDCL_IS;
% mDCAgLon, myeloid dendritic cells + Low antigen , cells
dydt(10,1) = + pars.DeltamDC*mDCL_IS - pars.k_trM_mDC*mDCAgLon_IS   - pars.BetamDC*mDCAgLon_IS - pars.XimDC*mDCAgLon_IS;
% mDCAgMon, myeloid dendritic cells + Medium antigen , cells
dydt(11,1) =  + pars.k_trM_mDC*mDCAgLon_IS - pars.k_trH_mDC*mDCAgMon_IS   - pars.BetamDC*mDCAgMon_IS - pars.XimDC*mDCAgMon_IS;
% mDCAgH, myeloid dendritic cells + High antigen , cells
dydt(12,1) =+ pars.k_trH_mDC*mDCAgMon_IS - pars.k_atrH_mDC*mDCAgH_IS- pars.BetamDC*mDCAgH_IS - pars.XimDC*mDCAgH_IS;
% mDCAgMoff, myeloid dendritic cells + Medium antigen , cells
dydt(13,1) = + pars.k_atrH_mDC*mDCAgH_IS - pars.k_atrM_mDC*mDCAgMoff_IS   - pars.BetamDC*mDCAgMoff_IS - pars.XimDC*mDCAgMoff_IS;
% mDCAgLoff, myeloid dendritic cells + Low antigen , cells
dydt(14,1) = + pars.k_atrM_mDC*mDCAgMoff_IS  - pars.k_atrL_mDC*mDCAgLoff_IS   - pars.BetamDC*mDCAgLoff_IS - pars.XimDC*mDCAgLoff_IS;
% mDCoff, mature plasmacitoid dendritic cells without antigen exposition, cells
dydt(15,1) = pars.k_atrL_mDC*mDCAgLoff_IS   - pars.BetamDC*mDCoff_IS - pars.XimDC*mDCoff_IS;

% pDC, plasmacytoid dendritic cells________________________________________
% pDC_IS basic on : birth (omeostasis) + recruitment - uptake - death - migration IS2BL
dydt(16,1) = pars.AlphapIDC_IS + pars.EtapDC*mRNA/(pars.K_mRNA+mRNA)*pDC_BL - pars.GammapDC*pars.SF_Gamma*mRNA*pDC_IS - pars.BetapIDC*pDC_IS  - pars.OmegapDC*pDC_IS;
% pDCL, plasmacytoid DC + LNP, cells
dydt(17,1) = + pars.GammapDC*pars.SF_Gamma*mRNA*pDC_IS - pars.DeltapDC*pDCL_IS - pars.BetapDC*pDCL_IS  - pars.MupDC*pDCL_IS;
% pDCAgLon, plasmacytoid dendritic cells + Low antigen , cells
dydt(18,1) = + pars.DeltapDC*pDCL_IS - pars.k_trM_pDC*pDCAgLon_IS   - pars.BetapDC*pDCAgLon_IS - pars.XipDC*pDCAgLon_IS;
% pDCAgMon, plasmacytoid dendritic cells + Medium antigen , cells
dydt(19,1) =  + pars.k_trM_pDC*pDCAgLon_IS - pars.k_trH_pDC*pDCAgMon_IS   - pars.BetapDC*pDCAgMon_IS - pars.XipDC*pDCAgMon_IS;
% pDCAgH, plasmacytoid dendritic cells + High antigen , cells
dydt(20,1) =+ pars.k_trH_pDC*pDCAgMon_IS - pars.k_atrH_pDC*pDCAgH_IS- pars.BetapDC*pDCAgH_IS - pars.XipDC*pDCAgH_IS;
% pDCAgMoff, plasmacytoid dendritic cells + Medium antigen , cells
dydt(21,1) = + pars.k_atrH_pDC*pDCAgH_IS - pars.k_atrM_pDC*pDCAgMoff_IS   - pars.BetapDC*pDCAgMoff_IS - pars.XipDC*pDCAgMoff_IS;
% pDCAgLoff, plasmacytoid dendritic cells + Low antigen , cells
dydt(22,1) = + pars.k_atrM_pDC*pDCAgMoff_IS  - pars.k_atrL_pDC*pDCAgLoff_IS   - pars.BetapDC*pDCAgLoff_IS - pars.XipDC*pDCAgLoff_IS;
% pDCoff, mature plasmacytoid dendritic cells without antigen exposition, cells
dydt(23,1) = pars.k_atrL_pDC*pDCAgLoff_IS   - pars.BetapDC*pDCoff_IS - pars.XipDC*pDCoff_IS;


%% LYMPH NODE______________________________________________________________

% NP, neutrophils__________________________________________________________
% NPL, Neutrophils + LNP, cells
dydt(24,1) = pars.MuNP*NPL_IS  - pars.DeltaNP*NPL_LN - pars.BetaNP*NPL_LN;
% NPAg, Neutrophils + Ag, cells
dydt(25,1) = pars.XiNP*NPAg_IS + pars.DeltaNP*NPL_LN - pars.BetaNP*NPAg_LN;

% MN, monocytes____________________________________________________________
% MNL, Monocytes + LNP, cells
dydt(26,1) = pars.MuMN*MNL_IS  - pars.DeltaMN*MNL_LN - pars.BetaMN*MNL_LN;
% MNAg, Neutrophils + Ag, cells
dydt(27,1) = pars.XiMN*MNAg_IS + pars.DeltaMN*MNL_LN - pars.BetaMN*MNAg_LN;

% LN_________________DC:reversible maturation process______________________

% mDC, myeloid dendritic cells_____________________________________________
% mDCL, myeloid DC + LNP, cells
dydt(28,1) = pars.MumDC*mDCL_IS  - pars.DeltamDC*mDCL_LN - pars.BetamDC*mDCL_LN;
% mDCAgLon, myeloid dendritic cells + Low antigen , cells
dydt(29,1) = pars.XimDC*mDCAgLon_IS + pars.DeltamDC*mDCL_LN  - pars.k_trM_mDC*mDCAgLon_LN   - pars.BetamDC*mDCAgLon_LN;
% mDCAgMon, myeloid dendritic cells + Medium antigen , cells
dydt(30,1) = pars.XimDC*mDCAgMon_IS  + pars.k_trM_mDC*mDCAgLon_LN - pars.k_trH_mDC*mDCAgMon_LN   - pars.BetamDC*mDCAgMon_LN;
% mDCAgH, myeloid dendritic cells + High antigen , cells
dydt(31,1) = pars.XimDC*mDCAgH_IS + pars.k_trH_mDC*mDCAgMon_LN - pars.k_atrH_mDC*mDCAgH_LN - pars.BetamDC*mDCAgH_LN ;
% mDCAgMoff, myeloid dendritic cells + Medium antigen , cells
dydt(32,1) = pars.XimDC*mDCAgMoff_IS  + pars.k_atrH_mDC*mDCAgH_LN - pars.k_atrM_mDC*mDCAgMoff_LN  - pars.BetamDC*mDCAgMoff_LN;
% mDCAgLoff, myeloid dendritic cells + Low antigen , cells
dydt(33,1) = pars.XimDC*mDCAgLoff_IS  + pars.k_atrM_mDC*mDCAgMoff_LN  - pars.k_atrL_mDC*mDCAgLoff_LN   - pars.BetamDC*mDCAgLoff_LN;
% mDC, plasmacitoid DC (naive), cells      %pDC_LN basic off
dydt(34,1) = pars.XimDC*mDCoff_IS + pars.k_atrL_mDC*mDCAgLoff_LN  - pars.BetamDC*mDC_LN;

% pDC, plasmacytoid dendritic cells________________________________________
% pDCL, plasmacytoid DC + LNP, cells
dydt(35,1) = pars.MupDC*pDCL_IS  - pars.DeltapDC*pDCL_LN - pars.BetapDC*pDCL_LN;
% pDCAgLon, plasmacytoid dendritic cells + Low antigen , cells
dydt(36,1) = pars.XipDC*pDCAgLon_IS + pars.DeltapDC*pDCL_LN  - pars.k_trM_pDC*pDCAgLon_LN   - pars.BetapDC*pDCAgLon_LN;
% pDCAgMon, plasmacytoid dendritic cells + Medium antigen , cells
dydt(37,1) = pars.XipDC*pDCAgMon_IS  + pars.k_trM_pDC*pDCAgLon_LN - pars.k_trH_pDC*pDCAgMon_LN   - pars.BetapDC*pDCAgMon_LN;
% pDCAgH, plasmacytoid dendritic cells + High antigen , cells
dydt(38,1) = pars.XipDC*pDCAgH_IS + pars.k_trH_pDC*pDCAgMon_LN - pars.k_atrH_pDC*pDCAgH_LN - pars.BetapDC*pDCAgH_LN ;
% pDCAgMoff, plasmacytoid dendritic cells + Medium antigen , cells
dydt(39,1) = pars.XipDC*pDCAgMoff_IS  + pars.k_atrH_pDC*pDCAgH_LN - pars.k_atrM_pDC*pDCAgMoff_LN   - pars.BetapDC*pDCAgMoff_LN;
% pDCAgLoff, plasmacytoid dendritic cells + Low antigen , cells
dydt(40,1) = pars.XipDC*pDCAgLoff_IS  + pars.k_atrM_pDC*pDCAgMoff_LN  - pars.k_atrL_pDC*pDCAgLoff_LN   - pars.BetapDC*pDCAgLoff_LN;
%pDC, plasmacytoid DC (naive), cells      %pDC_LN basic off
dydt(41,1) = pars.XipDC*pDCoff_IS + pars.k_atrL_pDC*pDCAgLoff_LN  - pars.BetapDC*pDC_LN;

% % T-cells__________________________________________________________________
% NT, naive helper T cells,  cells
dydt(42,1) = pars.BetaNT*(pars.NT0-NT) - pars.SigmaNT*D_N.*NT;
% AT_N, activated helper T cells derived from NT,  cells
dydt(43,1) = pars.SigmaNT*D_N.*NT + pars.RhoAT*E_N.*AT_N - pars.BetaAT*AT_N;
% dydt(43,1) = pars.SigmaNT*D_N.*NT + pars.RhoAT*E_NM.*AT_N - pars.BetaAT*AT_N;
% MT, memory T helper cell, cells
dydt(44,1) = pars.RhoAT*(1-E_N)*pars.f1.*AT_N + pars.RhoAT*(1-E_M)*pars.f1.*AT_M- pars.SigmaMT*D_M.*MT - pars.BetaMT*MT;
% dydt(44,1) = pars.RhoAT*(1-E_NM)*pars.f1.*AT_N + pars.RhoAT*(1-E_NM)*pars.f1.*AT_M- pars.SigmaMT*D_M.*MT - pars.BetaMT*MT;
% AT_M,	activated helper T cells derived from MT, cells
dydt(45,1) = pars.SigmaMT*D_M.*MT + pars.RhoAT*E_M.*AT_M - pars.BetaAT*AT_M;
% dydt(45,1) = pars.SigmaMT*D_M.*MT + pars.RhoAT*E_NM.*AT_M - pars.BetaAT*AT_M;
% FT, functional helper T cell, cells
dydt(46,1) = pars.RhoAT*(1-E_N)*(1-pars.f1).*AT_N + pars.RhoAT*(1-E_M)*(1-pars.f1).*AT_M -pars.BetaFT*FT;
% dydt(46,1) = pars.RhoAT*(1-E_NM)*(1-pars.f1).*AT_N + pars.RhoAT*(1-E_NM)*(1-pars.f1).*AT_M -pars.BetaFT*FT;

% B-cells__________________________________________________________________
% NB, naive B cells,   cells
dydt(47:46+pars.J,1) = pars.BetaNB*(pars.NB0 - NB) - pars.SigmaNB*F*P_N.*NB;
% AB_N, activated B cells derived from naive B cells,  cells
dydt(47+pars.J:46+2*pars.J,1) = pars.SigmaNB*G*P_N.*NB - pars.delayB*AB_N - pars.BetaAB*AB_N;
% AB_N2, activated B cells at a "mature" state,  cells
dydt(51+9*pars.J:50+10*pars.J,1) = pars.delayB*AB_N + pars.RhoAB_N*H*P_N.*GCB - pars.BetaAB*GCB;
% MB, memory B cells, cells
dydt(47+2*pars.J:46+3*pars.J,1) = pars.RhoAB_N*(1-H)*P_N*pars.g1.*GCB + pars.RhoAB_M*(1-H)*P_M*pars.g1.*AB_M - pars.SigmaMB*F*P_M.*MB - pars.BetaMB*MB;
% AB_M, activated B cells derived from memory B cells, cells
dydt(47+3*pars.J:46+4*pars.J,1) = pars.SigmaMB*G*P_M.*MB + pars.RhoAB_M*H*P_M.*AB_M - pars.BetaAB*AB_M;
% SP, short-lived plasma (antibody secreting) B cells,  cells
dydt(47+4*pars.J:46+5*pars.J,1) = pars.RhoAB_N*(1-H)*P_N*pars.g2.*GCB + pars.RhoAB_M*(1-H)*P_M*pars.g2.*AB_M - (pars.BetaSP + pars.LambdaSP)*SP;
% LP, long-lived (antibody secreting) B cells , cells
dydt(47+5*pars.J:46+6*pars.J,1) = pars.RhoAB_N*(1-H)*P_N*(1-pars.g1-pars.g2).*GCB...
    +pars.RhoAB_M*(1-H)*P_M*(1-pars.g1-pars.g2).*AB_M - (pars.BetaLP + pars.LambdaLP)*LP;

%% Blood------------------------------------------------------------------

% Innate immune system.....................................................
% NP, Neutrophils, cells in blood
dydt(47+6*pars.J) =  pars.AlphaNP_BL- pars.EtaNP*(mRNA/pars.K_mRNA+mRNA)*NP_BL - pars.BetaNP*NP_BL + pars.OmegaNP*NP_IS;
% MN, monocytes, cells in blood
dydt(48+6*pars.J) = pars.AlphaMN_BL - pars.EtaMN*(mRNA/pars.K_mRNA+mRNA)*MN_BL - pars.BetaMN*MN_BL + pars.OmegaMN*MN_IS;
% mDC, naive myeloid Dendritic Cells in blood
dydt(49+6*pars.J) = pars.AlphamIDC_BL - pars.EtamDC*(mRNA/pars.K_mRNA+mRNA)*mDC_BL - pars.BetamIDC*mDC_BL + pars.OmegamDC*mDC_IS;
% pDC, naive plasmacytoid Dendritic Cells in blood
dydt(50+6*pars.J) = pars.AlphapIDC_BL - pars.EtapDC*(mRNA/pars.K_mRNA+mRNA)*pDC_BL - pars.BetapIDC*pDC_BL + pars.OmegapDC*pDC_IS;


% B-cells..................................................................
% SP, short-lived plasma (antibody secreting) B cells,  cells in blood
dydt(51+6*pars.J:50+7*pars.J) = - pars.BetaSP*SP_BL + pars.LambdaSP*SP;
% LP, long-lived (antibody secreting) B cells , cells in blood
%dydt(51+7*pars.J:50+8*pars.J) = - pars.BetaLP*LP_BL + pars.LambdaLP*LP - pars.EpsilonLP*LP_BL;
dydt(51+7*pars.J:50+8*pars.J) = - pars.BetaLP*LP_BL + pars.LambdaLP*LP;

% Antibodies...............................................................
% A, total antibody, pmole in blood
dydt(51+8*pars.J:50+9*pars.J) = - pars.BetaA.*Ab_BL + pars.AlphaA*(SP_BL/pars.NA*1E12) + pars.AlphaA*(LP_BL/pars.NA*1E12);


end