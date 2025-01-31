function dydt = model_equations(~, X, pars)


%% define state variable___________________________________________________

%IS.......................................................................

mRNA = X(1); % mRNA encapsulated in the LNPs, pmol

% Neutrophils at IS, NP_IS-----------------------------------------------
NP_IS = X(2); % naive Monocytes at IS
NPL_IS = X(3); % NP+LNP at IS
NPAg_IS = X(4); % NP with expressed antigen at IS

% Monocytes at IS, MN_IS-------------------------------------------------
MN_IS = X(5); % naive Monocytes at IS
MNL_IS = X(6); % MN+LNP at IS
MNAg_IS = X(7); % MN with expressed antigen at IS

% Myeloid dendritc cells at IS, mDC_IS-----------------------------------
mDC_IS = X(8); % naive mDC at IS
mDCL_IS = X(9); % mDC+LNP at IS
mDCAg_IS = X(10); % mDC with antigen exposition at IS

% Plasmacytoid dendritc cells at IS, pDC_IS------------------------------
pDC_IS = X(11); % naive pDC at IS
pDCL_IS = X(12); % pDC+LNP at IS
pDCAg_IS = X(13); % pDC with  antigen exposition at IS


%draining LN...............................................................

% Neutrophils at LN, NP_LN-----------------------------------------------
NPL_LN = X(14); % NP+LNP at draining LN
NPAg_LN = X(15); % NP with expressed antigen at draining LN

% Monocytes at LN, MN_LN-------------------------------------------------
MNL_LN = X(16); % MN+LNP at draining LN
MNAg_LN = X(17); % MN with expressed antigen at draining LN

% Myeloid dendritc cells at LN, mDC_LN-----------------------------------
mDCL_LN = X(18);
mDCAg_LN = X(19); % mDC with antigen exposition at draining LN

% Plasmacytoid dendritc cells at LN, pDC_LN------------------------------
pDCL_LN = X(20);
pDCAg_LN = X(21); % pDC with antigen exposition (binding phase) at draining LN


%BL.......................................................................

% Innate immune system.....................................................
NP_BL = X(22); % NP, Neutrophils, cells in blood
MN_BL = X(23); % MN, monocytes, cells in blood
mDC_BL = X(24); %Naive myeloid Dendritic Cells in blood
pDC_BL = X(25); %Naive plasmacytoid Dendritic Cells in blood



%% Differential equations__________________________________________________

% mRNA encapsulated in the LNPs, pmol______________________________________
dydt(1,1) = - pars.GammapDC*mRNA*pDC_IS - pars.GammamDC*mRNA*mDC_IS - pars.GammaMN*mRNA*MN_IS...
            - pars.GammaNP*mRNA*NP_IS - pars.kdmrna*mRNA - pars.kdeg*mRNA;


% NP, neutrophils__________________________________________________________

% NP, Neutrophils, cells
dydt(2,1) = pars.AlphaNP_IS + pars.EtaNP*(mRNA/(pars.K_mRNA+mRNA))*NP_BL - pars.GammaNP*pars.SF_Gamma*mRNA*NP_IS - pars.BetaNP*NP_IS   - pars.OmegaNP*NP_IS;
% NPL, Neutrophils + LNP, cells
dydt(3,1) =  + pars.GammaNP*pars.SF_Gamma*mRNA*NP_IS - pars.DeltaNP*NPL_IS - pars.BetaNP*NPL_IS  - pars.MuNP*NPL_IS;
% NPAg, Neutrophils + Ag, cells
dydt(4,1) = + pars.DeltaNP*NPL_IS - pars.BetaNP*NPAg_IS - pars.XiNP*NPAg_IS;


% MN, monocytes____________________________________________________________

% MN, monocytes, cells
dydt(5,1) = pars.AlphaMN_IS + pars.EtaMN*(mRNA/(pars.K_mRNA+mRNA))*MN_BL - pars.GammaMN*pars.SF_Gamma*mRNA*MN_IS  - pars.BetaMN*MN_IS   - pars.OmegaMN*MN_IS;
% MNL, monocytes + LNP, cells
dydt(6,1) = + pars.GammaMN*pars.SF_Gamma*mRNA*MN_IS - pars.DeltaMN*MNL_IS - pars.BetaMN*MNL_IS  - pars.MuMN*MNL_IS;
% MNAg, monocytes + Ag, cells
dydt(7,1) =  + pars.DeltaMN*MNL_IS - pars.BetaMN*MNAg_IS - pars.XiMN*MNAg_IS;


% mDC, myeloid dendritic cells_____________________________________________

% mDC, plasmacitoid DC, cells   
dydt(8,1) = pars.AlphamIDC_IS + pars.EtamDC*(mRNA/(pars.K_mRNA+mRNA))*mDC_BL - pars.GammamDC*pars.SF_Gamma*mRNA*mDC_IS - pars.BetamIDC*mDC_IS  - pars.OmegamDC*mDC_IS;
% mDCL, myeloid DC + LNP, cells
dydt(9,1) = + pars.GammamDC*pars.SF_Gamma*mRNA*mDC_IS - pars.DeltamDC*mDCL_IS - pars.TaumDC*mDCL_IS  - pars.MumDC*mDCL_IS;
% mDCAg, myeloid dendritic cells + antigen , cells
dydt(10,1) = + pars.DeltamDC*mDCL_IS - pars.TaumDC*mDCAg_IS - pars.XimDC*mDCAg_IS;


% pDC, plasmacytoid dendritic cells________________________________________

% pDC, plasmacytoid DC, cells 
dydt(11,1) = pars.AlphapIDC_IS + pars.EtapDC*(mRNA/(pars.K_mRNA+mRNA))*pDC_BL - pars.GammapDC*pars.SF_Gamma*mRNA*pDC_IS - pars.BetapIDC*pDC_IS  - pars.OmegapDC*pDC_IS;
% pDCL, plasmacytoid DC + LNP, cells
dydt(12,1) = + pars.GammapDC*pars.SF_Gamma*mRNA*pDC_IS - pars.DeltapDC*pDCL_IS - pars.TaupDC*pDCL_IS  - pars.MupDC*pDCL_IS;
% pDCAg, plasmacytoid dendritic cells + antigen , cells
dydt(13,1) = + pars.DeltapDC*pDCL_IS - pars.TaupDC*pDCAg_IS - pars.XipDC*pDCAg_IS;



%% LYMPH NODE______________________________________________________________

% NP, neutrophils__________________________________________________________

% NPL, Neutrophils + LNP, cells
dydt(14,1) = pars.MuNP*NPL_IS  - pars.DeltaNP*NPL_LN - pars.BetaNP*NPL_LN;
% NPAg, Neutrophils + Ag, cells
dydt(15,1) = pars.XiNP*NPAg_IS + pars.DeltaNP*NPL_LN - pars.BetaNP*NPAg_LN;


% MN, monocytes____________________________________________________________

% MNL, Monocytes + LNP, cells
dydt(16,1) = pars.MuMN*MNL_IS  - pars.DeltaMN*MNL_LN - pars.BetaMN*MNL_LN;
% MNAg, Neutrophils + Ag, cells
dydt(17,1) = pars.XiMN*MNAg_IS + pars.DeltaMN*MNL_LN - pars.BetaMN*MNAg_LN;


% mDC, myeloid dendritic cells_____________________________________________

% mDCL, myeloid DC + LNP, cells
dydt(18,1) = pars.MumDC*mDCL_IS  - pars.DeltamDC*mDCL_LN - pars.TaumDC*mDCL_LN;
% mDCAg, myeloid dendritic cells +  antigen , cells
dydt(19,1) = pars.XimDC*mDCAg_IS + pars.DeltamDC*mDCL_LN  - pars.TaumDC*mDCAg_LN;


% pDC, plasmacytoid dendritic cells________________________________________

% pDCL, plasmacytoid DC + LNP, cells
dydt(20,1) = pars.MupDC*pDCL_IS  - pars.DeltapDC*pDCL_LN - pars.TaupDC*pDCL_LN;
% pDCAgLon, plasmacytoid dendritic cells + Low antigen , cells
dydt(21,1) = pars.XipDC*pDCAg_IS + pars.DeltapDC*pDCL_LN  - pars.TaupDC*pDCAg_LN;



%% Blood------------------------------------------------------------------

% Innate immune system.....................................................
% N, Neutrophils, cells in blood
dydt(22,1) =  pars.AlphaNP_BL- pars.EtaNP*(mRNA/(pars.K_mRNA+mRNA))*NP_BL - pars.BetaNP*NP_BL + pars.OmegaNP*NP_IS;
% MN, monocytes, cells in blood
dydt(23,1) = pars.AlphaMN_BL - pars.EtaMN*(mRNA/(pars.K_mRNA+mRNA))*MN_BL - pars.BetaMN*MN_BL + pars.OmegaMN*MN_IS;
%Naive myeloid Dendritic Cells in blood
dydt(24,1) = pars.AlphamIDC_BL- pars.EtamDC*(mRNA/(pars.K_mRNA+mRNA))*mDC_BL - pars.BetamIDC*mDC_BL + pars.OmegamDC*mDC_IS;
%Naive plasmacytoid Dendritic Cells in blood
dydt(25,1) = pars.AlphapIDC_BL- pars.EtapDC*(mRNA/(pars.K_mRNA+mRNA))*pDC_BL - pars.BetapIDC*pDC_BL + pars.OmegapDC*pDC_IS;

end