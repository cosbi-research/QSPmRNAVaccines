function pars_int = parameters_subcellular(mRNA0_g, Ag_MW)


%% General parameters------------------------------------------------------


% g of mRNA
pars_int.mRNA0_g = mRNA0_g;
% molecular weight of encoded antigen, g/mol
pars_int.Ag_MW = Ag_MW;
% NA: Avogadro constant    
pars_int.NA = 6.022e23;
% initial amount of mRNA encapsulated in the LNPs, pmol
pars_int.mRNA0_pmol = (10^12)*pars_int.mRNA0_g/(pars_int.Ag_MW);
% initial amount of Spike mRNA molecules, #molecules
pars_int.mRNA_init = pars_int.mRNA0_g*pars_int.NA/(pars_int.Ag_MW); 
% weight of a single molecule of mRNA, g
pars_int.mRNA_w = pars_int.mRNA0_g/pars_int.mRNA_init; 
% degradation rate of mRNA (half-life of 10h, converted into minutes)
pars_int.kdmrna = log(2)/(10*60);   


% L, cell volume
pars_int.Vc = 10^(-15); 
% maximum number of mRNa molecules per LNP
pars_int.nmRNAmolecules_max = 5; 
% L, Endosomal volume
pars_int.Ve = 0.04*10^(-15);
%initiation/termination rate for BNT162b2, term/min
pars_int.kp = 0.262*60; 
% min^(-1), Rate constant for internalization of MHC II
pars_int.kin = 0.01; 
% min^(-1), Rate constant for vesicle recycle to surface  of MHC II
pars_int.kout = 0.02; 
% endosomal escape rate (F.Vanoni's thesis), L/hours = L/(60*min)
pars_int.kesc = 5*10^(-4)*(1/60);
% 1/(M*min) - Molari(mol/L), Association constant for MHC II-peptides complexes
pars_int.ka = 100;
% min^(-1), Dissociation constant for MHC II-peptides complexes
pars_int.kd = 4.2*10^(-4); 
% min^(-1), Rate constant for degradation of native antigen
pars_int.kpept = 0.012;
% min^(-1), Rate constant for routing of molecules to lysosome
pars_int.kdAg = 0.01; 
%min^(-1), Rate constant for MHC II synthesis
pars_int.ksyn = 1*10^(-3); 
% #molecules, Number of MHC II in a cell
pars_int.MHC_tot=2*10^5; 
% #molecules of MHC II unbounded on the PM
pars_int.MHC_PM_un_0 = pars_int.MHC_tot;
% #molecules of MHC II unbounded inside the cell
pars_int.MHC_INT_un_0 = 0;
% molari, initial concentration of mRNA in the endosome
pars_int.mRNA0 = pars_int.nmRNAmolecules_max*pars_int.mRNA_w/(pars_int.Ag_MW*pars_int.Ve); 

end