function dydt = model_equations_subcellular(~,y,pars_int)
    

%% Variables:

mRNA_e = y(1); % mRNA in the endosome, concentration, M (molar)
mRNA_c = y(2); % concentration, M (molari) in the cell (cytosol)
Ag = y(3); % antigen protein concentration, M (molar) in the cell (cytosol)
P = y(4); % antigen peptides concentration, M (molar) in the cell (cytosol)
MHC_PM_un = y(5); % MHC II unbounded on the plasma membrane
MHC_INT_un = y(6); % MHC II unbounded internalized
MHC_INT_b = y(7); % MHC II bounded internalized
MHC_PM_b = y(8); % MHC II bounded on the plasma membrane



%% Equations:

% mRNA in the endosome
dydt(1,1) = - (pars_int.kesc/pars_int.Ve)*mRNA_e;
% mRNA escaped from the endosome in the cytosol 
dydt(2,1) = + (pars_int.kesc/pars_int.Vc)*mRNA_e - pars_int.kdmrna*mRNA_c;
% antigen protein concentration, M (molari) in the cell (cytosol)
dydt(3,1) = (pars_int.kp/(pars_int.Ve*pars_int.NA))*((pars_int.Ag_MW*pars_int.Vc)/pars_int.mRNA_w)*mRNA_c - (pars_int.kdAg+pars_int.kpept)*Ag;
% antigen peptides concentration, M (molari) in the cell (cytosol)
dydt(4,1) = + pars_int.kpept*Ag + (pars_int.kd/((pars_int.Ve)*pars_int.NA))*MHC_INT_b - (pars_int.ka/((pars_int.Ve)*pars_int.NA))*MHC_INT_un*P - pars_int.kdAg*P;
% MHC II unbounded on the plasma membrane
dydt(5,1) = + pars_int.kout*MHC_INT_un - pars_int.kin*MHC_PM_un + pars_int.kd*MHC_PM_b - pars_int.ksyn*MHC_PM_un;
% MHC II unbounded internalized
dydt(6,1) = + pars_int.kin*MHC_PM_un  - pars_int.kout*MHC_INT_un + pars_int.kd*MHC_INT_b - pars_int.ka*MHC_INT_un*P + pars_int.ksyn*(MHC_PM_un + MHC_PM_b);
% MHC II bounded internalized
dydt(7,1) = + pars_int.ka*MHC_INT_un*P - pars_int.kd*MHC_INT_b -pars_int.kout*MHC_INT_b;
% MHC II bounded on the plasma membrane
dydt(8,1) = + pars_int.kout*MHC_INT_b - pars_int.ksyn*MHC_PM_b - pars_int.kd*MHC_PM_b;

end