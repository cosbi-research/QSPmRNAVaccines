%% retrieve experimental data

load Liang_data.mat


%% Setting the product-specific parameters

dose_g = 50e-6;% g of mRNA in Liang dose protocol (50 micrograms)
load Liang_fit.mat % loading the parametrization optimized by means of antigen prenenting cells data (Liang et al.)
Ag_MW= 232340; % single strand molecular weight of MCitrin molecule, g/mol


%% general parameters for the simulation

t_end = 100; % simulation end-point (days)
Tolerance = 1e-10; % Tolerance for the integration (MATLAB default is 1e-6)
color = lines(6);


%% simulation of the dynamics

[t_record, y_record] = simulation(Liang_fit, Liang_data, dose_g, Ag_MW, t_end, Tolerance);



%% retrieve experimental data______________________________________________

% times
Texp = Liang_data.Time./24;

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
mDCpisexp = Liang_data.mDCISp;
mDCLNPlnexp = Liang_data.mDCLNLNP;
mDCplnexp = Liang_data.mDCLNp;
mDCtotexp = Liang_data.mDCT;

%pDC, plasmacytoid dendritic cells
pDCLNPisexp = Liang_data.pDCISLNP;
pDCpisexp = Liang_data.pDCISp;
pDCLNPlnexp = Liang_data.pDCLNLNP;
pDCplnexp = Liang_data.pDCLNp;
pDCtotexp = Liang_data.pDCT;


%% State variable outputs

%IS................................................................
% % Neutrophils at IS, NP_IS---------------------------------------
NP_IS = y_record(:,2); % naive Monocytes at IS
NPL_IS = y_record(:,3); % NP+LNP at IS
NPAg_IS = y_record(:,4); % NP with expressed antigen at IS

% % Monocytes at IS, MN_IS-----------------------------------------
MN_IS = y_record(:,5); % naive Monocytes at IS
MNL_IS = y_record(:,6); % MN+LNP at IS
MNAg_IS = y_record(:,7); % MN with expressed antigen at IS

% % Myeloid dendritc cells at IS, mDC_IS---------------------------
mDC_IS = y_record(:,8); % naive mDC at IS
mDCL_IS = y_record(:,9); % mDC+LNP at IS
mDCAg_IS = y_record(:,10); % mDC with antigen exposition at IS

% % Plasmacytoid dendritc cells at IS, pDC_IS----------------------
pDC_IS = y_record(:,11); % naive pDC at IS
pDCL_IS = y_record(:,12); % pDC+LNP at IS
pDCAg_IS = y_record(:,13); % pDC with  antigen exposition at IS


%draining LN.......................................................
% % Neutrophils at LN, NP_LN---------------------------------------
NPL_LN = y_record(:,14); % NP+LNP at draining LN
NPAg_LN = y_record(:,15); % NP with expressed antigen at draining LN

% % Monocytes at LN, MN_LN-----------------------------------------
MNL_LN = y_record(:,16); % MN+LNP at draining LN
MNAg_LN = y_record(:,17); % MN with expressed antigen at draining LN

% % Myeloid dendritc cells at LN, mDC_LN---------------------------
mDCL_LN = y_record(:,18);
mDCAg_LN = y_record(:,19); % mDC with antigen exposition at draining LN

% % Plasmacytoid dendritc cells at LN, pDC_LN----------------------
pDCL_LN = y_record(:,20);
pDCAg_LN = y_record(:,21); % pDC with antigen exposition at draining LN



%% figure plot ************************************************************

figname = 'Liang Data';
size = [25, 18];
figure('name',figname,'Units', 'centimeters', 'Position', [0, 0, size]);

subplot(4,4,1)
plot(t_record,NP_IS+NPL_IS+NPAg_IS,'LineWidth',1.5,'Color',color(4,:))
hold on
plot(Texp,NPtotexp,'^','MarkerSize',8,'Color',color(4,:),'LineWidth',1.5)
xlim([0 10])
set(gca,'FontSize',12)

subplot(4,4,2)
plot(t_record,MN_IS+MNL_IS+MNAg_IS,'LineWidth',1.5,'Color',color(4,:))
hold on
plot(Texp,MNtotexp,'^','MarkerSize',8,'Color',color(4,:),'LineWidth',1.5)
xlim([0 10])
set(gca,'FontSize',12)

subplot(4,4,3)
plot(t_record,mDC_IS+mDCL_IS+mDCAg_IS,'LineWidth',1.5,'Color',color(4,:))
hold on
plot(Texp,mDCtotexp,'^','MarkerSize',8,'Color',color(4,:),'LineWidth',1.5)
xlim([0 10])
set(gca,'FontSize',12)

subplot(4,4,4);
q=plot(t_record,pDC_IS+pDCL_IS+pDCAg_IS,'LineWidth',1.5,'Color',color(4,:));
hold on
w=plot(Texp,pDCtotexp,'^','MarkerSize',8,'Color',color(4,:),'LineWidth',1.5);
xlim([0 10])
set(gca,'FontSize',12)

subplot(4,4,5)
plot(t_record,NPL_IS,'LineWidth',1.5,'Color',color(1,:))
hold on
plot(t_record,NPAg_IS,'LineWidth',1.5,'Color',color(2,:))
plot(Texp,NPLNPisexp,'o','MarkerSize',8,'Color',color(1,:),'LineWidth',1.5)
plot(Texp,NPpisexp,'*','MarkerSize',8,'Color',color(2,:),'LineWidth',1.5)
xlim([0 10]) 
set(gca,'FontSize',12)

subplot(4,4,6)
plot(t_record,MNL_IS,'LineWidth',1.5,'Color',color(1,:))
hold on
plot(t_record,MNAg_IS,'LineWidth',1.5,'Color',color(2,:))
plot(Texp,MNLNPisexp,'o','MarkerSize',8,'Color',color(1,:),'LineWidth',1.5)
plot(Texp,MNpisexp,'*','MarkerSize',8,'Color',color(2,:),'LineWidth',1.5)
xlim([0 10])
set(gca,'FontSize',12)

subplot(4,4,7);
plot(t_record,mDCL_IS,'LineWidth',1.5,'Color',color(1,:))
hold on
plot(t_record,mDCAg_IS,'LineWidth',1.5,'Color',color(2,:))
plot(Texp,mDCLNPisexp,'o','MarkerSize',8,'Color',color(1,:),'LineWidth',1.5)
plot(Texp,mDCpisexp,'*','MarkerSize',8,'Color',color(2,:),'LineWidth',1.5)
xlim([0 10])
set(gca,'FontSize',12)

subplot(4,4,8)
plot(t_record,pDCL_IS,'LineWidth',1.5,'LineStyle','-','Color',color(1,:))
hold on
plot(t_record,pDCAg_IS,'LineWidth',1.5,'Color',color(2,:))
plot(Texp,pDCLNPisexp,'o','MarkerSize',8,'Color',color(1,:),'LineWidth',1.5)
plot(Texp,pDCpisexp,'*','MarkerSize',8,'Color',color(2,:),'LineWidth',1.5)
xlim([0 10])
set(gca,'FontSize',12)

subplot(4,4,9)
a=plot(t_record,NPL_LN,'LineWidth',1.5,'Color',color(1,:));
hold on
s=plot(t_record,NPAg_LN,'LineWidth',1.5,'Color',color(2,:));
d=plot(Texp,NPLNPlnexp,'o','MarkerSize',8,'Color',color(1,:),'LineWidth',1.5);
f=plot(Texp,NPplnexp,'*','MarkerSize',8,'Color',color(2,:),'LineWidth',1.5);
xlim([0 10])
set(gca,'FontSize',12)

subplot(4,4,10);
plot(t_record,MNL_LN,'LineWidth',1.5,'Color',color(1,:))
hold on
plot(t_record,MNAg_LN,'LineWidth',1.5,'Color',color(2,:))
plot(Texp,MNLNPlnexp,'o','MarkerSize',8,'Color',color(1,:),'LineWidth',1.5)
plot(Texp,MNplnexp,'*','MarkerSize',8,'Color',color(2,:),'LineWidth',1.5)
xlim([0 10])
set(gca,'FontSize',12)

subplot(4,4,11)
plot(t_record,mDCL_LN,'LineWidth',1.5,'Color',color(1,:))
hold on
plot(t_record,mDCAg_LN,'LineWidth',1.5,'Color',color(2,:))
plot(Texp,mDCLNPlnexp,'o','MarkerSize',8,'Color',color(1,:),'LineWidth',1.5)
plot(Texp,mDCplnexp,'*','MarkerSize',8,'Color',color(2,:),'LineWidth',1.5)
xlim([0 10])
set(gca,'FontSize',12)

subplot(4,4,12)
plot(t_record,pDCL_LN,'LineWidth',1.5,'Color',color(1,:))
hold on
plot(t_record,pDCAg_LN,'LineWidth',1.5,'Color',color(2,:))
plot(Texp,pDCLNPlnexp,'o','MarkerSize',8,'Color',color(1,:),'LineWidth',1.5)
plot(Texp,pDCplnexp,'*','MarkerSize',8,'Color',color(2,:),'LineWidth',1.5)
xlim([0 10]) 
set(gca,'FontSize',12)

hL = subplot(4,4,14.5);
poshL = get(hL,'position');     % Getting its position

lgd = legend([q,a,s,w,d,f], {'model prediction','model prediction','model prediction','Liang et al. (total cells)',...
    'Liang et al. (LNP-internalizing cells)','Liang et al. (Ag-expressing cells)'}, 'orientation', 'horizontal', 'NumColumns', 3, ...
    FontSize=12);

set(lgd,'position',poshL);      % Adjusting legend's position
axis(hL,'off');                 % Turning its axis off

saveas(gcf, strcat(strrep(figname, ' ', '_'), '.png'))
