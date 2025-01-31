%% Setting the vaccine-specific parameters (Moderna mRNA-1273)

%Moderna mRNA-1273 vaccine
dose_g = 100e-6; % g of mRNA in Moderna vaccine (100 micrograms)
t_2nddose = 28; % day at which the second dose is administered
load mRNA-1273_fit.mat % loading the parametrization optimized by means of Pfizer antibodies data
Ag_MW = 1377479.8; % (g/mol) molecular weight of BNT162b2 mRNA sequence 


%% Loading antigen prenenting cells data (Liang et al.) and fitted parameters

load Liang_data.mat % Liang data for initial conditions
load Liang_fit.mat % fitted parameters from Liang data

%% general parameters for the simulation

t_end = t_2nddose+42*7; % simulation end-point (38 weeks after the second dose)
Tolerance = 1e-10; % Tolerance for the integration (MATLAB default is 1e-6)
color = autumn(4);


%% simulation of the dynamics

[t_record25, y_record25, pars25] = simulation(Liang_fit, opt_param, Liang_data, dose_g/4, Ag_MW, t_2nddose, t_end, Tolerance);
[t_record50, y_record50, pars50] = simulation(Liang_fit, opt_param, Liang_data, dose_g/2, Ag_MW, t_2nddose, t_end, Tolerance);
[t_record, y_record, pars] = simulation(Liang_fit, opt_param, Liang_data, dose_g, Ag_MW, t_2nddose, t_end, Tolerance);


% IgG production

Ab_BL25 = y_record25(:,51+8*pars25.J:50+9*pars25.J);
Ab_tot_BL25=sum(Ab_BL25,2);
IgG25=((Ab_tot_BL25*10^(-12))*pars25.MW_Ab*10^9)/(pars25.V_BL*10^3); %ng/mL

Ab_BL50 = y_record50(:,51+8*pars50.J:50+9*pars50.J);
Ab_tot_BL50 = sum(Ab_BL50, 2);
IgG50 = ((Ab_tot_BL50*10^(-12))*pars50.MW_Ab*10^9)/(pars50.V_BL*10^3); %ng/mL

Ab_BL = y_record(:,51+8*pars.J:50+9*pars.J);
Ab_tot_BL = sum(Ab_BL, 2);
IgG = ((Ab_tot_BL*10^(-12))*pars.MW_Ab*10^9)/(pars.V_BL*10^3); %ng/mL


%% Data loading

load KeshavarzData.mat % data are in µg/mL

KeshavarzMod_times = KeshavarzDataModerna.VarName2(24:end);
KeshavarzMod_data = KeshavarzDataModerna.VarName4(24:end).*1000; % conversion from µg/mL to ng/mL
KeshavarzMod_data = KeshavarzMod_data(KeshavarzMod_data>0);
KeshavarzMod_times = KeshavarzMod_times(KeshavarzMod_data>0);
KeshMod_ord_uniq_times = unique(KeshavarzMod_times);
KeshMod_ord_uniq_data = zeros(length(KeshMod_ord_uniq_times),1);
KeshMod_mob_mean = zeros(length(KeshMod_ord_uniq_times),1);

for t = 1: length(KeshMod_ord_uniq_times)
    eq_times = KeshavarzMod_times == KeshMod_ord_uniq_times(t);
    KeshMod_ord_uniq_data(t) = geomean(KeshavarzMod_data(eq_times));
end

for t = 1: length(KeshMod_mob_mean)
    close_times = abs(KeshMod_ord_uniq_times - KeshMod_ord_uniq_times(t)) <= 6;
    KeshMod_mob_mean(t) = geomean(KeshMod_ord_uniq_data(close_times));
end


load JochumData.mat % data are in U/mL, but are said to be eq to BAU/mL

Phase1times = Phase1Trial.Day-1;
GMC_25ug = Phase1Trial.GMC_25ug;
low95CI_25ug = Phase1Trial.low95CI_25ug;
upp95CI_25ug = Phase1Trial.up95CI_25ug;
GMC_100ug = Phase1Trial.GMC_100ug;
GMR = GMC_25ug./GMC_100ug;
rel_data25 = KeshMod_mob_mean([3,12,19,33]-1).*GMR(2:end);

load KirsteData.mat

GMC_50ug = Phase2Trial.GMC_50ug;
low95CI_50ug = Phase2Trial.low95CI_50ug;
upp95CI_50ug = Phase2Trial.up95CI_50ug;
GMR2 = GMC_50ug./GMC_100ug;
rel_data50 = KeshMod_mob_mean([3,12,19,33]-1).*GMR2(2:end);



%% PLOTS **********************************************************

figname = 'Antibodies in blood';
size = [17, 15];
figure('name',figname,'Units', 'centimeters', 'Position', [0, 0, size]);


subplot(2,2,[1,2])
semilogy(t_record,IgG,'Color', color(1,:),'LineWidth',2);
hold on
semilogy(KeshMod_ord_uniq_times, KeshMod_ord_uniq_data, '.', 'Color', color(1,:),'LineWidth', 1.5, 'MarkerSize', 15)
xlim([0 t_end])
ylim([300 200000])
xlabel('days')
ylabel('RBD antibodies [ng/mL]')
ax = gca; % current axes
ax.FontSize = 12;
legend({'model prediction', 'Keshavarz et al. (100μg)'}, 'FontSize', 12, "Location","best");

subplot(2,2,3)
m = semilogy(t_record25,IgG25,'Color', color(2,:),'LineWidth',2);
hold on
j_1 = errorbar(Phase1times(2:end), rel_data25, rel_data25 - KeshMod_mob_mean([3,12,19,33]-1).*low95CI_25ug(2:end)./GMC_100ug(2:end), KeshMod_mob_mean([3,12,19,33]-1).*upp95CI_25ug(2:end)./GMC_100ug(2:end) - rel_data25, 'v', 'color', color(2,:), 'LineWidth',1.5);
j = semilogy(Phase1times(2:end), rel_data25, 'v', 'color', color(2,:), 'LineWidth',1.5, 'MarkerSize', 8);
xlim([0 t_end])
ylim([300 200000])
xlabel('days')
ylabel('RBD antibodies [ng/mL]')
ax = gca; % current axes
ax.FontSize = 12;
legend([m j], {'model prediction', 'Jochum et al. (25μg)'}, 'FontSize', 12, 'Location','best')   

subplot(2,2,4)
m = semilogy(t_record50,IgG50,'Color', color(3,:),'LineWidth',2);
hold on
k_1 = errorbar(Phase1times(2:end), rel_data50, rel_data50 - KeshMod_mob_mean([3,12,19,33]-1).*low95CI_50ug(2:end)./GMC_100ug(2:end), KeshMod_mob_mean([3,12,19,33]-1).*upp95CI_50ug(2:end)./GMC_100ug(2:end) - rel_data50, '^', 'color', color(3,:), 'LineWidth',1.5);
k = semilogy(Phase1times(2:end), rel_data50, '^', 'color', color(3,:), 'LineWidth',1.5, 'MarkerSize', 8);
xlim([0 t_end])
ylim([300 200000])
xlabel('days')
ylabel('RBD antibodies [ng/mL]')
ax = gca; % current axes
ax.FontSize = 12;
legend([m k], {'model prediction', 'Kirste et al. (50μg)'}, 'FontSize', 12, 'Location', 'best')

saveas(gcf, strcat(strrep(figname, ' ', '_'), '.png'))