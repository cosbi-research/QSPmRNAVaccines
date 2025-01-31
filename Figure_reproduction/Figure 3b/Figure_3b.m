%% Setting the vaccine-specific parameters 

%Pfizer BNT162b2 vaccine
dose_g = 30e-6; % g of mRNA in Pfizer vaccine (30 micrograms)
t_2nddose = 21; % day at which the second dose is administered
load BNT162b2_fit.mat % loading the parametrization optimized by means of Pfizer antibodies data
Ag_MW = 1377479.8; % (g/mol) molecular weight of BNT162b2 mRNA sequence 


%% Loading antigen prenenting cells data (Liang et al.) and fitted parameters

load Liang_data.mat % Liang data for initial conditions
load Liang_fit.mat % fitted parameters from Liang data


%% general parameters for the simulation

t_end = t_2nddose+38*7; % simulation end-point (38 weeks after the second dose)
Tolerance = 1e-10; % Tolerance for the integration (MATLAB default is 1e-6)

%% simulation of the dynamics

[t_record30, y_record30, pars30] = simulation(Liang_fit, opt_param, Liang_data, dose_g, Ag_MW, t_2nddose, t_end, Tolerance);
[t_record1, y_record1, pars1] = simulation(Liang_fit, opt_param, Liang_data, dose_g/30, Ag_MW, t_2nddose, t_end, Tolerance);

%A, total antibody, pmole in blood
Ab30 = y_record30(:,51+8*pars30.J:50+9*pars30.J);
Ab_tot_BL30 = sum(Ab30, 2);
Ab1 = y_record1(:,51+8*pars1.J:50+9*pars1.J);
Ab_tot_BL1 = sum(Ab1, 2);
% IgG production
IgG30 = ((Ab_tot_BL30*10^(-12))*pars30.MW_Ab*10^9)/(pars30.V_BL*10^3); %ng/mL
IgG1 = ((Ab_tot_BL1*10^(-12))*pars1.MW_Ab*10^9)/(pars1.V_BL*10^3); %ng/mL


%% Data loading

load SahinData.mat % data are in U/mL

timepoints_Pfizer = pfizerDataWithCI.day(8:13)-1;
points1ng = pfizerDataWithCI.geomean(2:6);
points30ng = pfizerDataWithCI.geomean(22:27);
upp951ng = pfizerDataWithCI.upper95ci(2:6);
upp9530ng = pfizerDataWithCI.upper95ci(22:27);
low951ng = pfizerDataWithCI.lower95ci(2:6);
low9530ng = pfizerDataWithCI.lower95ci(22:27);

LLOQ = 1.15;
points1ng(1) = LLOQ;
fact_of_conv = 10.708442657140134;


load KeshavarzData.mat % data are in µg/mL

Keshavarz_times = KeshavarzData.VarName2(24:end);
Keshavarz_data = KeshavarzData.VarName3(24:end).*1000; % conversion from µg/mL to ng/mL
Kesh_ord_uniq_times = unique(Keshavarz_times);
Kesh_ord_uniq_data = zeros(length(Kesh_ord_uniq_times),1);

for t = 1: length(Kesh_ord_uniq_times)
    eq_times = Keshavarz_times == Kesh_ord_uniq_times(t);
    Kesh_ord_uniq_data(t) = geomean(Keshavarz_data(eq_times));
end

load GoelData.mat % data are in µg/mL

Goel_times = Goeldata.Days(35:end);
Goel_data = Goeldata.antiRBDIgG(35:end).*1000;
Goel_ord_uniq_times = unique(Goel_times);
Goel_ord_uniq_data = zeros(length(Goel_ord_uniq_times),1);

for t = 1: length(Goel_ord_uniq_times)
    eq_times = Goel_times == Goel_ord_uniq_times(t);
    Goel_ord_uniq_data(t) = geomean(Goel_data(eq_times));
end

%% Plot
figname = 'Antibodies in blood';
size = [15, 10];
figure('name',figname,'Units', 'centimeters', 'Position', [0, 0, size]);

color = winter(7);

k = semilogy(Kesh_ord_uniq_times, Kesh_ord_uniq_data, '.', 'Color', color(1,:),'LineWidth', 1.5, 'MarkerSize', 10);
hold on
g = semilogy(Goel_ord_uniq_times, Goel_ord_uniq_data, '*', 'Color', color(1,:),'LineWidth', 1.5, 'MarkerSize', 10);
m30 = semilogy(t_record30, IgG30, 'Color', color(1,:), 'LineWidth', 2);
s30_1 = errorbar(timepoints_Pfizer, fact_of_conv.*points30ng, fact_of_conv.*points30ng - fact_of_conv.*low9530ng, fact_of_conv.*upp9530ng - fact_of_conv.*points30ng, 'o', 'Color', color(1,:), 'LineWidth', 1.5, 'MarkerSize', 10);
s30 = semilogy(timepoints_Pfizer, fact_of_conv.*points30ng, 'o', 'Color', color(1,:), 'LineWidth', 1.5, 'MarkerSize', 10);
m1 = semilogy(t_record1, IgG1, 'Color', color(7,:), 'LineWidth', 2);
s1_1 = errorbar(timepoints_Pfizer(1:end-1), fact_of_conv.*points1ng, fact_of_conv.*points1ng - fact_of_conv.*low951ng, fact_of_conv.*upp951ng - fact_of_conv.*points1ng, 'o', 'Color', color(7,:), 'LineWidth', 1.5, 'MarkerSize', 10);
s1 = semilogy(timepoints_Pfizer(1:end-1), fact_of_conv.*points1ng, 'o', 'Color', color(7,:), 'LineWidth', 1.5, 'MarkerSize', 10);
xlim([0 t_end])
ylim([1 246915])
xlabel('days')
ylabel('RBD antibodies [ng/mL]')
ax = gca; % current axes
ax.FontSize = 12;
legend([m1, m30, s1, s30, k, g], {'model prediction (1μg)', 'model prediction (30μg)', ...
    'Sahin et al. (1μg)', 'Sahin et al. (30μg)', 'Keshavarz et al. (30μg)', 'Goel et al. (30μg)'}, 'FontSize', 12, 'Location', 'south');

saveas(gcf, strcat(strrep(figname, ' ', '_'), '.png'))