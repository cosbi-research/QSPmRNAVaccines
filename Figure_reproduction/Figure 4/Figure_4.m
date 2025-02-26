%% Setting the vaccine-specific parameters 

%Pfizer BNT162b2 vaccine
dose_g = 30e-6; % g of mRNA in Pfizer vaccine (30 micrograms)
admin_times1 = [21 291]; % days at which the second and the booster doses are administered
admin_times2 = 70; % days at which the second is administered - extended interval
load BNT162b2_fit.mat % loading the parametrization optimized by means of Pfizer antibodies data
Ag_MW = 1377479.8; % (g/mol) molecular weight of BNT162b2 mRNA sequence 


%% Loading antigen prenenting cells data (Liang et al.) and fitted parameters

load Liang_data.mat % Liang data for initial conditions
load Liang_fit.mat % fitted parameters from Liang data


%% general parameters for the simulation

% simulation end-point
t_end1 = admin_times1(1)+38*7; % simulation end-point (38 weeks after the second dose)
t_end2 = admin_times1(2)+38*7; % simulation end-point (38 weeks after the booster dose)
Tolerance = 1e-10; % Tolerance for the integration (MATLAB default is 1e-6)


%% simulation of the dynamics

[t_record10, y_record10, pars10] = simulation(Liang_fit, opt_param, Liang_data, dose_g/3, Ag_MW, admin_times1(1), t_end1, Tolerance);
[t_record20, y_record20, pars20] = simulation(Liang_fit, opt_param, Liang_data, dose_g*2/3, Ag_MW, admin_times1(1), t_end1, Tolerance);
[t_record30, y_record30, pars30] = simulation_tris(Liang_fit, opt_param, Liang_data, dose_g, Ag_MW, admin_times1, t_end2, Tolerance);
[t_record10w, y_record10w, pars10w] = simulation(Liang_fit, opt_param, Liang_data, dose_g, Ag_MW, admin_times2, t_end1, Tolerance);

% IgG production, BL

Ab10 = y_record10(:,51+8*pars10.J:50+9*pars10.J);
Ab_tot_BL10 = sum(Ab10, 2);
IgG10 = ((Ab_tot_BL10*10^(-12))*pars10.MW_Ab*10^9)/(pars10.V_BL*10^3);

Ab20 = y_record20(:,51+8*pars20.J:50+9*pars20.J);
Ab_tot_BL20 = sum(Ab20, 2);
IgG20 = ((Ab_tot_BL20*10^(-12))*pars20.MW_Ab*10^9)/(pars20.V_BL*10^3);

Ab30 = y_record30(:,51+8*pars30.J:50+9*pars30.J);
Ab_tot_BL30 = sum(Ab30, 2);
IgG30 = ((Ab_tot_BL30*10^(-12))*pars30.MW_Ab*10^9)/(pars30.V_BL*10^3);

Ab10w = y_record10w(:,51+8*pars10w.J:50+9*pars10w.J);
Ab_tot_BL10w = sum(Ab10w, 2);
IgG10w = ((Ab_tot_BL10w*10^(-12))*pars10w.MW_Ab*10^9)/(pars10w.V_BL*10^3);



%% loading data :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

load SahinData.mat

timepoints_Pfizer = pfizerDataWithCI.day(8:13)-1;
points10ng = pfizerDataWithCI.geomean(8:13);
upp9510ng = pfizerDataWithCI.upper95ci(8:13);
low9510ng = pfizerDataWithCI.lower95ci(8:13);
points20ng = pfizerDataWithCI.geomean(15:20);
upp9520ng = pfizerDataWithCI.upper95ci(15:20);
low9520ng = pfizerDataWithCI.lower95ci(15:20);

fact_of_conv = 10.708442657140134; 

load NaaberData.mat

Naaber_times = NaaberData.times(2:end);
Naaber_data = NaaberData.SRBD_Ab(2:end);
Naaber_Q1 = [666; 13985; 8225; 3097; 893];
Naaber_Q3 = [2582; 36616; 17348; 6924; 2463];

load TakeuchiData.mat 

Takeuchi_times = Takeuchidata.day(2:end);
Takeuchi_data = Takeuchidata.IgGSRBD(2:end);
Takeuchi_Q1 = [277; 11372; 5168; 2854; 439; 301];
Takeuchi_Q3 = [1013; 29747; 13143; 8209; 1345; 1038];

% Papazisis data
Papazisis_times = [14; 21+14; 21+90; 21+180; 21+270; 291+30];
Papazisis_geomean = [460.57; 15585.50; 2945.32; 890.06; 445.85; 19871.33];
Papazisis_upp95CI = [665.53; 18196.66; 3515.33; 1069.96; 516.90; 22365.14];
Papazisis_low95CI = [318.73; 13349.04; 2467.74; 740.40; 384.56; 17665.59];

fact_of_conv_AU = 3.393458396112447; 

load PayneData.mat 

Payne_times = [0, 28, 70, 98, 160];

fact_of_conv_P = 0.4783;


%% PLOTS **********************************************************

% data of Sahin's paper (Pfizer clinical trial) are plotted with the 95%CI as error bar

figname = 'Antibodies in blood';
size = [25, 15];
figure('name',figname,'Units', 'centimeters', 'Position', [0, 0, size]);


color = winter(7);
color10w = cool(5);

subplot('Position',[0.1 0.65 0.35 0.30])
m = semilogy(t_record10,IgG10,'Color', color(2,:),'LineWidth',2);
hold on
s10_e = errorbar(timepoints_Pfizer(2:end), fact_of_conv.*points10ng(2:end), fact_of_conv.*points10ng(2:end) - fact_of_conv.*low9510ng(2:end), ...
    fact_of_conv.*upp9510ng(2:end) - fact_of_conv.*points10ng(2:end), 'o', 'Color', color(2,:), 'LineWidth', 1.5, 'MarkerSize', 8);
s10 = semilogy(timepoints_Pfizer(2:end), fact_of_conv.*points10ng(2:end), 'o', 'Color', color(2,:), 'LineWidth', 1.5, 'MarkerSize', 8);
xlim([0 t_end1])
ylim([11 247000])
xlabel('days')
ylabel('RBD antibodies [ng/mL]')
ax = gca; % current axes
ax.FontSize = 14;
legend([m s10], {'model prediction', 'Sahin et al. (10μg)'}, 'FontSize', 14, Location='south')   


subplot('Position',[0.55 0.65 0.35 0.30])
m = semilogy(t_record20,IgG20,'Color', color(3,:),'LineWidth',2);
hold on
s20_e = errorbar(timepoints_Pfizer(2:end), fact_of_conv.*points20ng(2:end), fact_of_conv.*points20ng(2:end) - fact_of_conv.*low9520ng(2:end), fact_of_conv.*upp9520ng(2:end) - fact_of_conv.*points20ng(2:end), 'o', 'Color', color(3,:), 'LineWidth', 1.5, 'MarkerSize', 8);
s20 = semilogy(timepoints_Pfizer(2:end), fact_of_conv.*points20ng(2:end), 'o', 'Color', color(3,:), 'LineWidth', 1.5, 'MarkerSize', 8);
xlim([0 t_end1])
ylim([11 247000])
xlabel('days')
ylabel('RBD antibodies [ng/mL]')
ax = gca; % current axes
ax.FontSize = 14;
legend([m s20], {'model prediction', 'Sahin et al. (20μg)'}, 'FontSize', 14, Location='south')   


subplot('Position',[0.1 0.15 0.45 0.40])
m = semilogy(t_record30,IgG30,'Color', color(4,:),'LineWidth',2);
hold on
k_1 = errorbar(Papazisis_times, fact_of_conv_AU.*Papazisis_geomean, fact_of_conv_AU.*(Papazisis_geomean - Papazisis_low95CI), fact_of_conv_AU.*(Papazisis_upp95CI - Papazisis_geomean), 'x', 'Color', color(4,:), 'LineWidth', 1.5, 'MarkerSize', 8);
k = semilogy(Papazisis_times, fact_of_conv_AU.*Papazisis_geomean, 'x', 'Color', color(4,:), 'LineWidth', 1.5, 'MarkerSize', 8);
n_1 = errorbar(Naaber_times, Naaber_data.*fact_of_conv_AU, fact_of_conv_AU.*(Naaber_data - Naaber_Q1), fact_of_conv_AU.*(Naaber_Q3 - Naaber_data), '*', 'Color', color(4,:),'LineWidth', 1.5, 'MarkerSize', 8);
n = semilogy(Naaber_times, Naaber_data.*fact_of_conv_AU, '*', 'Color', color(4,:),'LineWidth', 1.5, 'MarkerSize', 8);
t_1 = errorbar(Takeuchi_times, Takeuchi_data.*fact_of_conv_AU, fact_of_conv_AU.*(Takeuchi_data - Takeuchi_Q1), fact_of_conv_AU.*(Takeuchi_Q3 - Takeuchi_data), 'diamond', 'Color', color(4,:),'LineWidth', 1.5, 'MarkerSize', 8);
t = semilogy(Takeuchi_times, Takeuchi_data.*fact_of_conv_AU, 'diamond', 'Color', color(4,:),'LineWidth', 1.5, 'MarkerSize', 8);
xlim([0 t_end2])
ylim([11 247000])
xlabel('days')
ylabel('RBD antibodies [ng/mL]')
ax = gca; % current axes
ax.FontSize = 14;
legend([m k n t], {'model prediction', 'Kontopoulou et al. (3 doses)', 'Naaber et al. (2 doses)', 'Takeuchi et al. (2 doses)'},...
    'FontSize', 14, 'Location','southeast')


subplot('Position',[0.65 0.15 0.25 0.40])
m10w = semilogy(t_record10w, IgG10w,'Color', color10w(5,:), 'LineWidth', 1.5);
hold on
p_l_1 = errorbar(Payne_times(2:end), fact_of_conv_P.*Pl_medians, fact_of_conv_P.*Pl_1stquart, ...
    fact_of_conv_P.*Pl_3rdquart, 'pentagram', "Color", color10w(5,:),'LineWidth', 1, 'MarkerSize', 5);
p_l = plot(Payne_times(2:end), fact_of_conv_P.*Pl_medians, 'pentagram', "Color", color10w(5,:),'LineWidth', 1, 'MarkerSize', 5);
xlim([0 t_end1])
ylim([11 247000])
xlabel('days')
ylabel('RBD antibodies [ng/mL]')
ax = gca; % current axes
ax.FontSize = 14;
legend([m10w p_l], {'model prediction ', 'Payne et al. (long interval)'}, 'FontSize', 14, 'Location','southeast')


saveas(gcf, strcat(strrep(figname, ' ', '_'), '.png'))