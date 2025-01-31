%% Setting the vaccine-specific parameters 

%Pfizer BNT162b2 vaccine (over 60 years old population)
dose_g = 30e-6; % g of mRNA in Pfizer vaccine (30 micrograms)
t_2nddose = 21; % day at which the second dose is administered
load BNT162b2_fit_over60.mat % loading the parametrization optimized by means of Pfizer over 60 antibodies data
load BNT162b2_fit.mat % loading the parametrization optimized by means of Pfizer antibodies data
Ag_MW = 1377479.8; % (g/mol) molecular weight of BNT162b2 mRNA sequence 


%% Loading antigen prenenting cells data (Liang et al.) and fitted parameters

load Liang_data.mat % Liang data for initial conditions
load Liang_fit.mat % fitted parameters from Liang data


%% general parameters for the simulation

t_end = t_2nddose+40*7; % simulation end-point (40 weeks after the second dose)
Tolerance = 1e-10; % Tolerance for the integration (MATLAB default is 1e-6)
t_plot = 21;
color = winter(7);

%% simulation of the dynamics

[t_record, y_record, pars] = simulation(Liang_fit, opt_param, Liang_data, dose_g, Ag_MW, t_2nddose, t_end, Tolerance);
[t_record_o60, y_record_o60, pars_o60] = simulation(Liang_fit, opt_param_o60, Liang_data, dose_g, Ag_MW, t_2nddose, t_end, Tolerance);

% IgG production
Ab = y_record(:,51+8*pars.J:50+9*pars.J);
Ab_tot_BL = sum(Ab, 2);
IgG = ((Ab_tot_BL*10^(-12))*pars.MW_Ab*10^9)/(pars.V_BL*10^3); %ng/mL

MB = y_record(:,47+2*pars.J:46+3*pars.J);
SP_BL = y_record(:,51+6*pars.J:50+7*pars.J);
LP_BL = y_record(:,51+7*pars.J:50+8*pars.J);

MB6 = y_record_o60(:,47+2*pars.J:46+3*pars.J);
SP_BL6 = y_record_o60(:,51+6*pars.J:50+7*pars.J);
LP_BL6 = y_record_o60(:,51+7*pars.J:50+8*pars.J);


%% Data loading :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

load GoelData_over60.mat % data are in µg/mL

strat_data_age = stratifieddataGoel.Age;
strat_data_days = stratifieddataGoel.Days;
strat_data_IgG = stratifieddataGoel.antiRBDIgG.*1000;

isover60 = strat_data_age >= 60;
timesnot0 = strat_data_days > 0;
isover60andnot0 = and(isover60, timesnot0);
Goel_strat_days_over60 = strat_data_days(isover60andnot0);
Goel_strat_IgG_over60 = strat_data_IgG(isover60andnot0);


% Papazisis data ----------------------------------------------------------

Pap_strat_times = [21+270; 291+30];
Pap_9mo_over60_GMC = 261.98;
Pap_95CI_low = 190.83;
Pap_95CI_upp = 359.66;
fact_of_conv_Pap = 3.393458396112447;


% Roltgen data ----------------------------------------------------------

load RoltgenData_over60.mat

fact_of_conv_Rolt = 1.237560294520265;

logic1=Roltgenstrat1.variable == 'CoV2_RBD_IgG_Mean';
logic2=Roltgenstrat1.age_group == '> 60';
logic3=isnan(Roltgenstrat1.msd_AU); 
logic4=Roltgenstrat1.covid_ever == 'No';
logic = logic1 & logic2 & not(logic3) & logic4;
logic_tot = logic1 & not(logic3) & logic4;

Roltgen_RBD_tot = Roltgenstrat1.msd_AU(logic_tot);
Roltgen_cat_times_tot = Roltgenstrat1.timepoint(logic_tot);
Roltgen_times_tot = zeros(length(Roltgen_cat_times_tot),1);

for i=1:length(Roltgen_cat_times_tot)
    if Roltgen_cat_times_tot(i) == 'D0'
        Roltgen_times_tot(i) = 0;
    elseif Roltgen_cat_times_tot(i) == 'D7'
        Roltgen_times_tot(i) = 7;
    elseif Roltgen_cat_times_tot(i) == 'D21'
        Roltgen_times_tot(i) = 21;
    elseif Roltgen_cat_times_tot(i) == 'D28'
        Roltgen_times_tot(i) = 28;
    elseif Roltgen_cat_times_tot(i) == 'D42'
        Roltgen_times_tot(i) = 42;
    elseif Roltgen_cat_times_tot(i) == 'D90/120'
        Roltgen_times_tot(i) = 100;
    end
end

Roltgen_RBD_over60 = Roltgenstrat1.msd_AU(logic);
Roltgen_cat_times_over60 = Roltgenstrat1.timepoint(logic);
Roltgen_times_over60 = zeros(length(Roltgen_cat_times_over60),1);

for i=1:length(Roltgen_cat_times_over60)
    if Roltgen_cat_times_over60(i) == 'D0'
        Roltgen_times_over60(i) = 0;
    elseif Roltgen_cat_times_over60(i) == 'D7'
        Roltgen_times_over60(i) = 7;
    elseif Roltgen_cat_times_over60(i) == 'D21'
        Roltgen_times_over60(i) = 21;
    elseif Roltgen_cat_times_over60(i) == 'D28'
        Roltgen_times_over60(i) = 28;
    elseif Roltgen_cat_times_over60(i) == 'D42'
        Roltgen_times_over60(i) = 42;
    elseif Roltgen_cat_times_over60(i) == 'D90/120'
        Roltgen_times_over60(i) = 90;
    end
end

Rolt_uniq_times = [7; 21; 28; 42; 90];

Rolt_AU_over60_d7 = Roltgen_RBD_over60(Roltgen_times_over60==7);
Rolt_AU_over60_d21 = Roltgen_RBD_over60(Roltgen_times_over60==21);
Rolt_AU_over60_d28 = Roltgen_RBD_over60(Roltgen_times_over60==28);
Rolt_AU_over60_d42 = Roltgen_RBD_over60(Roltgen_times_over60==42);
Rolt_AU_over60_d90 = Roltgen_RBD_over60(Roltgen_times_over60==90);

geomean_Rolt_over60 = [geomean(Rolt_AU_over60_d7); ...
    geomean(Rolt_AU_over60_d21); geomean(Rolt_AU_over60_d28); ...
    geomean(Rolt_AU_over60_d42); geomean(Rolt_AU_over60_d90)];

length_Rolt_data = [length(Rolt_AU_over60_d7);...
    length(Rolt_AU_over60_d21); length(Rolt_AU_over60_d28);...
    length(Rolt_AU_over60_d42); length(Rolt_AU_over60_d90)];

std_Rolt_over60 = [std(log10(Rolt_AU_over60_d7));...
    std(log10(Rolt_AU_over60_d21)); std(log10(Rolt_AU_over60_d28));...
    std(log10(Rolt_AU_over60_d42)); std(log10(Rolt_AU_over60_d90))];

std_err_Rolt_over60 = std_Rolt_over60./sqrt(length_Rolt_data);

t_value = tinv(0.975, length_Rolt_data - 1);

Rolt_95CI_low = geomean_Rolt_over60./(10.^(t_value.*std_err_Rolt_over60));
Rolt_95CI_upp = geomean_Rolt_over60.*(10.^(t_value.*std_err_Rolt_over60));


%% PLOTS **********************************************************

figname = 'Antibodies in blood - over 60 years old';
size = [27, 25];
figure('name',figname,'Units', 'centimeters', 'Position', [0, 0, size]);

subplot(5,3,[1 6])
m = semilogy(t_record, IgG, 'Color', color(6,:), 'LineWidth', 2);
hold on
k = semilogy(Pap_strat_times(1), Pap_9mo_over60_GMC*fact_of_conv_Pap, 'x', 'color', color(6,:), 'LineWidth',1.5, 'MarkerSize', 7);
k1 = errorbar(Pap_strat_times(1), Pap_9mo_over60_GMC*fact_of_conv_Pap, (Pap_9mo_over60_GMC - Pap_95CI_low)*fact_of_conv_Pap, (- Pap_9mo_over60_GMC + Pap_95CI_upp)*fact_of_conv_Pap, 'x', 'color', color(6,:), 'LineWidth',1.5, 'MarkerSize', 7);
r = semilogy(Rolt_uniq_times, geomean_Rolt_over60.*fact_of_conv_Rolt, 'square', 'color', color(6,:), 'LineWidth',1.5, 'MarkerSize', 7);
r1 = errorbar(Rolt_uniq_times, geomean_Rolt_over60.*fact_of_conv_Rolt, (geomean_Rolt_over60-Rolt_95CI_low).*fact_of_conv_Rolt, (-geomean_Rolt_over60+Rolt_95CI_upp).*fact_of_conv_Rolt, 'square', 'color', color(6,:), 'LineWidth',1.5, 'MarkerSize', 7);
g = semilogy(Goel_strat_days_over60, Goel_strat_IgG_over60, '*', 'color', color(6,:), 'LineWidth',1.5, 'MarkerSize', 7);
xlim([0 t_end])
ylim([10 600000])
xlabel('days')
ylabel('RBD antibodies [ng/mL]')
ax = gca; % current axes
ax.FontSize = 12;
legend([m k r g], {'model prediction', 'Kontopoulou et al.', 'Roltgen et al.', 'Goel et al.'}, 'fontSize', 12, ...
    'Location','best');

subplot(5,3,[7 10])
plot(t_record,sum(MB,2), 'Color', color(1,:), 'LineWidth', 2);
hold on
plot(t_record_o60,sum(MB6,2), 'Color', color(6,:), 'LineWidth', 2);
xlim([0 t_plot])
xlabel('days')
ylabel('Memory B-cells [#cells]')
set(gca,'FontSize',12)

subplot(5,3,[8 11])
plot(t_record,sum(SP_BL,2), 'Color', color(1,:), 'LineWidth', 2);
hold on
plot(t_record_o60,sum(SP_BL6,2), 'Color', color(6,:), 'LineWidth', 2);
xlim([0 t_plot])
xlabel('days')
ylabel('Short living PCs [#cells]')
set(gca,'FontSize',12)

subplot(5,3,[9 12])
p1 = plot(t_record,sum(LP_BL,2), 'Color', color(1,:), 'LineWidth', 2);
hold on
p2 = plot(t_record_o60,sum(LP_BL6,2), 'Color', color(6,:), 'LineWidth', 2);
xlim([0 t_plot])
xlabel('days')
ylabel('Long living PCs [#cells]')
set(gca,'FontSize',12)

hL = subplot(5,3,14);
poshL = get(hL,'position');     % Getting its position

lgd = legend([p1, p2], {'general population', 'over 60 y.o. population'}, "FontSize", 12, 'NumColumns', 2);

set(lgd,'position',poshL);      % Adjusting legend's position
axis(hL,'off');                 % Turning its axis off

set(gca,'FontSize',12);

saveas(gcf, strcat(strrep(figname, ' ', '_'), '.png'))