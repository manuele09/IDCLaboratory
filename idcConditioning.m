%% Section 4 of report
clear
mkdir ./Figures/idcConditioning
load ./Data/idcConditioning.mat 

% capacitances contains the 11 capacitance values tested;
% volt_readings is a matrix of 5x11 
% (5 voltage measures for each capacitance)

% Sorts in ascending order the capacitances and the corresponding voltages
[capacitances, ind] = sort(capacitances, "ascend");
volt_readings = volt_readings(:, ind);

%% Mean and Standard Deviation (Figure 23 and Table 10)
volt_mean = mean(volt_readings);
volt_std = std(volt_readings);
volt_no_offset = volt_readings - volt_mean(1); %To remove the offset at 0 Farads
% Expressing the standard deviation in percentage with respect to the mean
volt_std_perc = 100*volt_std/mean(volt_mean); 
rep = mean(volt_std_perc);
disp("The repeatibility is: "+rep + "%")

% Figure 23
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 1, 1)
% Mean bar plot
cat = categorical(capacitances);
bar(cat, volt_mean)
ylabel("Voltage [V]")
xlabel("Capacitance [pF]")
set(gca,'FontSize',30)
title("Mean")
% Std bar plot
subplot(2, 1, 2)
bar(cat, volt_std_perc)
title("Std ")
xlabel("Capacitance [pF]")
ylabel("Voltage in % [V]")
set(gca,'FontSize',30)
saveas(fig, "Figures/idcConditioning/bar_cap_conditioning.png");

%% Getting the model

% Put all the data in two vectors
cap_vec = zeros(11*5, 1);
volt_vec = zeros(11*5, 1);
counter = 1;
for i=1:1:11
    for j=1:1:5
        cap_vec(counter) = capacitances(i);
        volt_vec(counter) = volt_no_offset(j, i);
        counter = counter + 1;
    end
end

% Model Fitting
linearModelThroughOrigin = fittype('a*x', 'independent', 'x', 'coefficients', 'a');
fitResult = fit(cap_vec, volt_vec, linearModelThroughOrigin);
coeffs = coeffvalues(fitResult);
a = coeffs(1);

% Residual computation
residuals = zeros(5, 11);
for i=1:1:11
    residuals(:, i) = abs(volt_no_offset(:, i) - fitResult(capacitances(i)));
end
residuals_vec = reshape(residuals, [55, 1]);
uncertainty = 2*max(std(residuals));

% Trasduction Diagram (Figure 24)
fig = figure('units','normalized','outerposition',[0 0 1 1]);
scatter(cap_vec, volt_vec, 500, rand(length(cap_vec), 3), 'x', 'LineWidth', 2.5, 'HandleVisibility','off');
hold on
plot(cap_vec, fitResult(cap_vec), "LineWidth", 2, "DisplayName", "Model", "Color", "red")
legend("show", "Location","northwest")
ylabel("Voltage (difference) [V]")
xlabel("Capacitance [pF]")
% title("Trasduction Diagram")
set(gca,'FontSize',40)
saveas(fig, "Figures/idcConditioning/trasduction.png");

% Calibration Diagram (Figure 25)
fig = figure('units','normalized','outerposition',[0 0 1 1]);
scatter(volt_vec, cap_vec, 500, rand(length(cap_vec), 3), 'x', 'LineWidth', 2.5,'HandleVisibility','off');
hold on
plot(volt_vec, volt_vec/a, "LineWidth", 2, "DisplayName", "Model", "Color", "red")
hold on
plot(volt_vec, volt_vec/a + uncertainty/a, "-", "LineWidth", 2, "Color", "green", "DisplayName","Uncertainty band")
hold on
plot(volt_vec, volt_vec/a - uncertainty/a,  "-", "LineWidth", 2, "Color", "green", "HandleVisibility", "off")
legend("show", "Location","northwest")
xlabel("Voltage (difference) [V]")
ylabel("Capacitance [pF]")
% title("Calibration Diagram")
set(gca,'FontSize',40)
saveas(fig, "Figures/idcConditioning/calibration.png");