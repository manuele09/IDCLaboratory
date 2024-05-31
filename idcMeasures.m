%% Section 3 of Report
clear
mkdir ./Figures/idcMeasures
load ./Data/idcMeasures.mat
% capacitances (F), finger_finger (Ohm), ribbon_ribbon (Ohm), ribbon_finger (Ohm)
% are matrices, where each column represents a sensor, and each row represents a measure;
% length_measures is a table, where each field contains 10 measures relative to each sensors.

% Mean and Standard Deviation
% Table 5
mean_rr = mean(ribbon_ribbon);
std_rr = std(ribbon_ribbon);
% Table 6
mean_rf = mean(ribbon_finger);
std_rf = std(ribbon_finger);
% Table 7
mean_ff = mean(finger_finger);
std_ff = std(finger_finger);
% Table 8
mean_c = mean(capacitances);
std_c = std(capacitances);
row_capacitors = reshape(capacitances(:,2:end), 90, 1);
dStrings = ["Device 1", "Device 2", "Device 3", "Device 4", "Device 5", "Device 6","Device 7", "Device 8", "Device 9", "Device 10"];
%% Resistance Points Plots
% Figure 9
fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(finger_finger, ".", "MarkerSize", 40);
% title("Finger Finger");
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
lgd = legend(dStrings, NumColumns=5);
lgd.FontSize = 16; % Set the desired font size for the legend
set(gca,'FontSize',40)
saveas(fig, "Figures/idcMeasures/finger_finger.png");

% Figure 5
fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(ribbon_ribbon, ".", "MarkerSize", 40);
% title("Ribbon Ribbon");
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
ylim([0.8, 1.5])
lgd = legend(dStrings, NumColumns=5);
lgd.FontSize = 16; % Set the desired font size for the legend
set(gca,'FontSize',40)
saveas(fig, "Figures/idcMeasures/ribbon_ribbon.png");

% Figure 7
fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(ribbon_finger, ".", "MarkerSize", 40);
% title("Ribbon Finger");
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
lgd = legend(dStrings, NumColumns=5);
lgd.FontSize = 16; % Set the desired font size for the legend
set(gca,'FontSize',40)
saveas(fig, "Figures/idcMeasures/ribbon_finger.png");
%% Resistance Bar Plots (mean and std)
% Figure 10
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 1, 1);
bar(mean_ff);
title("Mean")
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
set(gca, 'FontSize', 30)
disp("Mean finger finger (min; max): " + min(mean_ff)+"; "+max(mean_ff))
subplot(2, 1, 2);
bar(std_ff);
title("Standard Deviation")
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
set(gca, 'FontSize', 30)
disp("Std finger finger (min; max): " + min(std_ff)+"; "+max(std_ff))
saveas(fig, "Figures/idcMeasures/bar_finger_finger.png");

% Figure 6
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 1, 1);
bar(mean_rr);
title("Mean")
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
set(gca, 'FontSize', 30)
disp("Mean Ribbon Ribbon (min; max): " + min(mean_rr)+"; "+max(mean_rr))
subplot(2, 1, 2);
bar(std_rr);
title("Standard Deviation")
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
set(gca, 'FontSize', 30)
disp("Std Ribbon Ribbon (min; max): " + min(std_rr)+"; "+max(std_rr))
saveas(fig, "Figures/idcMeasures/bar_ribbon_ribbon.png");

% Figure 8
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 1, 1);
bar(mean_rf);
title("Mean")
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
set(gca, 'FontSize', 30)
disp("Mean Ribbon Finger (min; max): " + min(mean_rf)+"; "+max(mean_rf))

subplot(2, 1, 2);
bar(std_rf);
title("Standard Deviation")
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
set(gca, 'FontSize', 30)
disp("Std Ribbon Finger (min; max): " + min(std_rf)+"; "+max(std_rf))
saveas(fig, "Figures/idcMeasures/bar_ribbon_finger.png");

%% Capacitance Points Plot and Bar Plots (mean and std)
% Figure 12
fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(capacitances/1e-12, ".", "MarkerSize", 40);
ylabel("Capacitance [pF]")
ylim([0.6, 1.5])
xlabel("Measure Iterations")
lgd = legend(["1", "2", "3", "4", "5", "6","7", "8", "9", "10"], NumColumns=10);
lgd.FontSize = 27; % Set the desired font size for the legend
set(gca,'FontSize',40)
saveas(fig, "Figures/idcMeasures/capacitances.png");

% Figure 13
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 1, 1);
bar(mean_c/1e-12);
title("Mean")
ylabel("Capacitance [pF]")
xlabel("Measure Iterations")
set(gca, 'FontSize', 35)

subplot(2, 1, 2);
bar(std_c/1e-12);
title("Standard Deviation")
ylabel("Capacitance [pF]")
xlabel("Measure Iterations")
set(gca, 'FontSize', 35)
saveas(fig, "Figures/idcMeasures/bar_capacitance.png");

%% Inkjet printing camera measures (Table 9)
disp("Mean inkjet printing camera measures:")
disp([mean(length_measures.finger_lenght) ...
     mean(length_measures.finger_width) ...
     mean(length_measures.finger_spacing)])

%% Compare Capacitance (only PET layer) respect to the
% length measures done with the printing machine.

% (parameters of table 1 and 2)
l = 5e-3;           %overlapping finger length
n = 6;              %number of IDC finger pairs
b = 0.3e-3;         %finger width;
d = 0.3e-3;         %finger spacings;
h1 = 140e-6;        %layer 1 thickness;
h2 = 1e-5;          %layer 2 thickness;
h3 = 1e-5;          %layer 3 thickness;
eps1 = 3.5;         %dielectric permittivity of PET (substrate);
eps2 = 1;           %dielectric permittivity of Politopamina (sensitive layer);
eps3 = 1;           %dielectric permittivity of layer 3 (MUT);
lambda = 2*(b + d);

CKim = c_idc3k(eps1,eps2,eps3,h1,h2,h3,b,d,l,n);
disp("Capacitance (previous values): " + CKim/1e-12 + " pF")

l_ = mean(length_measures.finger_lenght); 
b_ = mean(length_measures.finger_width);
d_ = mean(length_measures.finger_spacing);

% using the printed parameters (during printing)
CKim_ = c_idc3k(eps1,eps2,eps3,h1,h2,h3,b_,d_,l_,n);
disp("Capacitance (measured lengths): " + CKim_/1e-12 + " pF")
disp("The difference is: " + abs(CKim_-CKim)/1e-12 + " pF")

% Matching model with the real measured capacitances
% Given the capacitance measures obtained with the LCR meter
% we find what h1 and eps1 parameters should be
% to match the real capacitances.
startH = 100e-6; 
stepH = 10e-6;
endH = 200e-6; 

startEps = 1.23;
stepEps = 0.2;
endEps = 7;

epsVector = zeros(floor((endEps-startEps)/stepEps), 1); % Vector containing the permettivitys of the second layer
hVector = zeros(floor((endH-startH)/stepH), 1); % Vector containing the thicknesses of the second layer
capacitanceVector = zeros(floor((endEps-startEps)/stepEps), 1); % Vector containing the obtained Capacitance values

% Figure 14
fig = figure('units','normalized','outerposition',[0 0 1 1]);
logs.coordinates = [];
logs.capacintance = [];
j = 1;
for newH=startH:stepH:endH
    hVector(j) = newH; % Save for plotting
    i = 1;
    for newEps=startEps:stepEps:endEps        
        capacitanceVector(i) = c_idc3k(newEps,eps2,eps3,newH,h2,h3,b,d,l,n);
        epsVector(i) = newEps; % Save for plotting
        logs.coordinates = [logs.coordinates, [newEps, newH]'];
        logs.capacintance = [logs.capacintance, capacitanceVector(i)];
        i = i+1;
    end
    plot(epsVector, capacitanceVector*1e12, ".-", "MarkerSize", 30, "LineWidth", 2, "DisplayName", strcat("h1= ",string(newH*1e6), " um"));
    hold on
    j = j+1;
end

for i=1:1:length(mean_c)
    err = abs(logs.capacintance - mean_c(i));
    [value, index] = min(err);
    hold on
    coords = logs.coordinates(:, index);
    plot(coords(1), mean_c(i)*1e12, ".", "MarkerSize", 50, "Color", "Black", "HandleVisibility","off");
end

legend("show", "FontSize", 20, NumColumns=10)
legend("Location","northwest")
xlabel("Permittivity [F/m]");
ylabel("Capacitance [pF]");
set(gca,'FontSize', 40)
saveas(fig, "Figures/idcMeasures/matching_capacitance.png");

% Matching model with the real measured capacitances
% Given the capacitance measures obtained with the LCR meter
% we find what b and d parameters should be
% to match the real capacitances.
startD = 0.2e-3; %Finger Spacing
stepD = 0.01e-3;
endD = 0.4e-3;

startB = 0.2e-3; %Finger Width
stepB = 0.01e-3;
endB = 0.4e-3;

DVector = zeros(floor((endB-startB)/stepB), 1); % Vector containing the permettivitys of the second layer
DVector = zeros(floor((endD-startD)/stepD), 1); % Vector containing the thicknesses of the second layer
capacitanceVector = zeros(floor((endB-startB)/stepB), 1); % Vector containing the obtained Capacitance values

% Figure 15
fig = figure('units','normalized','outerposition',[0 0 1 1]);
logs.coordinates = [];
logs.capacintance = [];
j = 1;
for newD=startD:stepD:endD
    DVector(j) = newD; % Save for plotting
    i = 1;
    for newB=startB:stepB:endB        
        capacitanceVector(i) = c_idc3k(eps1,eps2,eps3,h1,h2,h3,newB,newD,l,n);
        DVector(i) = newB; % Save for plotting
        logs.coordinates = [logs.coordinates, [newB, newD]'];
        logs.capacintance = [logs.capacintance, capacitanceVector(i)];
        i = i+1;
    end
    plot(DVector*1e6, capacitanceVector*1e12, ".-", "MarkerSize", 30, "LineWidth", 2, "DisplayName", strcat("spacing= ",string(newD*1e6), " um"));
    hold on
    j = j+1;
end


for i=1:1:length(mean_c)
    err = abs(logs.capacintance - mean_c(i));
    [value, index] = min(err);
    hold on
    eps = logs.coordinates(:, index);
    plot(eps(1)*1e6, mean_c(i)*1e12, ".", "MarkerSize", 50, "Color", "Black", "HandleVisibility","off");
end

legend("show", "FontSize", 15, NumColumns=10)
xlabel("Finger Width [um]");
ylabel("Capacitance [pF]");
set(gca,'FontSize', 40)
saveas(fig, "Figures/idcMeasures/matching_capacitance2.png");
