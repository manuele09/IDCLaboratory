clear
load ./Data/v.mat
% Each column (of each loaded variable)
% represents a sensor
% Example: plot(finger_finger(:, 1)) to plot just one of them
capacitors = capacitors*1e-12; %%%%%%%%%%%%
mean_ff = mean(finger_finger, 1);
std_ff = std(finger_finger);
mean_rr = mean(ribbon_ribbon, 1);
std_rr = std(ribbon_ribbon);
mean_rf = mean(ribbon_finger, 1);
std_rf = std(ribbon_finger);
mean_c = mean(capacitors, 1);
std_c = std(capacitors);
row_capacitors = reshape(capacitors(:,2:end), 90, 1);

%% Resistance Plots
dStrings = ["Device 1", "Device 2", "Device 3", "Device 4", "Device 5", "Device 6","Device 7", "Device 8", "Device 9", "Device 10"];
fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(finger_finger, ".", "MarkerSize", 40);
% title("Finger Finger");
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
lgd = legend(dStrings, NumColumns=5);
lgd.FontSize = 16; % Set the desired font size for the legend
set(gca,'FontSize',40)
saveas(fig, "Figures/finger_finger.png");

fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(ribbon_ribbon, ".", "MarkerSize", 40);
% title("Ribbon Ribbon");
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
ylim([0.8, 1.5])
lgd = legend(dStrings, NumColumns=5);
lgd.FontSize = 16; % Set the desired font size for the legend
set(gca,'FontSize',40)
saveas(fig, "Figures/ribbon_ribbon.png");

fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(ribbon_finger, ".", "MarkerSize", 40);
% title("Ribbon Finger");
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
lgd = legend(dStrings, NumColumns=5);
lgd.FontSize = 16; % Set the desired font size for the legend
set(gca,'FontSize',40)
saveas(fig, "Figures/ribbon_finger.png");
%% Resistance Bar Plots (mean and std)
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 1, 1);
bar(mean_ff);
title("Mean")
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
disp("Mean Min/Max finger finger: " + min(mean_ff)+"/"+max(mean_ff))

subplot(2, 1, 2);
bar(std_ff);
title("Standard Deviation")
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
fontsize(30, "points")
disp("Std Min/Max finger finger: " + min(std_ff)+"/"+max(std_ff))
saveas(fig, "Figures/bar_finger_finger.png");

fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 1, 1);
bar(mean_rr);
title("Mean")
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
disp("Mean Min/Max Ribbon Ribbon: " + min(mean_rr)+"/"+max(mean_rr))

subplot(2, 1, 2);
bar(std_rr);
title("Standard Deviation")
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
fontsize(30, "points")
disp("Std Min/Max Ribbon Ribbon: " + min(std_rr)+"/"+max(std_rr))
saveas(fig, "Figures/bar_ribbon_ribbon.png");


fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 1, 1);
bar(mean_rf);
title("Mean")
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
disp("Mean Min/Max Ribbon Finger: " + min(mean_rf)+"/"+max(mean_rf))

subplot(2, 1, 2);
bar(std_rf);
title("Standard Deviation")
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
fontsize(30, "points")
disp("Std Min/Max Ribbon Finger: " + min(std_rf)+"/"+max(std_rf))
saveas(fig, "Figures/bar_ribbon_finger.png");

%% Capacitance Bar Plots (mean and std)
fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(capacitors, ".", "MarkerSize", 40);
ylabel("Capacitance [F]")
ylim([0.6e-12, 1.5e-12])
xlabel("Measure Iterations")
lgd = legend(["1", "2", "3", "4", "5", "6","7", "8", "9", "10"], NumColumns=10);
lgd.FontSize = 27; % Set the desired font size for the legend
set(gca,'FontSize',40)
saveas(fig, "Figures/capacitances.png");

fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 1, 1);
bar(mean_c/1e-12);
title("Mean")
ylabel("Capacitance [pF]")
xlabel("Measure Iterations")

subplot(2, 1, 2);
bar(std_c/1e-12);
title("Standard Deviation")
ylabel("Capacitance [pF]")
xlabel("Measure Iterations")
fontsize(35, "points")
saveas(fig, "Figures/bar_capacitance.png");

%% Inkjet printing camera measures

% bar(X, [mean(misure.finger_lenght) ...
%      mean(misure.finger_width) ...
%      mean(misure.finger_spacing)]);
% 
% bar(X, [std(misure.finger_lenght) ...
%      std(misure.finger_width) ...
%      std(misure.finger_spacing)]);

%% Compare Capacitance (only PET layer) respect to the measures
% done with the printing machine.
l = 5e-3;           %overlapping finger length
n = 6;              %number of IDC finger pairs
b = 0.3e-3;         %finger width;
d = 0.3e-3;         %finger spacings;

h1 = 140e-6;        %layer 1 thickness;
h2 = 1e-5;          %layer 2 thickness;
h3 = 1e-5;          %layer 3 thickness;
eps1 = 3.5;         %dielectric permittivity of PET (substrate);
eps2 = 1.23;        %dielectric permittivity of Politopamina (sensitive layer);
eps3 = 1;           %dielectric permittivity of layer 3 (MUT);
lambda = 2*(b + d);

% using the setted parameters (during printing)
CKim = c_idc3k(eps1,eps2,eps3,h1,h2,h3,b,d,l,n);
disp("Setted Values Capacitance: " + CKim/1e-12 + " pF")

l = 4998.1e-6; 
b = 0.2997e-3;
d = 0.28785e-3;

% using the printed parameters (during printing)
CKim = c_idc3k(eps1,eps2,eps3,h1,h2,h3,b,d,l,n);
disp("Real Values (means) Capacitance: " + CKim/1e-12 + " pF")

%% Matching model with the real measured capacitances
% Given the capacitance measures obtained with the LCR meter
% we find what values (h1 and eps1) must assume the model
% to match the real capacitances.
fig = figure('units','normalized','outerposition',[0 0 1 1]);

startH = 100e-6; 
stepH = 10e-6;
endH = 200e-6; 

startEps = 1.23;
stepEps = 0.2;
endEps = 6;

epsVector = zeros(floor((endEps-startEps)/stepEps), 1); % Vector containing the permettivitys of the second layer
hVector = zeros(floor((endH-startH)/stepH), 1); % Vector containing the thicknesses of the second layer
capacitanceVector = zeros(floor((endEps-startEps)/stepEps), 1); % Vector containing the obtained Capacitance values

s.cord = [];
s.c = [];
j = 1;
for newH=startH:stepH:endH
    hVector(j) = newH; % Save for plotting
    i = 1;
    for newEps=startEps:stepEps:endEps        
        capacitanceVector(i) = c_idc3k(newEps,eps2,eps3,newH,h2,h3,b,d,l,n);
        epsVector(i) = newEps; % Save for plotting
        s.cord = [s.cord, [newEps, newH]'];
        s.c = [s.c, capacitanceVector(i)];
        i = i+1;
    end
    plot(epsVector, capacitanceVector*1e12, ".-", "MarkerSize", 30, "LineWidth", 2, "DisplayName", strcat("h1= ",string(newH*1e6), " um"));
    hold on
    j = j+1;
end
legend("show", "FontSize", 15, NumColumns=10)
xlabel("Permittivity [F/m]");
ylabel("Capacitance [pF]");
set(gca,'FontSize', 35)

mean_c_pico = mean(capacitors, 1);

for i=1:1:length(mean_c_pico)
    err = abs(s.c - mean_c_pico(i));
    [value, index] = min(err);
    hold on
    eps = s.cord(:, index);
    plot(eps(1), mean_c_pico(i)*1e12, ".", "MarkerSize", 50, "Color", "Black", "HandleVisibility","off");
end
saveas(fig, "Figures/matching_capacitance.png");
%% Matching model with the real measured capacitances
% Given the capacitance measures obtained with the LCR meter
% we find what values (h1 and eps1) must assume the model
% to match the real capacitances.
fig = figure('units','normalized','outerposition',[0 0 1 1]);

startD = 0.2e-3; %d=0.3e-3 Finger Spacing
stepD = 0.01e-3;
endD = 0.4e-3;

startB = 0.2e-3; %b=0.3e-3 Finger Width
stepB = 0.01e-3;
endB = 0.4e-3;

DVector = zeros(floor((endB-startB)/stepB), 1); % Vector containing the permettivitys of the second layer
DVector = zeros(floor((endD-startD)/stepD), 1); % Vector containing the thicknesses of the second layer
capacitanceVector = zeros(floor((endB-startB)/stepB), 1); % Vector containing the obtained Capacitance values

s.cord = [];
s.c = [];
j = 1;
for newD=startD:stepD:endD
    DVector(j) = newD; % Save for plotting
    i = 1;
    for newB=startB:stepB:endB        
        capacitanceVector(i) = c_idc3k(eps1,eps2,eps3,h1,h2,h3,newB,newD,l,n);
        DVector(i) = newB; % Save for plotting
        s.cord = [s.cord, [newB, newD]'];
        s.c = [s.c, capacitanceVector(i)];
        i = i+1;
    end
    plot(DVector*1e6, capacitanceVector*1e12, ".-", "MarkerSize", 30, "LineWidth", 2, "DisplayName", strcat("spacing= ",string(newD*1e6), " um"));
    hold on
    j = j+1;
end
legend("show", "FontSize", 15, NumColumns=10)
xlabel("Finger Width [um]");
ylabel("Capacitance [pF]");
set(gca,'FontSize', 35)

mean_c_pico = mean(capacitors, 1);

for i=1:1:length(mean_c_pico)
    err = abs(s.c - mean_c_pico(i));
    [value, index] = min(err);
    hold on
    eps = s.cord(:, index);
    plot(eps(1)*1e6, mean_c_pico(i)*1e12, ".", "MarkerSize", 50, "Color", "Black", "HandleVisibility","off");
end
saveas(fig, "Figures/matching_capacitance.png");
