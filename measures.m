clear
load ./Data/v.mat
% Ogni colonna (di ogni variabile caricata) 
% rappresenta un sensore
% Es: plot(finger_finger(:, 1)) per plottarne solo uno

mean_ff = mean(finger_finger, 1);
std_ff = std(finger_finger);
mean_rr = mean(ribbon_ribbon, 1);
std_rr = std(ribbon_ribbon);
mean_rf = mean(ribbon_finger, 1);
std_rf = std(ribbon_finger);
mean_c = mean(capacitors, 1);
std_c = std(capacitors);
row_capacitors = reshape(capacitors(:,2:end), 90, 1);

C1 = 4.13;
C2 = 3.55;
% (120 - x) / 120 = R
disp("The Ratio is: " + C2/C1);
R = C2/C1;
R = R - (1-R)*1;
disp("New ratio is: " + R);
x = 150-150*R;
disp("Maximum potentiometer: " + x);
%%
% (120 + x) / 120 = R
disp("The greatest ratio is: " + max(row_capacitors)/min(row_capacitors));
R = max(row_capacitors)/min(row_capacitors);
R = (R-1)*4+R;
disp("New ratio is: " + R);
x = R*120 - 120;
disp("Maximum potentiometer: " + x);
dStrings = ["Device 1", "Device 2", "Device 3", "Device 4", "Device 5", "Device 6","Device 7", "Device 8", "Device 9", "Device 10"];
%% Resistance Plots
fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(finger_finger, "LineWidth", 6);
% title("Finger Finger");
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
lgd = legend(dStrings, NumColumns=5);
lgd.FontSize = 16; % Set the desired font size for the legend
set(gca,'FontSize',40)
saveas(fig, "Figures/finger_finger.png");

fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(ribbon_ribbon, "LineWidth", 6);
% title("Ribbon Ribbon");
ylabel("Resistance [Ohm]")
xlabel("Measures Iterations")
lgd = legend(dStrings, NumColumns=5);
lgd.FontSize = 16; % Set the desired font size for the legend
set(gca,'FontSize',40)
saveas(fig, "Figures/ribbon_ribbon.png");

fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(ribbon_finger, "LineWidth", 6);
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

%% Capacitance Plots
fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(capacitors, "LineWidth", 6);
ylabel("Capacitance [F]")
xlabel("Measure Iterations")
lgd = legend(["1", "2", "3", "4", "5", "6","7", "8", "9", "10"], NumColumns=10);
lgd.FontSize = 27; % Set the desired font size for the legend
set(gca,'FontSize',40)

saveas(fig, "Figures/capacitances.png");
%% Capacitance Bar Plots (mean and std)

fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 1, 1);
bar(mean_c);
title("Mean")
ylabel("Capacitance [F]")
xlabel("Measure Iterations")

subplot(2, 1, 2);
bar(std_c);
title("Standard Deviation")
ylabel("Capacitance [F]")
xlabel("Measure Iterations")
fontsize(35, "points")
saveas(fig, "Figures/bar_capacitance.png");

%% Inkjet printing camera measures
%Rivedere misura lunghezza
fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(misure.finger_lenght, ".-", "MarkerSize",60, "LineWidth", 6);
disp("Finger Length: " + mean(misure.finger_lenght) + "um")
fontsize(35, "points")
ylabel("Finger Length [um]")
xlabel("Measure Iterations")
saveas(fig, "Figures/finger_length.png");


fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(misure.finger_width, ".-", "MarkerSize",60, "LineWidth", 6);
disp("Finger width: " + mean(misure.finger_width) + "um")
fontsize(35, "points")
ylabel("Finger width [um]")
xlabel("Measure Iterations")
saveas(fig, "Figures/finger_width.png");

fig = figure('units','normalized','outerposition',[0 0 1 1]);
plot(misure.finger_spacing, ".-", "MarkerSize",60, "LineWidth", 6);
disp("Finger Spacing: " + mean(misure.finger_spacing) + "um")
fontsize(35, "points")
ylabel("Finger Spacing [um]")
xlabel("Measure Iterations")
saveas(fig, "Figures/finger_spacing.png");


% labels = ["IDC 1", "IDC 2", "IDC 3", "IDC 4", "IDC 5", "IDC 6", "IDC 7", "IDC 8", "IDC 9", "IDC 10"];

%% Compare Capacitance only Pet 
clear
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

l = 4998.1e-6; % DA RIVEDERE
b = 0.2997e-3;
d = 0.28785e-3;

% using the printed parameters (during printing)
CKim = c_idc3k(eps1,eps2,eps3,h1,h2,h3,b,d,l,n);
disp("Real Values (means) Capacitance: " + CKim/1e-12 + " pF")

%% Matching model with the real measured capacitances
fig = figure('units','normalized','outerposition',[0 0 1 1]);
startH = 0.1e-3; 
stepH = 0.05e-3;
endH = lambda;

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
    plot(epsVector, capacitanceVector, ".-", "MarkerSize", 30, "LineWidth", 2, "DisplayName", strcat("h1= ",string(newH*1e6), " um"));
    hold on
    j = j+1;
end
legend("show", "FontSize", 15, NumColumns=10)
xlabel("Permittivity [F/m]");
ylabel("Capacitance [F]");
set(gca,'FontSize', 35)

load ./Data/v.mat
mean_c = mean(capacitors*1e-12, 1);

for i=1:1:length(mean_c)
    err = abs(s.c - mean_c(i));
    [value, index] = min(err);
    hold on
    eps = s.cord(:, index);
    plot(eps(1), mean_c(i), ".", "MarkerSize", 50, "Color", "Black", "HandleVisibility","off") %"DisplayName", None); %strcat("IDC", string(i)));
    % disp("Real C: "+ mean_c(i)+ ", eps: " + eps(1));
    % disp("Simulato: " + s.c(index))
end
saveas(fig, "Figures/matching_capacitance.png");
