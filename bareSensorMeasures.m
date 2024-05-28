clear
load ./Data/v.mat
% Each column (of each loaded variable)
% represents a sensor
% Example: plot(finger_finger(:, 1)) to plot just one of them

mean_ff = mean(finger_finger, 1);
std_ff = std(finger_finger);
mean_rr = mean(ribbon_ribbon, 1);
std_rr = std(ribbon_ribbon);
mean_rf = mean(ribbon_finger, 1);
std_rf = std(ribbon_finger);
mean_c = mean(capacitors, 1);
std_c = std(capacitors);
row_capacitors = reshape(capacitors(:,2:end), 90, 1);

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

fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 1, 1);
X = categorical({'Finger Length', 'Finger Width', 'Finger Spacing'});
bar(X, [mean(misure.finger_lenght) ...
     mean(misure.finger_width) ...
     mean(misure.finger_spacing)]);

title("Mean")
ylabel("Capacitance [F]")
xlabel("Measure Iterations")

subplot(2, 1, 2);
bar(X, [std(misure.finger_lenght) ...
     std(misure.finger_width) ...
     std(misure.finger_spacing)]);
title("Standard Deviation")
ylabel("Capacitance [F]")
xlabel("Measure Iterations")
fontsize(35, "points")
% saveas(fig, "Figures/bar_capacitance.png");
%% Compare Capacitance (only PET layer) respect to the measures
% done with the printing machine.
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

mean_c_pico = mean(capacitors*1e-12, 1);

for i=1:1:length(mean_c_pico)
    err = abs(s.c - mean_c_pico(i));
    [value, index] = min(err);
    hold on
    eps = s.cord(:, index);
    plot(eps(1), mean_c_pico(i), ".", "MarkerSize", 50, "Color", "Black", "HandleVisibility","off");
end
saveas(fig, "Figures/matching_capacitance.png");
