load ./Data/elettronica.mat
index = elettronica(:, 1);
cap_val = elettronica(:, 2);
[val, ind] = sort(cap_val, "ascend");
volt = elettronica(:, 3:end);
volt = volt(ind, :);
cap_val = cap_val(ind, :);
volt = volt';

volt_mean = mean(volt);
volt = volt - volt_mean(1);
volt_std = std(volt);
m = mean(volt_mean);
volt_std_perc = 100*volt_std/m;
rep = mean(volt_std_perc);
disp("The repeatibility is: "+rep + "%")


fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 1, 1)

cat = categorical(cap_val);
bar(cat, volt_mean)
ylabel("Voltage [V]")
xlabel("Capacitance [pF]")
title("Mean")
subplot(2, 1, 2)
bar(cat, volt_std_perc)
title("Std ")
xlabel("Capacitance [pF]")
ylabel("Voltage in % [V]")
saveas(fig, "Figures/bar_cap_conditioning.png");
fontsize(30, "points");
%% Trasduction
cap_vec = zeros(11*5, 1);
volt_vec = zeros(11*5, 1);
counter = 1;
for i=1:1:11
    for j=1:1:5
        cap_vec(counter) = cap_val(i);
        volt_vec(counter) = volt(j, i);
        counter = counter + 1;
    end
end
%%

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:1:5*11
    plot(cap_vec(i), volt_vec(i), "x", "MarkerSize", 20, "LineWidth", 2.5, "HandleVisibility","off");
    hold on
end
% plot(cap_vec, volt_vec, "x", "MarkerSize", 10, "LineWidth", 5, "HandleVisibility","off");
linearModelThroughOrigin = fittype('a*x', 'independent', 'x', 'coefficients', 'a');
% fitResult = fit(cap_vec, volt_vec, 'poly1');
fitResult = fit(cap_vec, volt_vec, linearModelThroughOrigin);
hold on
plot(cap_vec, fitResult(cap_vec), "LineWidth", 2, "DisplayName", "Model", "Color", "red")
legend("show")
% Residual computation
residuals = zeros(5, 11);
for i=1:1:11
    residuals(:, i) = abs(volt(:, i) - fitResult(cap_val(i)));
end

std_res = 3*max(std(residuals));

% hold on
% plot(cap_val, fitResult(cap_val) + std_res, "-", "LineWidth", 2, "Color", "green", "DisplayName","Uncertainty band")
% hold on
% plot(cap_val, fitResult(cap_val) - std_res, "-", "LineWidth", 2, "Color", "green", "HandleVisibility", "off")

ylabel("Voltage (difference) [V]")
xlabel("Capacitance [pF]")
% title("Trasduction Diagram")
fontsize(30, "points")

%% Calibration Diagram

coeffs = coeffvalues(fitResult);
a = coeffs(1);
% b = coeffs(2);
% c = coeffs(3);

% inverseQuadratic = @(y_val) (-b + sqrt(b^2 - 4*a*(c - y_val))) / (2*a); %, ...
                             % (-b - sqrt(b^2 - 4*a*(c - y_val))) / (2*a)];
inverse = @(y) (y/a);
figure
for i=1:1:5*11
    plot(volt_vec(i), cap_vec(i), "x", "MarkerSize", 20, "LineWidth", 2.5, "HandleVisibility","off");
    hold on
end

hold on
plot(volt_vec, inverse(volt_vec), "LineWidth", 2, "DisplayName", "Model", "Color", "red")
legend("show")
hold on
plot(volt_vec, inverse(volt_vec) + std_res/a, "-", "LineWidth", 2, "Color", "green", "DisplayName","Uncertainty band")
hold on
plot(volt_vec, inverse(volt_vec) - std_res/a,  "-", "LineWidth", 2, "Color", "green", "HandleVisibility", "off")

% hold on
% plot(volt_vec, inverse(volt_vec + std_res), ".-", "MarkerSize", 25, "LineWidth", 2, "Color", "black", "DisplayName","Uncertainty band")
% hold on
% plot(volt_vec, inverse(volt_vec- std_res),  ".-", "MarkerSize", 25, "LineWidth", 2, "Color", "black", "HandleVisibility", "off")

% hold on
% plot(fitResult(cap_val) + std_res', cap_val, ".", "MarkerSize", 35, "Color", "black", "DisplayName","Uncertainty band")
xlabel("Voltage (difference) [V]")
ylabel("Capacitance [pF]")
% title("Calibration Diagram")
fontsize(30, "points")
