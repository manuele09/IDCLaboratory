%% Section 2 of Report
clear
mkdir ./Figures/idcSimulation
%% Only PET Capacitance (Kim and Matlab model)
%Definition of IDC parameters (Table 1 and 2)
l = 5e-3;           %overlapping finger length
n = 6;              %number of IDC finger pairs
b = 0.3e-3;         %finger width;
d = 0.3e-3;         %finger spacings;

h1 = 140e-6;        %layer 1 thickness;
h2 = 10e-6;         %layer 2 thickness;
h3 = 10e-6;         %layer 3 thickness;
eps1 = 3.5;         %dielectric permittivity of PET (substrate);
eps2 = 1;           %dielectric permittivity of Polidopamina (sensitive layer);
eps3 = 1;           %dielectric permittivity of layer 3 (MUT);
lambda = 2*(b + d);


% Calculation of IDC Capacitance using Kim Model
CKim = c_idc3k(eps1,eps2,eps3,h1,h2,h3,b,d,l,n);
disp("Capacitance Only PET (KIM): " + CKim/1e-12 + " pF")

% Interdigital Matlab Model parameters (same params as above) 
idcObj = interdigitalCapacitor();
idcObj.NumFingers = n*2;
idcObj.FingerLength = l;
idcObj.FingerWidth = b;
idcObj.FingerSpacing = d;
% Additional parameters (Table 3)
idcObj.FingerEdgeGap = 0.5e-3;
idcObj.TerminalStripWidth = 1e-3;
idcObj.PortLineWidth = 1e-3;
idcObj.PortLineLength = 9e-3; % Ribbon Length
idcObj.Conductor = metal("silver");
idcObj.GroundPlaneWidth = 6.9e-3; % To neglect the Ground Plane
idcObj.Height = 1;

idcObj.Substrate = dielectric(Name="PET", EpsilonR=eps1, LossTangent=0.00100, Thickness=h1);
CModel = idcObj.capacitance(20e3);
disp("Capacitance Only PET (Matlab Model) = " + CModel/1e-12 + "pF");
disp("The mismatch is: " + string(abs(CKim-CModel)/1e-12) + "pF"); 
%% Capacitance with second layer (From now only Kim Model)
% Parameters in Table 4 (others are unchanged)
h2=100e-6;
eps2=1.23;

CKim = c_idc3k(eps1,eps2,eps3,h1,h2,h3,b,d,l,n);
disp("Capacitance With Second Layer: " + CKim/1e-12 + " pF")
%% Sensitivity with respect to h2 (Fig 2)
% h2 will vary from StartH to endH, with a step stepH.
startH = 3e-6;
stepH = 3e-6;
endH = lambda;

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for newH=startH:stepH:endH
    newC = c_idc3k(eps1,eps2,eps3,h1,newH,h3,b,d,l,n);
    plot(newH / 1e-6, newC/1e-12, ".-", "MarkerSize", 24);
    hold on
end
xlabel("h2 [um]");
ylabel("Capacitance [pF]");
% title("Capacitance Sensitivity analisys")
set(gca,'FontSize',40)
saveas(fig, "Figures/idcSimulation/h2Sensitivity.png");
%% Sensitivity with respect to h2 and eps2 (Fig 3)
startH =  3e-6; 
stepH = 10e-6;
endH = 200e-6;

startEps = 0.7;
stepEps = 0.2;
endEps = 4;

% Vector containing the permettivitys of the second layer
epsVector = zeros(floor((endEps-startEps)/stepEps), 1);
% Vector containing the thicknesses of the second layer
hVector = zeros(floor((endH-startH)/stepH), 1); 
% Vector containing the obtained Capacitance values (always overwritten)
capacitanceVector = zeros(floor((endEps-startEps)/stepEps), 1); 

fig = figure('units','normalized','outerposition',[0 0 1 1]);
j = 1;
for newH=startH:stepH:endH
    hVector(j) = newH; % Save for plotting
    i = 1;
    for newEps=startEps:stepEps:endEps
        capacitanceVector(i) = c_idc3k(eps1,newEps,eps3,h1,newH,h3,b,d,l,n);
        epsVector(i) = newEps; % Save for plotting
        i = i+1;
    end
    plot(epsVector, capacitanceVector/1e-12, ".-",...
        "LineWidth",2, "MarkerSize", 24, ...
        "DisplayName", strcat("h2= ",string(newH*1e3), " mm"));
    hold on
    j = j+1;
end
legend("show")
xlabel("Permittivity [F/m]");
ylabel("Capacitance [pF]");
% title("Capacitance Sensitivity analisys")
set(gca,'FontSize',40)
lgd = legend;
lgd.FontSize = 15; 
lgd.NumColumns = 8;
saveas(fig, "Figures/idcSimulation/h2eps2Sensitivity.png");