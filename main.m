%Definition of IDC parameters
clear
l = 5e-3;           %overlapping finger length
n = 6;               %number of IDC finger pairs
b = 0.3e-3;         %finger width;
d = 0.3e-3;         %finger spacings;

h1 = 140e-6;        %layer 1 thickness;
h2 = 1e-5;          %layer 2 thickness;
h3 = 1e-5;          %layer 3 thickness;
eps1 = 3.5;         %dielectric permittivity of PET (substrate);
eps2 = 1.23;        %dielectric permittivity of Polidopamina (sensitive layer);
eps3 = 1;           %dielectric permittivity of layer 3 (MUT);
lambda = 2*(b + d);

%In the upper case the Second and Third Layer capacitance are
%zeros because of h2 and h3 (Read below).

%Notes: 
% -to neglect the Second layer (C2) put it's thickness to
% a value smaller than 1e-5;
% -to neglect the Third layer (C3) put the sum (h2+h3) to
% a value smaller than 2e-5;
% -C2 is always ZERO when eps2==eps1 (second layer fuses 
% with the third layer).

%% Calculation of IDC Capacitance using Kim Model

CKim = c_idc3k(eps1,eps2,eps3,h1,h2,h3,b,d,l,n, true);
disp("IDC Capacitance Only PET: " + CKim/1e-12 + " pF")

%% Calculation of IDC Capacitance using Interdigital Model (RF and Mixed Signal Matlab module)
idcObj = interdigitalCapacitor();
idcObj.NumFingers = n*2;
idcObj.FingerLength = l;
idcObj.FingerWidth = b;
idcObj.FingerSpacing = d;

idcObj.FingerEdgeGap = 0.5e-3;
idcObj.TerminalStripWidth = 1e-3;
idcObj.PortLineWidth = 1e-3;
idcObj.PortLineLength = 9e-3; % Ribbon Length
idcObj.Conductor = metal("silver");

% Trascuriamo il piano di massa
idcObj.GroundPlaneWidth = 6.9e-3; 
idcObj.Height = 1;

idcObj.Substrate = dielectric(Name="PET", EpsilonR=eps1, LossTangent=0.00100, Thickness=h1);
CIdc = idcObj.capacitance(20e3);
disp("IDC Capacitance Interdigital Model = " + CIdc/1e-12 + "pF");
%% Capacitance with second layer
h2=1e-3;
CKim = c_idc3k(eps1,eps2,eps3,h1,h2,h3,b,d,l,n, true);
disp("IDC Capacitance With Second Layer: " + CKim/1e-12 + " pF")
h2=1e-5; %Restore the previous value

%% Test al variare di h2
% Verifichiamo la saturazione se h2 supera lambda/2
figure

startH = 3e-6;
stepH = 3e-6;
endH = lambda;

for newH=startH:stepH:endH
    newC = c_idc3k(eps1,eps2,eps3,h1,newH,h3,b,d,l,n);
    plot(newH / 1e-6, newC/1e-12, ".-", "MarkerSize", 24);
    hold on
end
xlabel("h2 [um]");
ylabel("Capacitance [pF]");
% title("Capacitance Sensitivity analisys")
set(gca,'FontSize',40)

%Note:
% -h2 non può avere valori inferiori a 1e-6 perchè il
% modello utilizzato genera un eccezione;
% -per 3um<h2<15um il dispositivo è insensibile. Questo
% è dovuto al fatto che la capacità del secondo layer è nulla.
% -operativamente vediamo che per mantenere un andamento lineare
% l'altezza massima h2 deve essere 200um.

%% Test for different thickness and permittivity of the second layer
figure

startH =  3e-6; 
stepH = 10e-6;
endH = 200e-6;

startEps = 0.7;
stepEps = 0.2;
endEps = 4;

epsVector = zeros(floor((endEps-startEps)/stepEps), 1); % Vector containing the permettivitys of the second layer
hVector = zeros(floor((endH-startH)/stepH), 1); % Vector containing the thicknesses of the second layer
capacitanceVector = zeros(floor((endEps-startEps)/stepEps), 1); % Vector containing the obtained Capacitance values

j = 1;
for newH=startH:stepH:endH
    hVector(j) = newH; % Save for plotting
    i = 1;
    for newEps=startEps:stepEps:endEps
        capacitanceVector(i) = c_idc3k(eps1,newEps,eps3,h1,newH,h3,b,d,l,n);
        epsVector(i) = newEps; % Save for plotting
        i = i+1;
    end
    plot(epsVector, capacitanceVector/1e-12, ".-", "MarkerSize", 24, "DisplayName", strcat("h2= ",string(newH*1e3), " mm"));
    hold on
    j = j+1;
end
legend("show")
xlabel("Permittivity [F/m]");
ylabel("Capacitance [pF]");
% title("Capacitance Sensitivity analisys")
set(gca,'FontSize',40)
lgd = legend;
lgd.FontSize = 12; % Set the desired font size for the legend


%% Test al variare di h1
% Questo è solo un mio test personale: anche qui si osserva la saturazione
figure

startH = 3e-6;
stepH = 3e-6;
endH = lambda;

for newH=startH:stepH:endH
    newC = c_idc3k(eps1,eps2,eps3,newH,h2,h3,b,d,l,n);
    plot(newH / 1e-6, newC, ".-", "MarkerSize", 24);
    hold on
end
xlabel("h1 [um]");
ylabel("Capacitance [F]");
title("Capacitance Sensitivity analisys")
set(gca,'FontSize',14)

%% Test for different thickness and permittivity of the first layer
figure

startH =  3e-6; 
stepH = 10e-6;
endH = 200e-6;

startEps = 0.7;
stepEps = 0.2;
endEps = 4;

epsVector = zeros(floor((endEps-startEps)/stepEps), 1); % Vector containing the permettivitys of the second layer
hVector = zeros(floor((endH-startH)/stepH), 1); % Vector containing the thicknesses of the second layer
capacitanceVector = zeros(floor((endEps-startEps)/stepEps), 1); % Vector containing the obtained Capacitance values

j = 1;
for newH=startH:stepH:endH
    hVector(j) = newH; % Save for plotting
    i = 1;
    for newEps=startEps:stepEps:endEps
        capacitanceVector(i) = c_idc3k(newEps,eps2,eps3,newH,h2,h3,b,d,l,n);
        epsVector(i) = newEps; % Save for plotting
        i = i+1;
    end
    plot(epsVector, capacitanceVector/1e-12, ".-", "MarkerSize", 24, "DisplayName", strcat("h2= ",string(newH*1e3), " mm"));
    hold on
    j = j+1;
end
legend("show")
xlabel("Permittivity [F/m]");
ylabel("Capacitance [pF]");
% title("Capacitance Sensitivity analisys")
set(gca,'FontSize',40)
lgd = legend;
lgd.FontSize = 12; % Set the desired font size for the legend
