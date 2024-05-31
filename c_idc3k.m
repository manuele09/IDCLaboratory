% CAPACITANCE OF A 3-LAYERS IDC
% KIM MODEL

% INPUT PARAMETERS 
% l = overlapping finger length
% n = number of IDC finger pairs
% b= finger width;
% d= finger spacings;
% a= half finger width; 
% h1= layer 1 thickness;
% h2= layer 2 thickness;
% h3= layer 3 thickness;
% eps1= dielectric permittivity of layer 1;
% eps2= dielectric permittivity of layer 2;
% eps3= dielectric permittivity of layer 3;
% display = true if you want to display all 3 partial capacitances

% OUTPUT PARAMETERS 
% Cidc= Overall IDC Capacitance [F]

function Cidc=c_idc3k(eps1,eps2,eps3,h1,h2,h3,b,d,l,n, display)
        if ~exist('display','var')
            display=false;    
        end
        eps0=8.85*10^(-12); % dielectric permittivity of vacuum, in F/m;
        a=b/2;
        
        % Computation of C1
        F1=sinh(pi*b/(4*h1));
        G1=sinh(pi/(2*h1)*(b/2+d));
        H1=sinh(pi/(2*h1)*(b/2+d+a));

        k1_primo=F1/G1*sqrt((H1^2-G1^2)/(H1^2-F1^2));
        k1=sqrt(1-k1_primo^2);
        
        if (0.707<=k1)&&(k1<=1)
        Kk1_primo_Kk1=(2/pi)*(log(2*sqrt((1+k1)/(1-k1))));
        elseif (0<=k1)&&(k1<=0.707)
        Kk1_primo_Kk1=pi/2/(log(2*sqrt((1+k1_primo)/(1-k1_primo))));
        end 
        
        C1=2*eps0*eps1/Kk1_primo_Kk1; % [F]
       
        % Computation of C3 
        F3=sinh(pi*b/(4*(h2+h3)));
        G3=sinh(pi/(2*(h2+h3))*(b/2+d));
        H3=sinh(pi/(2*(h2+h3))*(b/2+d+a));
        k3_primo=F3/G3*sqrt((H3^2-G3^2)/(H3^2-F3^2));
        k3=sqrt(1-k3_primo^2);
        
        if (0.707<=k3)&&(k3<=1)
        Kk3_primo_Kk3=(2/pi)*(log(2*sqrt((1+k3)/(1-k3))));
        elseif (0<=k3)&&(k3<=0.707)
        Kk3_primo_Kk3=pi/2/(log(2*sqrt((1+k3_primo)/(1-k3_primo))));
        end 
        
        C3=2*eps0*eps3/Kk3_primo_Kk3; % [F]
          
        % Computation of C2        
        F2=sinh(pi*b/(4*h2));
        G2=sinh(pi/(2*h2)*(b/2+d));
        H2=sinh(pi/(2*h2)*(b/2+d+a));
        k2_primo=F2/G2*sqrt((H2^2-G2^2)/(H2^2-F2^2));
        k2=sqrt(1-k2_primo^2);

        if (0.707<=k2)&&(k2<=1)
        Kk2_primo_Kk2=(2/pi)*(log(2*sqrt((1+k2)/(1-k2))));
        elseif (0<=k2)&&(k2<=0.707)
        Kk2_primo_Kk2=pi/2/(log(2*sqrt((1+k2_primo)/(1-k2_primo))));
        end
            
        C2=2*eps0*(eps2-eps3)/Kk2_primo_Kk2; % [F]
        
        % Computation of Cidc     
        Cidc=n*l*(C1+C2+C3); % [F]
        if (display)
        disp("C1 Kim Model [F] = " + (C1*n*l)/1e-12 + "pF");
        disp("C2 Kim Model [F] = " + (C2*n*l)/1e-12 + "pF");
        disp("C3 Kim Model [F] = " + (C3*n*l)/1e-12 + "pF");
        end

        
end
        
