function [xkplusone] = BCR_Discrete(xk,Dk,pck,ugk,Ts,fub);



nmodel = 20; %Number of internal nodes
tspan = [0,Ts]; %Time of operation
N = nmodel+1;     % number of reactor discretization points
ns = 4;   %species with space variation (COl,CO2l, COg,CO2g)
nsvar = N*ns;
%% Operation parameters + others 

% Set operating conditions
L = 1.06;                         % reactor length, m 
ug = ugk ;                        % constant gas velocity, m/h
ul = -158.5;                      % Liquid recyle rate, m/h
D = Dk; %0.05                         % dillution rate,1/h
eg = 0.309;                       % gas holdup, %
Area = 0.002436;                  % reactor cross-sectional area, m^2
zs = L/(N-1);                     % spatial step size, m
dl = 4.5;                         % liquid dispersion coefficient, 0.25 m^2/h
tr = 320.15;                      % temperature, C
pL = 1.013e5;                     % pressure at top of column, Pa
pc = pck;                         % partial pressure of CO %(60+2*(mm-1))/100;
pc2 = (fub-pck);                        % partial pressure of CO2

Hc = 8.0e-4;                      % Henry's constant for CO in water, mol/L*atm
Hc2 = 2.5e-2;                     % Henry's constant for CO2 in water, mol/L*atm
Qmedia = L*Area*D;

% Set mass transfer coefficients
klac = 523;                       % CO gas-liquid mass transfer coefficient, h-1
klac2 = klac;                     % CO2 gas-liquid mass transfer coefficient, h-1

% Calculate gas and liquid holdups
el = 1-eg;
po = pL+1000*9.81*L*el;
cgi = pc*po/8.314/tr;           % CO feed concentration, mol/m3 = mmol/L
c2gi = pc2*po/8.314/tr;         % CO2 feed concentration, mol/m3 = mmol/L


%{
% Set operating conditions
L = 10;                         % reactor length, m 
Area = 3;                       % reactor cross-sectional area, m^2
dR = sqrt(4*Area/pi);           % reactor diameter, m  (0.6308 m)
zs = L/(N-1);                   % spatial step size, m
dl = 0.25;                      % liquid dispersion coefficient, m^2/h
tr = 310.15;                    % temperature, K
pL = 1.013e5;                   % pressure at top of column, Pa
pc = 0.5;   %or:0.5                    % partial pressure of CO %(60+2*(mm-1))/100;
Hc = 8.0e-4;                    % Henry's constant for CO in water, mol/(L*atm)
ph = 0.2;                       % partial pressure of H2
Hh = 6.6e-4;                    % Henry's constant for H2 in water, mol/(L*atm) (NIST)
pn = 0.3;   %or:0.3                     % partial pressure of N2
pc2 = 0;                        % partial pressure of CO2
Hc2 = 2.5e-2;                   % Henry's constant for CO2 in water, mol/(L*atm)
g = 9.81;                       % gravitation, m/s^2
densityL = 993.34;              % density of liquid phase (roughly water at 310K), kg/m^3
viscosityL = 0.0009242;         % viscosity of liquid phase (roughly water at T = 310K P = 100-300 kPa), Pa*s
D = 0.03;                       % dilution rate
Qmedia = D*L*Area;              % volumetric flowrate
p_atm = 1;                      % 1 atmospheric pressure

% Feed gas composition (inlet condition)
ul = -50;                       % superficial liquid velocity
ug = 101.8090;                   % superficial gas velocity
db = 1.2251;                    % mm
ub = 0.33*(g^0.76)*((densityL/viscosityL)^0.52)*(((db/2)*1e-3).^1.28);
eg = 0.1246;                    % gas phase holdup
el = 1-eg;                      % liquid phase holdup
p0 = pL+1000*9.81*L*el;         % Pa
cgi = pc*p0/8.314/tr;           % CO feed concentration, mol/m3 = mmol/L
hgi = ph*p0/8.314/tr;           % H2 feed concentration, mol/m3 = mmol/L
c2gi = pc2*p0/8.314/tr;         % CO2 feed concentration, mol/m3 = mmol/L

% Set mass transfer coefficients
klac = (1e-4)*3600;              % CO gas-liquid mass transfer coefficient, m/h  
klah = klac;                      % H2 gas-liquid mass transfer coefficient, m/h
klac2 = klac;                     % CO2 gas-liquid mass transfer coefficient, m/h
klo= 0.54;                      % O2 gas-liquid mass transfer coeff m/h
%}

% Form condition vector
condit = [klac,Hc,klac2,Hc2,tr,zs,N,ns,ug,ul,dl,eg,el,cgi,c2gi,Area,Qmedia,D];

% Set uptake parameters
vcm = 40; %35  
Kmc = 0.4; %0.2
Kic = 2.2; %1.6
param = [
     vcm                 % maximum CO uptake rate, mmol/g/h; 
     Kmc                 % CO uptake saturation constant, mmol/L;
     vcm                 % maximum CO2 uptake rate, mmol/g/h 
     Kmc                 % CO2 uptake saturation constant, mmol/L;                 
     Kic                 % CO inhibition constant, mmol/L
     ];

 yo = [];
 

po = pL + 1000*9.81*zs*(N-i)*el;   
%     cgii = pc*po/8.314/tr;           
%     c2gii = pc2*po/8.314/tr; 
%     clsi = cgii*8.314*tr*Hc*1000/1.013e5;
%     c2lsi = c2gii*8.314*tr*Hc2*1000/1.013e5;  
cgii = xk(1:N);
c2gii = xk(N+1:2*N);
clsi = xk(2*N+1:3*N);
c2lsi = xk(3*N+1:4*N);


yz = [cgii c2gii clsi c2lsi];
yo = [yo yz];

size(yo)

%y0 vector : co_g, c2_g, co_l, co2_l, biomass, ace, eth
yo = [yo, xk(nsvar+1),xk(nsvar+2),xk(nsvar+3)];


[T,X] = ode15s(@(T,X) bcr_model(T,X,condit,param),tspan,yo); 
% zn = linspace(0,L,N)

xkplusone = X(end,:);
% figure(1)
% hold on
% plot(T, X(:,nsvar+1))
% plot(T, X(:,nsvar+2))
% plot(T, X(:,nsvar+3))
% legend('Biomass','Acetate','EtOH')