%% I/O Run: Generate data based on the physics-based model

clear all
clearvars
Ts = 1.0;
close all
Nsteps =40;
D=zeros(Nsteps,1); 
Pc = zeros(Nsteps,1);
Ug= zeros(Nsteps,1);

nsim=50;

Zin = []
Zout= []

for p=1:nsim

for i=1:Nsteps
    D(i)= (0.25 - 0.01).*rand(1,1) + 0.01;
    %D(i) = 0.2;
end

fub = 0.9;
for i=1:Nsteps
    Pc(i) = (fub-0.0).*rand(1,1) + 0.0;
    %Pc(i) = 0.5;
 
end

for i=1:Nsteps
     %Ug(i)= (100 - 40).*rand(1,1) + 40;
     ug(i) = 82.3;
end

%You will notice that CSTBR is 16 states..
%We do not need them all as the last ones are some "penalty"
%states that are used in the optimization problem of fba. 
N=21;
yk = zeros(Nsteps+1,4*N+3);
yksim = zeros(Nsteps+1,4*N+3);

%========================
%%Assign Initial Condition
yo = [];
eg = 0.309;                       % gas holdup, %
el = 1-eg;
tr = 320.15;                      % temperature, C
Area = 0.002436;                  % reactor cross-sectional area, m^2

L = 1.06;   % reactor length, m 
zs = L/(N-1);                     % spatial step size, m
pL = 1.013e5;                     % pressure at top of column, Pa
pc = Pc(1);                         % partial pressure of CO %(60+2*(mm-1))/100;
Hc = 8.0e-4;                      % Henry's constant for CO in water, mol/L*atm
pc2 = fub-Pc(1);                        % partial pressure of CO2
Hc2 = 2.5e-2;      
for i = 1:N
    po = pL + 1000*9.81*zs*(N-i)*el;   
    cgii = pc*po/8.314/tr;           
    c2gii = pc2*po/8.314/tr; 
    clsi = cgii*8.314*tr*Hc*1000/1.013e5;
    c2lsi = c2gii*8.314*tr*Hc2*1000/1.013e5;  
    yz = [cgii c2gii clsi c2lsi];
    yo = [yo yz];
 end
%y0 vector : co_g, c2_g, co_l, co2_l, biomass, ace, eth
yo = [yo, 0.3, 0, 0];
yk(1,:)=yo;
yksim(1,:) = yo; 
%=========================
save('init_cond.mat','yo')

%ykplusone = yk(1,:);
for j=1:Nsteps
Tin=(j-1)*Ts;
Tout=j*Ts;
tic
fprintf('Sampling Step %d',j)
[yk(j+1,:)] = BCR_Discrete(yk(j,:),D(j),Pc(j),Ug(j),Ts,fub);
[yksim(j+1,:)] = BCR_Discrete_Simulated(yk(j,:),D(j),Pc(j),Ug(j),Ts,fub);
Zin=[Zin ; [D(j),Pc(j),yksim(j,85),yksim(j,86),yksim(j,87) ]];
Zout=[Zout;[yksim(j,85),yksim(j,86),yksim(j,87)]];
toc

end

figure(1)
hold on
plot(yk(:,85),'r-o')
hold on
plot(yksim(:,85),'k--')

legend('X')

figure(2)
hold on
plot(yk(:,86),'g-o')
plot(yksim(:,86),'k--')

hold on
plot(yk(:,87),'b-o')
plot(yksim(:,87),'k--')

legend('A','E')

figure(3)
hold on
plot(yk(:,87)./yk(:,86),'r-o')
hold on
plot(yksim(:,87)./yksim(:,86),'k--')

legend('Selectivity')


Qmedia = L*Area*D;

n_acetate(p) = sum(Qmedia.*yk(1:Nsteps,86));
n_ethanol(p) = sum(Qmedia.*yk(1:Nsteps,87));

end

figure(9000)
plot(n_ethanol,'-o')
hold on
plot(n_acetate,'-o')
legend('E','A')

%save('Data.mat','Zin','Zout')


