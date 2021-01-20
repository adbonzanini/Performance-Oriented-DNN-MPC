function dy = bcr_model(t, y,condit,param)



klac = condit(1);
Hc = condit(2);
klac2 = condit(3);
Hc2 = condit(4);
tr = condit(5);
zs = condit(6);
N = condit(7);
ns = condit(8);
ug = condit(9);
ul = condit(10);
dl = condit(11);
eg = condit(12);
el = condit(13);
cgi = condit(14);
c2gi = condit(15);
Area = condit(16);
Qmedia = condit(17);
D = condit(18);

Amax = 0.15;
Emax = 0.2;
v_cm = param(1);
Kc   = param(2);
v_c2m = -0.1*param(3); %%
Kc2   = param(4);
Kic = param(5);

nmodel = N-1;
% Set other constants

Ma = 60.0/1000;       % Ma molecular weight g/mmol
Me = 46/1000;       % CO molecular weight g/mmol
Mbdo = 90/1000;     % 23BDO molecular weight g/mmol
Mlac = 90/1000;     % Lactate molecular weight g/mmol

% Define extracellular state variables

for i=1:N
    cg(i) = y(1+(i-1)*ns);   % gas CO concentration, mmol/L
    c2g(i) = y(2+(i-1)*ns);  % gas CO2 concentration, mmol/L
    cl(i) = y(3+(i-1)*ns);
    c2l(i) = y(4+(i-1)*ns);
end

X = y(N*ns+1);
A = y(N*ns+2);
E = y(N*ns+3);


%Define nodal rate terms 

ainh =max(0,(1-A/Amax));
einh =max(0,(1-E/Emax));

for i=1:N
 
co_term(i) = (cl(i))/(Kc + cl(i) + cl(i)^2/Kic);
co2_term(i) = c2l(i)/(Kc2 + c2l(i));

vc(i) = -v_cm*co_term(i)*co2_term(i); %adding minus due to cons. 
vc2(i) = -v_c2m*co_term(i)*co2_term(i);
end


% 
% for i = 1:nmodel
%     mu(i) = flux(i,1);      % growth rate, h-1
%     vc(i) = flux(i,2);      % CO uptake rate, mmol/g/h
%     vc2(i) = flux(i,3);     % CO2 uptake rate, mmol/g/h  
%     va(i) = flux(i,4);      % acetate synthesis rate, mmol/g/h 
%     ve(i) = flux(i,5);      % ethanol synthesis rate, mmol/g/h
%     vbdo(i) = flux(i,6);
%     vlac(i) = flux(i,7);
% end

%Some parameters
mu_max = 9e-1; %%
amax = 2.4; %%
emax = 0.8; %%

aco = mean(co_term);
aco2 = mean(co2_term);
avc  = mean(vc);
avc2 = mean(vc2);

mu = mu_max*aco2*aco;
vac = amax*aco*ainh;
vet = emax*aco*einh;


% pause(1.5)

%% Dynamics

dy = [];

i=1:N;
    
% Saturation gas concentrations
    
    cls(i) = cg(i)*8.314*tr*Hc/1.013e5*1000;
    c2ls(i) = c2g(i)*8.314*tr*Hc2/1.013e5*1000;

% Extracellular balances
        
i=1;
dcg(i) = 0; 
dc2g(i) = 0; 
cld(i) = (cl(i+1)-cl(i))/zs; %no flux
c2ld(i) = (c2l(i+1)-c2l(i))/zs; %no flux 
cld2(i) = (cl(i+1)-cl(i))/zs^2;  
c2ld2(i) = (c2l(i+1)-c2l(i))/zs^2;
dcl(i) = vc(i)*X + klac*(cls(i)-cl(i))/el-(ul-Qmedia/Area)*cld(i)/el+dl*cld2(i);    
dc2l(i) = vc2(i)*X + klac2*(c2ls(i)-c2l(i))/el-(ul-Qmedia/Area)*c2ld(i)/el+dl*c2ld2(i);

i=2;
cgd(i) = (cg(i)-cgi)/zs;
c2gd(i) = (c2g(i)-c2gi)/zs;  
dcg(i) = -klac*(cls(i)-cl(i))/eg-ug*cgd(i)/eg;   
dc2g(i) = -klac2*(c2ls(i)-c2l(i))/eg-ug*c2gd(i)/eg;  

cld(i) = (cl(i+1)-cl(i))/zs;
c2ld(i) = (c2l(i+1)-c2l(i))/zs;
cld2(i) = (cl(i+1)-2*cl(i)+cl(i-1))/zs^2;
c2ld2(i) = (c2l(i+1)-2*c2l(i)+c2l(i-1))/zs^2;
dcl(i) = vc(i)*X + klac*(cls(i)-cl(i))/el-(ul-Qmedia/Area)*cld(i)/el+dl*cld2(i);
dc2l(i) = vc2(i)*X + klac2*(c2ls(i)-c2l(i))/el-(ul-Qmedia/Area)*c2ld(i)/el+dl*c2ld2(i);

i=3;
cgd(i) = (2*cg(i+1)+3*cg(i)-6*cg(i-1)+cgi)/(6*zs);  
c2gd(i) = (2*c2g(i+1)+3*c2g(i)-6*c2g(i-1)+c2gi)/(6*zs);
dcg(i) = -klac*(cls(i)-cl(i))/eg-ug*cgd(i)/eg;    
dc2g(i) = -klac2*(c2ls(i)-c2l(i))/eg-ug*c2gd(i)/eg;

cld(i) = (cl(i+1)-cl(i))/zs;
c2ld(i) = (c2l(i+1)-c2l(i))/zs;
cld2(i) = (cl(i+1)-2*cl(i)+cl(i-1))/zs^2;
c2ld2(i) = (c2l(i+1)-2*c2l(i)+c2l(i-1))/zs^2;
dcl(i) = vc(i)*X + klac*(cls(i)-cl(i))/el-(ul-Qmedia/Area)*cld(i)/el+dl*cld2(i);
dc2l(i) = vc2(i)*X + klac2*(c2ls(i)-c2l(i))/el-(ul-Qmedia/Area)*c2ld(i)/el+dl*c2ld2(i);

i=4:N-2;
cgd(i) = (2*cg(i+1)+3*cg(i)-6*cg(i-1)+cg(i-2))/(6*zs);           
c2gd(i) = (2*c2g(i+1)+3*c2g(i)-6*c2g(i-1)+c2g(i-2))/(6*zs);
dcg(i) = -klac*(cls(i)-cl(i))/eg-ug*cgd(i)/eg;     
dc2g(i) = -klac2*(c2ls(i)-c2l(i))/eg-ug*c2gd(i)/eg;

cld(i) = (cl(i+1)-cl(i))/zs;
c2ld(i) = (c2l(i+1)-c2l(i))/zs;
cld2(i) = (cl(i+1)-2*cl(i)+cl(i-1))/zs^2;
c2ld2(i) = (c2l(i+1)-2*c2l(i)+c2l(i-1))/zs^2;
dcl(i) = vc(i)*X + klac*(cls(i)-cl(i))/el-(ul-Qmedia/Area)*cld(i)/el+dl*cld2(i);
dc2l(i) = vc2(i)*X + klac2*(c2ls(i)-c2l(i))/el-(ul-Qmedia/Area)*c2ld(i)/el+dl*c2ld2(i);

i=N-1;
cgd(i) = (2*cg(i+1)+3*cg(i)-6*cg(i-1)+cg(i-2))/(6*zs);           
c2gd(i) = (2*c2g(i+1)+3*c2g(i)-6*c2g(i-1)+c2g(i-2))/(6*zs);
dcg(i) = -klac*(cls(i)-cl(i))/eg-ug*cgd(i)/eg;     
dc2g(i) = -klac2*(c2ls(i)-c2l(i))/eg-ug*c2gd(i)/eg;
        
clr = cl(1)*abs(ul)*Area/(abs(ul)*Area+Qmedia);  
c2lr = c2l(1)*abs(ul)*Area/(abs(ul)*Area+Qmedia); 

cln = ((ul-Qmedia/Area)*zs*clr/el/dl-cl(i))/((ul-Qmedia/Area)*zs/el/dl-1);    
c2ln = ((ul-Qmedia/Area)*zs*c2lr/el/dl-c2l(i))/((ul-Qmedia/Area)*zs/el/dl-1); 
coeff=[cln, c2ln];

cld(i) = (cln-cl(i))/zs;
c2ld(i) = (c2ln-c2l(i))/zs;
cld2(i) = (cln-2*cl(i)+cl(i-1))/zs^2;
c2ld2(i) = (c2ln-2*c2l(i)+c2l(i-1))/zs^2;
dcl(i) = vc(i)*X + klac*(cls(i)-cl(i))/el-(ul-Qmedia/Area)*cld(i)/el+dl*cld2(i);
dc2l(i) = vc2(i)*X + klac2*(c2ls(i)-c2l(i))/el-(ul-Qmedia/Area)*c2ld(i)/el+dl*c2ld2(i);

i=N;
cgd(i) = (cg(i)-cg(i-1))/zs;            
c2gd(i) = (c2g(i)-c2g(i-1))/zs;
dcg(i) = -klac*(cls(i)-cln)/eg-ug*cgd(i)/eg;     
dc2g(i) = -klac2*(c2ls(i)-c2ln)/eg-ug*c2gd(i)/eg;

dcl(i) = 0;
dc2l(i) = 0;
      
i=1:N;
dyTmp = [dcg(i); dc2g(i); dcl(i); dc2l(i)];
dyTmp1 = reshape(dyTmp,[],1);
dy=[dy dyTmp1'];
 


dX =  mu*X - D*X;
dA =  Ma*vac*X - D*A;               
dE =  Me*vet*X - D*E; 

dy = [dy dX dA dE];

% % dylast=0;
% % for i=1:nmodel
% %     dylast = dylast + penalty(i);
% % end
% % dy = [dy dylast];
dy = dy';

end