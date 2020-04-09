% Immune response, Reynolds et al. 2006

function f = immune (t,x,pars)

persistent b_learned  % flag to indicate whether immune system has learned about covid-19
if isempty(b_learned)  % first time this function is called
    b_learned = false;
end

k_pm   = 0.6;  % /M-units/h
k_mp   = 0.01; % /P-units/h
s_m    = 0.005; % M-units/h
mu_m   = 0.002; % /h
k_pg   = 0.23;   % /h, varies, 0.21-2.44/h
p_inf  = 20e6;  % /cc
k_pn   = 1.8;  % /N*-units/h
k_np   = 0.1;  % /P-units/h
k_nn   = 0.01; % /N*-units/h
s_nr   = 0.08; % NR-units/h
mu_nr  = 0.12; % /h
mu_n   = 0.05; % /h
k_nd   = 0.02; % /D-unit/h
k_dn   = 0.35; % D-units/h
x_dn   = 0.06; % N*-units
mu_d   = 0.02; % /h
s_c    = 0.0125; % CA-units/h
k_cn   = 0.04; % CA-units/h
k_cnd  = 48;  % N*-units/D-units
mu_c   = 0.1;  % /h

% HTN simulation, increase NR generation rate
s_nr = s_nr;

% new coupling constants with RAS
kAT1R = 0.0035;
kACE2 = 0.021;
kAng17= 0.010;
kAT2R = 0.0051;

% retrieve inputs
P  = x(1);
Ns = x(2);
D  = x(3);
cA = x(4);

at1rang2  = pars(1);
ace2      = pars(2);
ang17     = pars(3);
at2rang2  = pars(4);
covid19   = pars(5);  % covid-19 combining with ACE2m and entering cells
cell2plasma = pars(6);  % cells burst releasing covid-19 back to plasma

% change some of the parameters for covid simulations
% if D < 0.04 & not(b_learned) % immune system doesn't know to respond yet
% if D < 0.01 & not(b_learned) % immune system doesn't know to respond yet
if P < 2 & not(b_learned) % immune system doesn't know to respond yet
% if P < 1 & not(b_learned) % immune system doesn't know to respond yet
    k_pm   = 0.1*k_pm; 
    k_pn   = 0.1*k_pn; 
else  %  enough shit happens, better do something
    b_learned = true;
    % no change in k_pm and k_pn
    k_pm   = 1.0*k_pm; 
    k_pn   = 1.0*k_pn;     
%     k_pm   = 2*k_pm; 
%     k_pn   = 2*k_pn; 
end
% p_inf  = 1e10;   % tons of resources

fNs   = fCA(Ns,cA);
dPdt  = k_pg*P*(1-P/p_inf) - (k_pm*s_m*P)/(mu_m+k_mp*P) - k_pn*fNs*P;
  dPdt   = dPdt + covid19;  % covid-19 combining with ACE2m and entering cells
  dPdt   = dPdt - cell2plasma;
%   [k_pg*P*(1-P/p_inf), - (k_pm*s_m*P)/(mu_m+k_mp*P), - k_pn*fNs*P, covid19]
%   pause
R     = fCA(k_nn*Ns+k_np*P+k_nd*D,cA);
dNsdt = s_nr*R/(mu_nr+R) - mu_n*Ns;
  dNsdt  = dNsdt + kAT1R*at1rang2;
fs    = fNs^6/(x_dn^6+fNs^6);
dDdt  = k_dn*fs - mu_d*D;
fNsD  = fCA(Ns+k_cnd*D,cA);
dCAdt = s_c + k_cn*fNsD/(1+fNsD) - mu_c*cA;
  dCAdt  = dCAdt + (kACE2*ace2 + kAng17*ang17 + kAT2R*at2rang2);

% form return vector, first convert from hour to minute
f = [dPdt; dNsdt; dDdt; dCAdt]/60;