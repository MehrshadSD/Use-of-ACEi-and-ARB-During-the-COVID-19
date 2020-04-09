close all
clear all

% load and set initial conditions and parameters
load ss_soln
x0_all = x_ss;
% initialize plasma covid and the immune system model away from steady state

% below IC used with critical P = 2
% normotensive
%
% x0_all(12:16) = [1 0 0.0726 0.0071 0.5196];  % low exposure, recovery
x0_all(12:16) = [10 0 0.0726 0.0071 0.5196];  % high exposure, failure to clear virus, permanet tissue damage

% HTN, mild inflammation, s_nr x 2
% recovery even at high exposure, if no drug ?!
% dead with ACEi or ARB
% but if ACEi increase ACE2 by only 50%, it lives
% dead with ARB even if no ACE2 regulation, but if kAT1R x 5 to increase Ns
% response to (lower) AT1R-AngII it lives
%
% x0_all(12:16) = [1 0 0.0754 0.0084 0.5237];  % HTN, s_nr x 2

% HTN, high inflammation, s_nr x 5
% dead with no drug, even at low exposure
% recovery with ACEi or ARB, with or without ACE2 upregulation
%
% x0_all(12:16) = [1 0 0.0859 0.0137 0.5399];  % HTN, s_nr x 5 (dead w/o RAS inhibitors)

% drug simulations
sim_ACEi = false;   % ACE inhibition
sim_ARB  = true;  % ARB 
sim_ADAM = false;   % ADAM17 activator

pars_all = [pars_all, sim_ACEi, sim_ARB, sim_ADAM];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advance the combined model

[t,x] = ode15s(@(t,x) ras_ace2_immune(t,x,pars_all), [0 30*24*60], x0_all);
t = t/(24*60); % convert t from min to days

% retrieve and plot RAS variables
agt  = x(:,1);
prc  = x(:,2);
angI = x(:,3);
angII= x(:,4);
at1r = x(:,5);
at2r = x(:,6);
ang17= x(:,7);
angIV= x(:,8);

figure(1), clf
subplot(2,2,1), plot(t,angI);
ylabel('Ang I');
xlabel('days');
subplot(2,2,2), plot(t,angII);
ylabel('Ang II');
xlabel('days');
subplot(2,2,3), plot(t,at1r);
ylabel('AT1R-Ang II');
xlabel('days');
subplot(2,2,4), plot(t,at2r);
ylabel('AT2R-Ang II');
xlabel('days');

% retrieve and plot ACE2 variables
ace2m = x(:,9);
ace2p = x(:,10);
ace2i = x(:,11);
covid19p = x(:,12);

figure(2), clf
subplot(2,2,1), plot(t,ace2m);
ylabel('ACE2_m');
xlabel('days');
subplot(2,2,2), plot(t,ace2p);
ylabel('ACE_p');
xlabel('days');
subplot(2,2,3), plot(t,ace2i);
ylabel('ACE_i');
xlabel('days');
subplot(2,2,4), plot(t,covid19p);
ylabel('plasma covid-19');
xlabel('days');

% retrieve and plot immune variables
P  = x(:,13);
Ns = x(:,14);
D  = x(:,15);
cA = x(:,16);

figure(3), clf
subplot(2,2,1), plot(t,P);
ylabel('Pathogen (P)');
xlabel('days');
subplot(2,2,2), plot(t,Ns);
ylabel('Activated phagocytes (N^*)');
xlabel('days');
subplot(2,2,3), plot(t,D);
ylabel('Tissue damage (D)');
xlabel('days');
subplot(2,2,4), plot(t,cA);
ylabel('Anti-inflammatory mediator (C_A)');
xlabel('days');

% save the steady state solution and parameters
% x_ss = x(end,:);   pars_ss = pars_all;
% save 'ss_soln' x_ss pars_all
