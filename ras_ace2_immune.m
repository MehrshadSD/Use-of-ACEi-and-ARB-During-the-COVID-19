function f = ras_ace2_immune (t,x,pars)

% get variables and parameters, last parameter is membrane [ACE2]
x_ras    = x(1:8);  
x_ace2   = x(9:12);
x_immune = x(13:end);

% pick out some that we need below
ace2m    = x_ace2(1);   % membrane ACE2
at1rang2 = x_ras(5);  % AT1R-Ang II
at2rang2 = x_ras(6);  % AT2R-Ang II
ang17    = x_ras(7);
covidc   = x_immune(1); % covid-19 in cells

% compute additional parameters

pars_ras    = pars(1:20);
pars_ace2   = pars(21:26);
pars_immune = pars(27:30);
sim_ACEi = pars(31);   % ACE inhibition
sim_ARB  = pars(32);  % ARB 
sim_ADAM = pars(33);   % ADAM17 activator

% pars_immune, when passed in, contains reference values
% replace by current values normalized by reference
at1rang2_ref = pars_immune(1);
ace2_ref     = pars_immune(2);
ang17_ref    = pars_immune(3);
at2rang2_ref = pars_immune(4);

% plasma ACE2 binding to covid19
ace2pcovid19_rate = 0.0001;  % made up 

pars_immune(1) = at1rang2/at1rang2_ref;
pars_immune(2) = ace2m/ace2_ref;
pars_immune(3) = ang17/ang17_ref;
pars_immune(4) = at2rang2/at2rang2_ref;

% % drug simulations
if sim_ACEi
    kappa_ACEi = 0.78;  % ACE inhibition
%     kappa_ACEi = 0.1;  % ACE inhibition
    kappa_ARB  = 0;
    ace_scaling= 5;
    pars_ace2(1) = pars_ace2(1)*ace_scaling;  % c_covid x5 as per ACEI, Ferrario et al. 2005, 
    pars_ras(16) = pars_ras(16)*ace_scaling;  % c_ACE2
    ace2pcovid19_rate = ace2pcovid19_rate * ace_scaling;  % plasma ACE2 also scaled proportionally
elseif sim_ARB
    kappa_ACEi = 0; 
    kappa_ARB = 0.67; % ARB
    ace_scaling= 3;
    pars_ace2(1) = pars_ace2(1)*ace_scaling;%3;  % ARB, Ferrario et al. 2005
    pars_ras(16) = pars_ras(16)*ace_scaling;%*3;  % c_ACE2
    ace2pcovid19_rate = ace2pcovid19_rate * ace_scaling;%3;  % plasma ACE2 also scaled proportionally
else
    kappa_ACEi = 0; % no drug
    kappa_ARB = 0; % no drug
    pars_ace2(1) = 1;  % determine ACE2-covid interactions, affected by ACE2 upregulation
end

pars_ras = [pars_ras ace2m kappa_ACEi kappa_ARB];
f_ras = RAS_human(t,x_ras,pars_ras);

% rate at which cells burst, releasing covid-19 back to plasma
% rate constant made up
cell2plasma = 1e-5 * covidc;

if sim_ADAM
    c_adam = pars_ace2(3);
    c_adam = c_adam * 40;
    pars_ace2(3) = c_adam;
end

pars_ace2 = [pars_ace2 cell2plasma ace2pcovid19_rate];
pars_ace2(5) = x_ras(5);  % use current AT1R-Ang II value
[f_ace2 covint] = ace2(t,x_ace2,pars_ace2);

pars_immune = [pars_immune covint cell2plasma];
f_immune = immune(t,x_immune,pars_immune);

f = [f_ras; f_ace2; f_immune];