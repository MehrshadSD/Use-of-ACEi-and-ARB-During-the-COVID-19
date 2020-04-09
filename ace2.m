% ACE2 dynamics

function [f covint] = ace2 (t,x,pars)

% retrieve parameters
covid190  = pars(1);
s_ace2m   = pars(2);
c_adam    = pars(3);
c_at1r    = pars(4);
at1rang2  = pars(5);
hace2     = pars(6);
cell2plasma = pars(7);  % cells burst releasing covid-19 to plasma
ace2pcovid19_rate = pars(8);

% retrieve variables
ace2m = x(1);
ace2p = x(2);
ace2i = x(3);
covp  = x(4);

% form the covid-19 coupling term, ACE2m, ACE2p, and plasma covid-19
Km = 1e9;
covid19 = covid190 * ace2m*covp/(covp+Km+1e4*ace2p);
% [covid190, ace2m,covp,(covp+Km+1e4*ace2p)]
% pause

% equations
dACE2m = s_ace2m - (c_adam + c_at1r*at1rang2)*ace2m - covid19 - (log(2)/hace2)*ace2m;
dACE2p = c_adam*ace2m - (log(2)/hace2)*ace2p;
dACE2i = c_at1r*at1rang2*ace2m + covid19 - (log(2)/hace2)*ace2i;

% membrane (ACE2m) to cell amount conversion, needs to be thought out
memb2cell = 1e3;
covint = memb2cell*covid19;  % amount of covid-19 entering cell (from plasma)

% plasma ACE2 binding to covid19
% ace2pcovid19_rate = 0.0001;  % made up 
ace2pcovid19 = ace2pcovid19_rate * ace2p * covp;

dcovp  = cell2plasma - covint - ace2pcovid19;
% [cell2plasma, - covint, - ace2pcovid19]

f = [dACE2m; dACE2p; dACE2i; dcovp];
% covid-19 internalization rate, this term will be used in the immune
% module, the scaling in front needs to be thought through
