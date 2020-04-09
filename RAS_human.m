%RAS for the human
function f = RAS_human (t,x,pars)


%% Retrieve parameters by name.
sex_par          = pars( 1);
N_rs             = pars( 2);
X_PRCPRA         = pars( 3);
h_renin          = pars( 4);
h_AGT            = pars( 5);
h_AngI           = pars( 6);
h_AngII          = pars( 7);
h_Ang17          = pars( 8);
h_AngIV          = pars( 9);
h_AT1R           = pars(10);
h_AT2R           = pars(11);
k_AGT            = pars(12);
c_ACE            = pars(13);
c_Chym           = pars(14);
c_NEP            = pars(15);
c_ACE2           = pars(16);
c_IIIV           = pars(17);
c_AT1R           = pars(18);
c_AT2R           = pars(19);
AT1R_eq          = pars(20);
ACE2m            = pars(21);
kappa_ACEi       = pars(22);
kappa_ARB        = pars(23);

% ACE2-mediated reaction can now depend on membrane [ACE2]
c_ACE2 = c_ACE2 * ACE2m;

%% Retrieve variables by name.

AGT           = x(1); 
PRC           = x(2); 
AngI          = x(3); 
AngII         = x(4); 
AT1R          = x(5); 
AT2R          = x(6); 
Ang17         = x(7); 
AngIV         = x(8); 

% nu_AT1
nu_AT1 = ( 10^(0.0102 - 0.95 * log10(AT1R / AT1R_eq)) );
% R_sec
R_sec = ( N_rs * nu_AT1 );
% PRA
PRA = ( PRC * X_PRCPRA );

% AGT
AGT_p = ( k_AGT - PRA - log(2)/h_AGT * AGT );
% PRC
PRC_p = ( R_sec - log(2)/h_renin * PRC );
% AngI
AngI_p = ( PRA - ((1-kappa_ACEi) * c_ACE + c_Chym + c_NEP) * AngI - log(2)/h_AngI * AngI );
% AngII
AngII_p = ( ((1-kappa_ACEi) * c_ACE + c_Chym) * AngI - (c_ACE2 + c_IIIV + (1-kappa_ARB) * c_AT1R + c_AT2R) * AngII - log(2)/h_AngII * AngII );
% AT1R
AT1R_p = ( (1-kappa_ARB) * c_AT1R * AngII - log(2)/h_AT1R * AT1R );
% AT2R
AT2R_p = ( c_AT2R * AngII - log(2)/h_AT2R * AT2R );
% Ang17
Ang17_p = ( c_NEP * AngI + c_ACE2 * AngII - log(2)/h_Ang17 * Ang17 );
% AngIV
AngIV_p = ( c_IIIV * AngII - log(2)/h_AngIV * AngIV );

f = [AGT_p; PRC_p; AngI_p; AngII_p; AT1R_p; AT2R_p; Ang17_p; AngIV_p];









