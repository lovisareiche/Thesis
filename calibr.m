function calibr(d,bg,bd,sigEG, sigED,phig,phid)
% function to calibrate the model
% this allows me to use the same baseline mod file and change the
% calibration

% function saves all parameters of the model in param.mat
% input:    d: delta (in baseline 0.01)
%           piG: PI_Gbar (in baseline 0)
%           b: bias (in baseline 0)
%           sigE: standard error of the shock to expectations

beta=0.99;
eta=1;
psi_A=1;
psi_B=1;
epsilon_G=6;
epsilon_D=6;
phi_G=phig;
phi_D=phid;
phi_PI=1.2;
rhoA=0.9;
rhoSi=0.5;
rhoC=0.9;
PI_Gbar=0;
PI_Dbar= 0;

delta = d;
alpha=0.8;
bias_G =bg;
bias_D = bd;

sigma_epsE_G=sigEG;
sigma_epsE_D=sigED;

% calculate Rbar
Si=1;
E=1;
PI_G=PI_Gbar;
PI_D=PI_Dbar;
E_PI_G = (1+PI_G+bias_G)*E -1;
E_PI_D = (1+PI_D+bias_D)*E -1;

PIhash_G=(((1+PI_G)^(1-epsilon_G)-phi_G)/(1-phi_G))^(1/(1-epsilon_G)) -1;
PIhash_D=(((1+PI_D)^(1-epsilon_D)-phi_D)/(1-phi_D))^(1/(1-epsilon_D)) -1;

mc_G=((1+PIhash_G)/(1+PI_G)*(epsilon_G-1)/epsilon_G) * (1-phi_G*beta*(1+PI_G)^epsilon_G)/(1-phi_G*beta*(1+PI_G)^(epsilon_G-1));
mc_D=((1+PIhash_D)/(1+PI_D)*(epsilon_D-1)/epsilon_D) * (1-phi_D*beta*(1+PI_D)^epsilon_D)/(1-phi_D*beta*(1+PI_D)^(epsilon_D-1));

q=mc_G/mc_D;
E_q = q*(1+E_PI_D)/(1+E_PI_G);
Rbar = q*(1+E_PI_D)/(beta*E_q);

clear Si E PI_G PI_D E_PI_G PIhash_G PIhash_D mc_G mc_D q E_q

save param.mat
end