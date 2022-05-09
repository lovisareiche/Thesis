%-------------------------------------------------------------------------
% Baseline mod file for durable goods NK
%-------------------------------------------------------------------------

% This file implements the baseline durable goods model with two 
% sectors. The model features log-utility, cobb douglas aggregation and the
% assumption that A and B optimise together (as one agent). There is a wage for 
% each sector but the wages must equal to induce work in both sectors.

@#define oneshock=1 % must manually select which model


% 1. Specify variables and parameters
%-------------------------------------------------------------------------

var 
@# if oneshock==1
    A       (long_name='technology')
@# elseif oneshock==0
    A_G     (long_name='technology G')
    A_D     (long_name='technology D')
@# endif
    Si      (long_name='Interest rate shock')
@# if oneshock==1
    C       (long_name='Cost push shock')
@# elseif oneshock==0
    C_G     (long_name='Cost push shock in G')
    C_D     (long_name='Cost push shock in D')
@# endif
    R       (long_name='Nominal Interest rate')
    Y_G     (long_name='aggregate output of groceries')
    Y_D     (long_name='aggregate output of durables')
    mc_G    (long_name='real marginal costs evaluated at grocery prices')
    mc_D    (long_name='real marginal costs evaluated at durable prices')
    w_G     (long_name='real wage in G in grocery prices')
    w_D     (long_name='real wage in D in grocery prices')
    PI_G    (long_name='Inflation of groceries')
    PI_D    (long_name='Inflation of durables')
    N_A     (long_name='labour supply of member A')
    N_B     (long_name='labour supply of member B')
    q       (long_name='relation durables prices to grocery prices')
    X_G     (long_name='Output gap groceries')
    X_D     (long_name='Output gap durables')
    GDP     (long_name='GDP aggregate')
    PI      (long_name='Price index')

    G       (long_name='consumption of groceries')
    D       (long_name='consumption of durables')
    vp_G    (long_name='price dispersion term of groceries')
    vp_D    (long_name='price dispersion term of durables')
    PIhash_G(long_name='reset price inflation for groceries')
    PIhash_D(long_name='reset price inflation for durables')
    x1_G    (long_name='variable 1 for recursive formulation of price setting of groceries')
    x2_G    (long_name='variable 2 for recursive formulation of price setting of groceries')
    x1_D    (long_name='variable 1 for recursive formulation of price setting of durables')
    x2_D    (long_name='variable 2 for recursive formulation of price setting of durables')
    Y_G_f   (long_name='flexible price output of groceries')
    Y_D_f   (long_name='flexible price output of durables')
    U_G     (long_name='marginal utility groceries')
    U_D     (long_name='marginal utility durables')
    ;
    
varexo 
    epsSi   (long_name='Innovation interest rate')
@# if oneshock==1
    epsA     (long_name='technology shock')
@# elseif oneshock==0
    epsA_G   (long_name='technology shock G')
    epsA_D   (long_name='technology shock D')
@# endif
@# if oneshock==1
    epsC     (long_name='cost push shock')
@# elseif oneshock==0
    epsC_G   (long_name='cost push shock G')
    epsC_D   (long_name='cost push shock D')
@# endif
    ;


parameters beta         (long_name='discount factor')
           alpha        (long_name='relative weight of two goods in utility')
           delta        (long_name='depreciation rate of durable goods')
           epsilon_G    (long_name='elasticity of substitution between goods varieties of groceries')
           epsilon_D    (long_name='elasticity of substitution between goods varieties of durables')
           eta          (long_name='labor disutility parameter')
           psi_A        (long_name='relative weight of labour disutility agent A')
           psi_B        (long_name='relative weight of labour disutility agent B')
           phi_G        (long_name='Calvo parameter groceries')
           phi_D        (long_name='Calvo parameter durables')
           phi_PI       (long_name='policy coefficient in Taylor rule')
           rhoA         (long_name='autocorrelation technology')
           rhoSi        (long_name='autocorrelation cost push shock')
           rhoC         (long_name='autocorrelation interest rate shock')
           PI_Gbar      (long_name='steady state groceries inflation')
           PI_Dbar      (long_name='steady state durables inflation')
           Rbar         (long_name='Steady state nominal interest')
           ;

load param;

set_param_value('alpha',alpha);
set_param_value('beta',beta);
set_param_value('delta',delta);
set_param_value('epsilon_D',epsilon_D);
set_param_value('epsilon_G',epsilon_G);
set_param_value('eta',eta);
set_param_value('phi_D',phi_D);
set_param_value('phi_G',phi_G);
set_param_value('phi_PI',phi_PI);
set_param_value('PI_Dbar',PI_Dbar);
set_param_value('PI_Gbar',PI_Gbar);
set_param_value('psi_A',psi_A);
set_param_value('psi_B',psi_B);
set_param_value('rhoA',rhoA);
set_param_value('rhoSi',rhoSi);
set_param_value('rhoC',rhoC);
Rbar = (1+PI_Gbar)/(beta);


% 2. Model equations
%-------------------------------------------------------------------------

model; 

[name='marg utility G']
U_G=alpha/G;
[name='marg utility D']
U_D=(1-alpha)/D;
[name='Euler durables groceries']
U_D = U_G*q -beta*(1-delta)*U_G(+1)*q(+1);
[name='Euler labour']
psi_A*N_A^eta /w_G = U_G;
psi_B*N_B^eta /w_D = U_G;
[name='intratemporal Euler']
U_G = beta*R/(1+PI_G(+1)) * U_G(+1);

@# if oneshock==1
    [name='real marginal costs groceries']
    mc_G=w_G/A;
    [name='real marginal costs durables']
    mc_D=w_D/(q*A);
@# elseif oneshock==0
    [name='real marginal costs groceries']
    mc_G=w_G/A_G;
    [name='real marginal costs durables']
    mc_D=w_D/(q*A_D);
@# endif

[name='market clearing']
G=Y_G;
(D-(1-delta)*D(-1)) = Y_D;

@# if oneshock==1
    [name='aggregate production of groceries']
    Y_G=(A*N_A)/(vp_G);
    [name='aggregate production of durables']
    Y_D=(A*N_B)/(vp_D);
@# elseif oneshock==0
    [name='aggregate production of groceries']
    Y_G=(A_G*N_A)/(vp_G);
    [name='aggregate production of durables']
    Y_D=(A_D*N_B)/(vp_D);
@# endif

[name='price dispersion law of motion groceries']
vp_G=(1-phi_G)*((1+PIhash_G)^(-epsilon_G))*((1+PI_G)^(epsilon_G)) + ((1+PI_G)^(epsilon_G))*phi_G*vp_G(-1);
[name='price dispersion law of motion durables']
vp_D=(1-phi_D)*((1+PIhash_D)^(-epsilon_D))*((1+PI_D)^(epsilon_D)) + ((1+PI_D)^(epsilon_D))*phi_D*vp_D(-1);
[name='inflation groceries']
(1+PI_G)^(1-epsilon_G)=(1-phi_G)*(1+PIhash_G)^(1-epsilon_G) + phi_G;
[name='inflation durables']
(1+PI_D)^(1-epsilon_D)=(1-phi_D)*(1+PIhash_D)^(1-epsilon_D) + phi_D;

@# if oneshock==1
[name='reset inflation groceries']
(1+PIhash_G)=(epsilon_G/(epsilon_G-1))*C*(x1_G/x2_G)*(1+PI_G);
[name='reset inflation durables']
(1+PIhash_D)=(epsilon_D/(epsilon_D-1))*C*(x1_D/x2_D)*(1+PI_D);
@# elseif oneshock==0
[name='reset inflation groceries']
(1+PIhash_G)=(epsilon_G/(epsilon_G-1))*C_G*(x1_G/x2_G)*(1+PI_G);
[name='reset inflation durables']
(1+PIhash_D)=(epsilon_D/(epsilon_D-1))*C_D*(x1_D/x2_D)*(1+PI_D);
@# endif

[name='variable 1 law of motion groceries']
x1_G=mc_G*Y_G*U_D^(-1) + phi_G*beta*(1+PI_G(+1))^(epsilon_G)*x1_G(+1);
[name='variable 2 law of motion groceries']
x2_G=Y_G*U_D^(-1) + phi_G*beta*(1+PI_G(+1))^(epsilon_G)*x2_G(+1);
[name='variable 1 law of motion durables']
x1_D=mc_D*Y_D*U_D^(-1) + phi_D*beta*(1+PI_D(+1))^(epsilon_D)*x1_D(+1);
[name='variable 2 law of motion durables']
x2_D=Y_D*U_D^(-1) + phi_D*beta*(1+PI_D(+1))^(epsilon_D)*x2_D(+1);

[name='Taylor rule']
R/Rbar = (((1+PI_G)^alpha)*((1+PI_D)^(1-alpha))/((1+PI_Gbar)^alpha)*((1+PI_Dbar)^(1-alpha)))^(phi_PI) * Si;

[name='total GDP']
GDP=Y_G+Y_D;

[name='inflation index']
1+PI=(1+PI_G)^alpha*(1+PI_D)^(1-alpha);

@# if oneshock==1
    [name='Technology']
    log(A)=rhoA*log(A(-1))+epsA; 
@# elseif oneshock==0
    [name='Technology G']
    log(A_G)=rhoA*log(A_G(-1))+epsA_G; 
    [name='Technology D']
    log(A_D)=rhoA*log(A_D(-1))+epsA_D; 
@# endif
[name='Interest rate shock']
log(Si)=rhoSi*log(Si(-1))+epsSi; 
@# if oneshock==1
    [name='Technology']
    log(C)=rhoC*log(C(-1))+epsC; 
@# elseif oneshock==0
    [name='Cost push shock G']
    log(C_G)=rhoC*log(C_G(-1))+epsC_G; 
    [name='Cost push shock D']
    log(C_D)=rhoC*log(C_D(-1))+epsC_D; 
@# endif


[name='definition of q']
q=q(-1)*(1+PI_D)/(1+PI_G);


@# if oneshock==1
    [name='flexible output']
    Y_G_f = (alpha/psi_A * (epsilon_G-1)/epsilon_G * A)^(1/(eta+1)) * A; 
    Y_D_f = ((1-alpha)/(psi_B*delta) * (epsilon_G-1)/epsilon_G * A * 1/((epsilon_G-1)/epsilon_G*epsilon_D/(epsilon_D-1)) * 1/(1-(1-delta)*(1+PI_D(+1))/R))^(1/(eta+1))*A;
@# elseif oneshock==0
    [name='flexible output']
    Y_G_f = (alpha/psi_A * (epsilon_G-1)/epsilon_G * A_G)^(1/(eta+1))*(A_G); 
    Y_D_f = ((1-alpha)/(psi_B*delta) * (epsilon_G-1)/epsilon_G * A_G * 1/((epsilon_G-1)/epsilon_G*A_G/A_D*epsilon_D/(epsilon_D-1)) * 1/(1-(1-delta)*(1+PI_D(+1))/R))^(1/(eta+1))*A_D;
@# endif

[name='output gap']
X_G = Y_G/Y_G_f;
X_D = Y_D/Y_D_f;

end;


% 3. Shocks
%-------------------------------------------------------------------------

shocks;
@# if oneshock==1
    var epsA; stderr 0.01;
@# elseif oneshock==0
    var epsA_G; stderr 0.01;
    var epsA_D; stderr 0.01;
@# endif
var epsSi; stderr 0.01;
@# if oneshock==1
    var epsC; stderr 0.001;
@# elseif oneshock==0
    var epsC_G; stderr 0.001;
    var epsC_D; stderr 0.001;
@# endif


end;


% 4. Steady state caculation
%-------------------------------------------------------------------------

steady_state_model;

@# if oneshock==1
    A=1;
    C=1;
@# elseif oneshock==0
    A_G=1;
    A_D=1;
    C_G=1;
    C_D=1;
@# endif
Si=1;
PI_G=PI_Gbar;
PI_D=PI_Dbar;

PIhash_G=(((1+PI_G)^(1-epsilon_G)-phi_G)/(1-phi_G))^(1/(1-epsilon_G)) -1;
PIhash_D=(((1+PI_D)^(1-epsilon_D)-phi_D)/(1-phi_D))^(1/(1-epsilon_D)) -1;

vp_G=((1-phi_G)*((1+PIhash_G)^(epsilon_G)))/((1+PI_G)^epsilon_G*(1-((1+PI_G)^(epsilon_G))*phi_G));
vp_D=((1-phi_D)*((1+PIhash_D)^(epsilon_D)))/((1+PI_D)^epsilon_D*(1-((1+PI_D)^(epsilon_D))*phi_D));

mc_G=((1+PIhash_G)/(1+PI_G)*(epsilon_G-1)/epsilon_G) * (1-phi_G*beta*(1+PI_G)^epsilon_G)/(1-phi_G*beta*(1+PI_G)^(epsilon_G-1));
mc_D=((1+PIhash_D)/(1+PI_D)*(epsilon_D-1)/epsilon_D) * (1-phi_D*beta*(1+PI_D)^epsilon_D)/(1-phi_D*beta*(1+PI_D)^(epsilon_D-1));

w_G=mc_G;
R = (1+PI_G)/(beta);

Y_G=(alpha*w_G/(psi_A))^(1/(1+eta))*(1/vp_G)^(eta/(1+eta));
N_A=(alpha*w_G/(Y_G*psi_A))^(1/eta);

w_D = (psi_B*Y_G/alpha *(vp_D * Y_G*mc_D/alpha * delta*(1-alpha)/(1-beta*(1-delta)))^eta)^(1/(1+eta));
q = w_D/mc_D;

Y_D=Y_G/(q*alpha) * delta*(1-alpha)/(1-beta*(1-delta));
N_B=Y_D*vp_D;

Y_G_f = (alpha/psi_A * (epsilon_G-1)/epsilon_G)^(1/(eta+1)); 
//Y_D_f = ((1-alpha)/psi_B * (epsilon_D-1)/epsilon_D * ((epsilon_G-1)/epsilon_G)/((epsilon_D-1)/epsilon_D))^(1/(eta+1)); 
Y_D_f = ((1-alpha)/(psi_B*delta) * ((epsilon_G-1)/epsilon_G) * 1/((epsilon_G-1)/epsilon_G*epsilon_D/(epsilon_D-1)) * 1/(1-beta*(1-delta)))^(1/(eta+1)); 
X_G = Y_G/Y_G_f;
X_D = Y_D/Y_D_f;

G=Y_G;
D = Y_D/(1-(1-delta));

U_G=alpha/G;
U_D=(1-alpha)/D;

x1_G=mc_G*Y_G*U_D^(-1) /(1- phi_G*beta*(1+PI_G)^(epsilon_G));
x2_G=Y_G*U_D^(-1) /(1- phi_G*beta*(1+PI_G)^(epsilon_G));
x1_D=mc_D*Y_D*U_D^(-1) /(1- phi_D*beta*(1+PI_D)^(epsilon_D));
x2_D=Y_D*U_D^(-1) /(1- phi_D*beta*(1+PI_D)^(epsilon_D));

GDP=Y_G+Y_D;
PI=(1+PI_G)^alpha*(1+PI_D)^(1-alpha)-1;


end;





steady;

//stoch_simul(order=1,irf=40,periods=200,nograph);
stoch_simul(order=1,irf=40,nograph);



