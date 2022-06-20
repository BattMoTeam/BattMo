mrstModule add ad-core

an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';

% jsonstruct = parseBattmoJson('ProtonicMembrane/protonmembrane.json');
paramobj = ProtonicMembraneCellInputParams([]);

%%  Constants

R    = 8.3145;         % gas constant
F    = 96485;          % Faraday constant
d    = 30*micro*meter; % thickness in micron
E_0  = 1.253;          % standard potential
Na   = 6.02e23;        % Avogadro
kb   = 8.62e-5;        % Boltzmann (eV/K)
kb_J = 1.38e-23;       % Boltzmann (eV/K)
e_c  = 1.602e-19;

%%  Input parameters - Operating conditions

t_p_O2   = 0.5;          % tr.nr holes in 1bar humidified oxygen
flow_H2O = 100;          % mL/min
T        = 600 + 273.15; % K

f = F./(R.*T); % exp-faktor i kinetikk-uttrykk

pH2O_in = 3;

SU       = 0.2;
pH2O     = pH2O_in*(1 - SU);
pO2      = 1;
pH2O_neg = 0.027;
pH2      = pH2O_in;
pO2_flow = (pO2/(pH2O))*flow_H2O;

%%  Input parameters - Transport parameters

R_ct     = 0.2; % Ohm*cm^2
i_lc     = 10;  % limiting current cathodic
i_la_0   = -2;  % limiting current anodic
i_la     = i_la_0.*pH2O;
sigma_n0 = 1e-5;

%%

E_a_OCV = E_0 - 0.00024516.*T - R.*T./(2*F)*log(pH2O./(pO2^(1/2)));
E_c_OCV = -(1./f).*log(pH2); % cathode half-cell potential

Y       = 0.16;
Y_mol   = Y/((4.22e-8)^3*Na);
dH_hyd  = -115;       % kJ/mol
dS_hyd  = -119;       % J / mol / K
D0_prot = 0.021*38/T; % pre-exp proton diffusion
Ea_prot = 48;         % Activation energy proton diffusion
D_prot  = D0_prot*exp(-(Ea_prot*1000)/R/T);
K_hyd   = exp((dS_hyd/R)-dH_hyd*1000/(T*R)); % equation (8) in reference [1]
K_H     = K_hyd*pH2O;
K_H_neg = K_hyd*pH2O_neg;
OH_pos  = ((3*K_H - sqrt(K_H*(9*K_H-6*K_H*Y+K_H*Y^2+24*Y-4*Y^2)))/(K_H-4));  % Per formula unit
OH_neg  = ((3*K_H_neg - sqrt(K_H_neg*(9*K_H_neg-6*K_H_neg*Y+K_H_neg*Y^2+24*Y-4*Y^2)))/(K_H_neg-4));

sigma_prot = (F*OH_pos*D_prot); % Proton conductivity

K_H_p = K_hyd*0.027;
p_ref = ((3*K_H_p - sqrt(K_H_p*(9*K_H_p-6*K_H_p*Y+K_H_p*Y^2+24*Y-4*Y^2)))/(K_H_p-4));

sigma_p0 = (t_p_O2/(1 - t_p_O2))*(F*p_ref*D_prot);


return


couplingTerms = {};

coupterm = couplingTerm('Anode-Electrolyte', {an, elyte});
coupterm.couplingcells = [1, 1];
coupterm.couplingfaces = [1, 1];
couplingTerms{end + 1} = coupterm;

coupterm = couplingTerm('Cathode-Electrolyte', {an, elyte});
coupterm.couplingcells = [1, G.cells.num];
coupterm.couplingfaces = [1, G.faces.num];
couplingTerms{end + 1} = coupterm;

paramobj.couplingTerms = couplingTerms;

model = ProtonicMembraneCell(paramobj);



mols_to_mLmin = R*10*273.15*60;


