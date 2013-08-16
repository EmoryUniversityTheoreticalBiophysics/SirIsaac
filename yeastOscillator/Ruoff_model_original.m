
function [T,yintegrated,yDerivs,params] = Ruoff_model_original(t,y0,Temperature)

% BCD 8.19.2011 ******************** 
% (temperature dependent) rate constants
% k1,k2 from "Fit of model to Hemker et al. data" 
% at end of Fig.3 caption; rest from Table 1
% in mM/min (except for K1, in mM?)
x_temp = [ 2.5,
          25.,
          2.0,
          16.0,
          100.0,
          1.28,
          12.0,
          1.8,
          13.0,
          0.52];
Tref = 286.5; % K
R = 0.0083144; % kJ/K/mol

% Energy barriers
% Table 4 from the paper, in kJ/mol
EJ0   = 16.2;
Ek3   = 44.9;
Ek4   = 58.7;
Ek1   = 13.8;
Ek2   = 60.7;
Ek5   = 41.2;
Ek6   = 31.4;
Ekappa= 15.9;
Ek    = 24.3;
EK1   = 47.0;
Ebarriers = [ EJ0, 
             Ek1, 
             Ek2, 
             Ek3, 
             Ek4, 
             Ek5, 
             Ek6, 
             Ek, 
             Ekappa, 
             EK1 ];
prefactors = x_temp.*exp(Ebarriers/(R*Tref));

x_temp_new = prefactors.*exp(-Ebarriers/(R*Temperature));

% Rate constants
% Concentration units are in mM and time unit is minutes 
% Table 1 from the paper
J0    = x_temp_new(1);
k1    = x_temp_new(2); 
k2    = x_temp_new(3);
k3    = x_temp_new(4);
k4    = x_temp_new(5);  
k5    = x_temp_new(6); 
k6    = x_temp_new(7);
k     = x_temp_new(8); 
kappa = x_temp_new(9); 
q     = 4.0; 
K1    = x_temp_new(10);
psi   = 0.1; 

params = [J0,k1,k2,k3,k4,k5,k6,k,kappa,q,K1,psi];
% end BCD 8.19.2011 ********************

[ T, yintegrated ] = ode45( @Ruoff_model_derivs, t, y0 );

% 9.2.2011 also return derivatives at each timepoint
yDerivs = yintegrated;
s = size(yDerivs);
for i=[1:s(1)] % length(yDerivs)
    yDerivs(i,:) = Ruoff_model_derivs(0.,yintegrated(i,:));
end


%function dy = Ruoff_model_derivs_noTime(y)
%dy = Ruoff_model_derivs(0.,y)
%end
                   
                   

function dy = Ruoff_model_derivs(t,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the model described in P. Ruoff et al. Biophysical
% Chemistry 106 (2003) pages 179-192
% Author:  Abhishek Soni, CFDRC
% Date: October 23, 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%x_temp=[0.719,
%14.029,
%1.341,
%4.906,
%17.809,
%0.6635,
%5.361,
%0.5214,
%6.079,
%0.3258];

%% Pooled Variable totals %%
N     = 1.0;  %% NADH pool      %%
A     = 4.0;  %% Adenosine pool %%

%% Mapping concentrations %%
S1 = y(1);
S2 = y(2); 
S3 = y(3); 
S4 = y(4); 
N2 = y(5); 
A3 = y(6); 
S5 = y(7); 

%% Computing Individual Reaction Rates %%
v1 = k1*S1*A3/(1+(A3/K1)^q);
v2 = k2*S2*(N-N2); 
v3 = k3*S3*(A-A3); 
v4 = k4*S4*N2; 
v5 = k5*A3; 
v6 = k6*S2*N2;
v7 = k*S5; 
J  = kappa*(S4-S5);

%% Computing Species Rates %%
dS1dt = J0 - v1; 
dS2dt = 2*v1 - v2 - v6; 
dS3dt = v2 - v3; 
dS4dt = v3 - v4 - J; 
dN2dt = v2 - v4 - v6; 
dA3dt = -2*v1 + 2*v3 - v5; 
dS5dt = psi*J - v7 ;

%% Mapping onto output variable %%
dy = zeros(7,1); 
dy(1) = dS1dt; 
dy(2) = dS2dt;
dy(3) = dS3dt;
dy(4) = dS4dt;
dy(5) = dN2dt;
dy(6) = dA3dt;
dy(7) = dS5dt; 

end

end
