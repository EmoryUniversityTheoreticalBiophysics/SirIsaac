%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the model described in P. Ruoff et al. Biophysical
% Chemistry 106 (2003) pages 179-192
% Author:  Abhishek Soni, CFDRC
% Date: October 23, 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setting the time range %%
%t = [0:0.01:100];
t = [0:.01:25]; % BCD 8.19.2011

% Set initial conditions as per Table 2 in the paper. 
y0 = zeros(7,1) ; 
y0(1) = 1.187;  % S1; 
y0(2) = 0.193;  % S2;
y0(3) = 0.050;  % S3;
y0(4) = 0.115;  % S4;
y0(5) = 0.077;  % N2;
y0(6) = 2.475;  % A3;
y0(7) = 0.077;  % S4ex;

%% Running the original model %%
% BCD 8.19.2011 changing to make a function of temperature
%[ T, Y_orig ] = ode45( @Ruoff_model_original, t, y0 );
temperature = 293.; %286.5;
[ T, Y_orig ] = Ruoff_model_original(t,y0,temperature);
% end BCD 8.19.2011

% Plot the results
figure(1);
subplot(3,3,1);
plot(T,Y_orig(:,5),'g');
xlabel('Time (min)'); 
ylabel('NADH');
subplot(3,3,2);
plot(T,Y_orig(:,1),'g');
xlabel('Time (min)'); 
ylabel('Glucose');
subplot(3,3,3);
plot(T,Y_orig(:,2),'g');
xlabel('Time (min)'); 
ylabel('GAP/DAP');
subplot(3,3,4);
plot(T,Y_orig(:,3),'g');
xlabel('Time (min)'); 
ylabel('BPG');
subplot(3,3,5);
plot(T,Y_orig(:,4),'g');
xlabel('Time (min)'); 
ylabel('Pyruvate/Acetaldehyde');
subplot(3,3,6);
plot(T,Y_orig(:,6),'g');
xlabel('Time (min)'); 
ylabel('ATP');
subplot(3,3,7);
plot(T,Y_orig(:,7),'g');
xlabel('Time (min)'); 
ylabel('Ext. Pyruvate/Acetaldehyde');
% BCD 8.19.2011
% stuff from Ruoff_model_original.m
x_temp=[0.719,
14.029,
1.341,
4.906,
17.809,
0.6635,
5.361,
0.5214,
6.079,
0.3258];
J0    = x_temp(1);
k1    = x_temp(2); 
k2    = x_temp(3);
k3    = x_temp(4);
k4    = x_temp(5);  
k5    = x_temp(6); 
k6    = x_temp(7);
k     = x_temp(8); 
kappa = x_temp(9); 
q     = 4.0; 
K1    = x_temp(10);
psi   = 0.1; 
subplot(3,3,8);
S1 = Y_orig(:,1);
A3 = Y_orig(:,6);
v1 = k1*S1.*A3./(1+(A3/K1).^q);
plot(T,v1,'g');
xlabel('Time (min');
ylabel('v1');

figure(5);
plot(T,Y_orig(:,5),'g');
xlabel('Time (min)'); 
ylabel('NADH');

% end BCD

figure(2);
index_1 = 1;
index_2 = 5;
hold off;
plot(Y_orig(:,index_1),Y_orig(:,index_2), 'g' )
xlabel('Glucose'); 
ylabel('NADH');
hold off;

figure(3);
index_1 = 5;
index_2 = 6;
hold off;
plot(Y_orig(:,index_1),Y_orig(:,index_2), 'g' )
xlabel('NADH'); 
ylabel('ATP');
hold off;

figure(4);
index_1 = 1;
index_2 = 6;
hold off;
plot(Y_orig(:,index_1),Y_orig(:,index_2), 'g' )
xlabel('Glucose'); 
ylabel('ATP');
hold off;



