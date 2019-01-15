%% run_ruoff_model_python_interface.m
%%
%% Bryan Daniels
%% 8.19.2011
%%
%% To interface Python (see simulateYeastOscillator.py) with Soni's MATLAB code.
%%
%% 

function [] = run_Ruoff_model_python_interface(pidNo,temperature,y01,y02,y03,y04,y05,y06,y07)

    timesFilename = [ 'times_temp' pidNo '.data' ];
    outputFilename = [ ...
        'run_ruoff_model_python_interface_output' ... 
        pidNo '.txt' ];
    outputParamsFilename = [ ...
        'run_ruoff_model_python_interface_output_params' ...
        pidNo '.txt' ];

    t = importdata(timesFilename)';

    %initial conditions
    y0 = zeros(7,1);
    y0(1) = y01; 
    y0(2) = y02;
    y0(3) = y03;
    y0(4) = y04;
    y0(5) = y05;
    y0(6) = y06;
    y0(7) = y07;
    
    [timesData,ydata,yderivs,params] = ... 
        Ruoff_model_original (t, y0, temperature);

    dlmwrite([outputFilename],[timesData,ydata,yderivs], ' ');
    dlmwrite([outputParamsFilename],[params], ' ');

end