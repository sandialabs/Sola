%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

% The automatic_differentiation module depends on adigator.
% Set "test_automatic_differentiation=true" to enable the test

% The python_adapter module depends on specific versions of python
% and Matlab. Set "test_python_adapter=true" to enable the test

test_automatic_differentiation = false;
test_bayesian_inversion = true;
test_linear_algebra_tools = true;
test_model_discrepancy = true;
test_model_reduction = true;
test_optimal_experimental_design = true;
test_optimization = true;
test_optimization_under_uncertainty = true;
test_pseudo_time_continuation = true;
test_python_adapter = false;

if ~test_automatic_differentiation
    fprintf([ ...
             '     \n', ...
             '     I am not testing the source code under sola/src/automatic_differentiation. \n', ...
             '     If ADiGator has been installed, modify the boolean variable test_automatic_differentiation to enable the test. \n', ...
             '     \n']);
end

if ~test_python_adapter
    fprintf([ ...
             '     I am not testing the source code under sola/src/python_adapter. \n', ...
             '     If the necessary Matlab and Python versions are installed, modify the boolean variable test_python_adapter to enable the test. \n', ...
             '     \n']);
end

save('Test_Settings.mat');

if test_automatic_differentiation
    cd automatic_differentiation/;
    run Execute_Tests.m;
    cd ../;
end

clear;
load('Test_Settings.mat');
if test_bayesian_inversion
    cd bayesian_inversion/;
    run Execute_Tests.m;
    cd ../;
end

clear;
load('Test_Settings.mat');
if test_linear_algebra_tools
    cd linear_algebra_tools/;
    run Execute_Tests.m;
    cd ../;
end

clear;
load('Test_Settings.mat');
if test_model_discrepancy
    cd model_discrepancy/;
    run Execute_Tests.m;
    cd ../;
end

clear;
load('Test_Settings.mat');
if test_model_reduction
    cd model_reduction/;
    run Execute_Tests.m;
    cd ../;
end

clear;
load('Test_Settings.mat');
if test_optimal_experimental_design
    cd optimal_experimental_design/;
    run Execute_Tests.m;
    cd ../;
end

clear;
load('Test_Settings.mat');
if test_optimization
    cd optimization/;
    run Execute_Tests.m;
    cd ../;
end

clear;
load('Test_Settings.mat');
if test_optimization_under_uncertainty
    cd optimization_under_uncertainty/;
    run Execute_Tests.m;
    cd ../;
end

clear;
load('Test_Settings.mat');
if test_pseudo_time_continuation
    cd pseudo_time_continuation/;
    run Execute_Tests.m;
    cd ../;
end

clear;
load('Test_Settings.mat');
if test_python_adapter
    cd python_adapter/;
    run Execute_Tests.m;
    cd ../;
end

clear;
delete('Test_Settings.mat');
