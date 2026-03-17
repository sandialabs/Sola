clear;
close all;
clc;

% The automatic_differentiation module depends on adigator.
% Set "test_automatic_differentiation=false" to disable the test
% if adigator is not available.

% The python_adapter module depends on specific versions of python
% and Matlab. Set "test_python_adapter=false" to disable the test
% if the correct versions are not available.

test_automatic_differentiation = true;
test_bayesian_inversion = true;
test_linear_algebra_tools = true;
test_model_discrepancy = true;
test_model_reduction = true;
test_optimal_experimental_design = true;
test_optimization = true;
test_optimization_under_uncertainty = true;
test_parametric_sensitivities = true;
test_python_adapter = true;

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
if test_parametric_sensitivities
    cd parametric_sensitivities/;
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
