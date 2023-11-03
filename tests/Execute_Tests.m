clear;
close all;
clc;

% Output to screen will identify failed tests, no output means all tests passed

test_optimization = true;
if test_optimization
    cd optimization/;
    Execute_Tests;
    cd ..;
end

test_model_discrepancy = true;
if test_model_discrepancy
    cd model_discrepancy/;
    Execute_Tests;
    cd ..;
end

test_bayesian_inversion = true;
if test_bayesian_inversion
    cd bayesian_inversion/;
    Execute_Tests;
    cd ..;
end

test_automatic_differentiation = true;
if test_automatic_differentiation
    cd automatic_differentiation/;
    Execute_Tests;
    cd ..;
end

clear;
