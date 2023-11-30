clear;
close all;
clc;

% Output to screen will identify failed tests, no output means all tests passed

test_optimization = true;
if test_optimization
    cd optimization/;
    run Execute_Tests;
    cd ..;
end

test_model_reduction = true;
if test_model_reduction
    cd model_reduction/;
    run Execute_Tests;
    cd ../;
end

test_model_discrepancy = true;
if test_model_discrepancy
    cd model_discrepancy/;
    run Execute_Tests;
    cd ../;
end

test_bayesian_inversion = true;
if test_bayesian_inversion
    cd bayesian_inversion/;
    run Execute_Tests;
    cd ../;
end

test_automatic_differentiation = true;
if test_automatic_differentiation
    cd automatic_differentiation/;
    run Execute_Tests;
    cd ../;
end

test_linear_algebra_tools = true;
if test_linear_algebra_tools
    cd linear_algebra_tools/;
    Execute_Tests;
    cd ..;
end

clear;
