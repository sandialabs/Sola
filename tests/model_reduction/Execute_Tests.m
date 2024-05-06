clear;
close all;
clc;

% Output to screen will identify failed tests, no output means all tests passed

test_POD = true;
if test_POD
    run Test_POD_Basis;
end

test_operators = true;
if test_operators
    run Test_Operators;
end

test_roms = true;
if test_roms
    run Test_OpInf_ROM_Constraint;
    run Test_OpInf_ROM_Constraint_ddts;
end

disp('All model_reduction tests passed!');

clear;
