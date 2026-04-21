%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '\nExecuting tests for model_reduction:\n');
% Output to screen will identify failed tests, no output means all tests passed

test_POD = true;
if test_POD
    run Test_POD_Basis;
end

test_operators = true;
if test_operators
    run Test_Operators;
    run Test_Operators_Multi;
end

test_roms = true;
if test_roms
    run Test_OpInf_ROM_Constraint;
    run Test_OpInf_ROM_Constraint_ddts;
end

fprintf(1, '\nModel reduction tests passed.\n');

clear;
