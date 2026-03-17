fprintf(1, '\nExecuting tests for model_discrepancy:\n');

cd hyperparam_2D/;
Test_1;
cd ..;

cd lumped_mass_unit_test/;
Test_1;
cd ..;

cd PDE_Test_Problem/;
Test_1;
cd ..;

cd synthetic_test/;
Test_1;
cd ..;

cd synthetic_test_control_vec/;
Test_1;
cd ..;

cd synthetic_test_elliptic_prior/;
Test_1;
cd ..;

cd synthetic_test_hessian_gevp/;
Test_1;
cd ..;

cd synthetic_test_hyperparam_1D/;
Test_1;
cd ..;

cd synthetic_test_lumped_mass/;
Test_1;
cd ..;

cd synthetic_test_multi_state/;
Test_1;
cd ..;

cd synthetic_test_transient/;
Test_1;
cd ..;

cd transient_control_synthetic_test/;
Test_1;
cd ..;

cd transient_multi_state_synthetic_test/;
Test_1;
cd ..;

cd Transient_Test_Problem/;
Test_1;
Test_2;
cd ..;
