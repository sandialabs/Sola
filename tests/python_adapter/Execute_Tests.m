fprintf(1, '\nExecuting tests for python_adapater:\n');

cd model_discrepancy/synthetic_test/;
Test_1;
cd ../../;

cd model_discrepancy/synthetic_test_with_gsvd;
Test_1;
cd ../../;

cd model_discrepancy/synthetic_test_with_hyperparam_auto_1D;
Test_1;
cd ../../;
