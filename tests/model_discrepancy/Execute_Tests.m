cd analytic_laplacian_2D/;
Driver_HDSA;
cd ..;

cd hyperparam_auto_2D/;
% Missing test
cd ..;

cd model_discrepancy_oed_unit_test/;
Driver;
cd ..;

cd model_discrepancy_synthetic_test/;
Driver;
cd ..;

cd model_discrepancy_synthetic_test_with_analytic_laplacian/;
% Need to look at this a bit closer
% Driver;
cd ..;

cd model_discrepancy_synthetic_test_with_gsvd/;
Driver;
cd ..;

cd model_discrepancy_synthetic_test_with_hessian_gevp/;
Driver;
cd ..;

cd model_discrepancy_synthetic_test_with_hyperparam_auto_1D/;
Driver;
cd ..;

cd PDE_Test_Problem/;
Driver_HDSA;
cd ..;

cd Transient_Hyperparameter_auto_Test_Problem/;
% Missing test
cd ..;

cd Transient_Test_Problem/;
Driver_Unit_Test_1;
cd ..;
