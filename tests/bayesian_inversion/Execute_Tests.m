fprintf(1,'\nExecuting tests for bayesian_inversion:\n');

cd Example_1/;
Test_1;
cd ..;

cd Poisson/;
Test_1;
cd ..;

cd Adv_Diff/;
Test_1;
cd ..;
