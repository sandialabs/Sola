fprintf(1,'\nExecuting tests for automatic_differentiation:\n');

cd Thermal/;
Test_1;
cd ..;

cd Transient_Thermal/;
Test_1;
cd ..;
