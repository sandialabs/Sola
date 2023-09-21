clear
close all
clc

test_optimization = true;

% Output to screen will identify failed tests, no output means all tests passed

if test_optimization
    
    cd optimization/
    
    cd Example_1/
    Driver_Example_1
    cd ..
    
    cd Example_2/
    Driver_Example_2
    cd ..
    
    cd Example_3/
    Driver_Example_3
    cd ..
    
    cd Example_4/
    Driver_Example_4
    cd ..
    
    cd Example_5/
    Driver_Example_5
    cd ..
    
    cd ..
    
end

test_model_discrepancy = true;

if test_model_discrepancy
    
    cd model_discrepancy/
    
    cd Hessian_Randomzied_GEVP/
    Driver
    cd ..
    
    cd model_discrepancy_synthetic_test/
    Driver
    cd ..
    
    cd model_discrepancy_synthetic_test_with_hessian_gevp/
    Driver
    cd ..
    
    cd model_discrepancy_synthetic_test_with_gsvd/
    Driver
    cd ..
    
    cd PDE_Test_Problem/
    Driver_HDSA
    cd ..
    
    cd Randomized_GEVP/
    Driver
    cd ..
    
    cd ..
    
end

clear