clear
close all
clc

lofi_mms_test = false;
hifi_mms_test = true;
test_fun_to_z_mapping = false;

if lofi_mms_test
    h = 2.^(-1:-1:-5);
    N = length(h);
    error = zeros(N,1);
    for k = 1:N
        mesh = PDE_Meshing(h(k));
        x = mesh.x;
        y = mesh.y;
        M = mesh.M;
        
        diff_react_lofi = Diff_React_Lofi(mesh);
        
        u_true = cos(pi*x).*cos(pi*y);
        z = (1+2*pi^2)*cos(pi*x).*cos(pi*y);
        u = diff_react_lofi.State_Solve(z);
        
        diff = u_true-u;
        error(k) = sqrt(diff'*M*diff);
    end
    
    for k = 1:N
        disp(['h = ',num2str(h(k)),' and log10(error) = ',num2str(error(k))])
    end
end

if hifi_mms_test
    h = 2.^(-1:-1:-4);
    N = length(h);
    error = zeros(N,1);
    for k = 1:N
        mesh = PDE_Meshing(h(k));
        x = mesh.x;
        y = mesh.y;
        M = mesh.M;
        
        diff_react_hifi = Diff_React_Hifi(mesh);
        
        control_fun = @(x,y) (2*pi^2)*cos(pi*x).*cos(pi*y) + (1+.2*sin(pi*x).^2).*(cos(pi*x).*cos(pi*y)).^3;
        u_true = cos(pi*x).*cos(pi*y);
        u = diff_react_hifi.State_Solve(control_fun);
        
        diff = u_true-u;
        error(k) = sqrt(diff'*M*diff);
    end
    
    for k = 1:N
        disp(['h = ',num2str(h(k)),' and log10(error) = ',num2str(error(k))])
    end
end

if test_fun_to_z_mapping
    mesh = PDE_Meshing(.1);
    x = mesh.x;
    y = mesh.y;
    M = mesh.M;
    diff_react_hifi = Diff_React_Hifi(mesh);
    
    z = x + y;
    control_fun = diff_react_hifi.Map_z_to_Control_Fun(z);
    z_tmp = diff_react_hifi.Map_Control_Fun_to_z(control_fun);
    if norm(z-z_tmp)~=0
       disp('Error in mapping z to Control_Fun') 
    end
end
