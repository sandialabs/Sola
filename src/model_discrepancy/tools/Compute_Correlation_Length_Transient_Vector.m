function [corr_length] = Compute_Correlation_Length_Transient_Vector(t,z,initial_guess)

    n_z = length(z);
    n_t = length(t);
    n_q = n_z/n_t;
    z = reshape(z,n_q,n_t)';
    corr_lengths = zeros(n_q,1);
    for k = 1:n_q
        corr_lengths(k) = Compute_Correlation_Length_1D(t,z(:,k),initial_guess);
    end
    corr_length = mean(corr_lengths);

end