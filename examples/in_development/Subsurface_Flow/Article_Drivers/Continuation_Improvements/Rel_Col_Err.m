%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errs = Rel_Col_Err(A, B)
    % Column-wise relative error; NaN for columns of A that are not finite.
    errs = nan(size(A, 2), 1);
    for ii = 1:size(A, 2)
        if all(isfinite(A(:, ii)))
            errs(ii) = Rel_Err(A(:, ii), B(:, ii));
        end
    end
end
