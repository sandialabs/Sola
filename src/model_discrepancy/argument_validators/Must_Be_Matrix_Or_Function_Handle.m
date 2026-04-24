function Must_Be_Matrix_Or_Function_Handle(A)
    if isa(A, 'function_handle')
        return
    end

    if (isnumeric(A) || islogical(A)) && ismatrix(A)
        return
    end

    error('Input must be either a numeric/logical matrix or a function handle.');
end
