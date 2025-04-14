function Must_Be_Cell_of_MD_u_Hyperparameter_Interface(inputCell)
    if ~all(cellfun(@(cell)isa(cell,'MD_u_Hyperparameter_Interface'),inputCell))
        error('All entries in the cell array must be of type MD_u_Hyperparameter_Interface.');
    end
end