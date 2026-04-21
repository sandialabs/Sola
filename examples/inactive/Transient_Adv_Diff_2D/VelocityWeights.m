%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [weights] = VelocityWeights(model, eps)
    % Create a weight matrix to ensure that a velocity field has proper
    % behavior at the boundaries.
    %
    % Parameters
    % ----------
    % model
    %   pde toolkit object, created with createpde().
    %
    % Returns
    % -------
    % weights : (num_nodes, 1) vector
    %     Weights to multiply the vector field. Zero at the boundary, etc.

    arguments
        model
        eps {mustBePositive} = 1000
    end

    boundaryNodes = findNodes(model.Mesh, 'region', 'edge', 1:model.Geometry.NumEdges);
    coordinates = model.Mesh.Nodes(:, boundaryNodes);

    % Sanity check: plot the boundary nodes.
    % figure;
    % scatter(coordinates(1, :), coordinates(2, :));

    % Get the distance from each mesh node to the mesh boundary.
    num_elements = size(model.Mesh.Nodes, 2);
    distances = zeros(num_elements, 1);
    for i = 1:num_elements
        node = model.Mesh.Nodes(:, i);
        nearest_index = dsearchn(coordinates', model.Mesh.Nodes(:, i)');
        nearest_Bnode = coordinates(:, nearest_index);
        distances(i) = sum((nearest_Bnode - node).^2);
    end

    weights = 1 - exp(-eps .* distances);

end
