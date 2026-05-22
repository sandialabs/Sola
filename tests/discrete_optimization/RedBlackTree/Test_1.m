clear;
close all;

N = 5;
max_sample = 100;
verbose = 1;

passed = 1;

for t = 1:3
    for s = 1:5
        tree = RedBlackTree();

        l = randperm(max_sample, N);

        for i = 1:N
            tree.Insert(l(i), 1);
        end

        height_ub = 2 * log2(N + 1);
        passed = passed & (tree.Min().key == min(l));
        passed = passed & (tree.Max().key == max(l));
        passed = passed & (Get_Height(tree) <= height_ub);
        if verbose
            fprintf('Min: %d/%d; Max: %d/%d; Height: %d/%.4f\n', ...
                    tree.Min().key, min(l), ...
                    tree.Max().key, max(l), ...
                    Get_Height(tree), height_ub);
        end
    end
    N = N * 10;
    max_sample = 10 * max_sample;
end

if passed
    fprintf(1, '\ndiscrete_optimization/RedBlackTree passed.\n');
else
    fprintf(2, '\ndiscrete_optimization/RedBlackTree failed.\n');
end

%% Functions for computing the height of the tree
function h = Height_Helper(tree, node)
    if node == tree.nil
        h = 0;
    else
        l_h = Height_Helper(tree, node.left);
        r_h = Height_Helper(tree, node.right);
        h = 1 + max(l_h, r_h);
    end
end

function h = Get_Height(tree)
    h = Height_Helper(tree, tree.root);
end
