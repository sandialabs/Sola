clear;
close all;

N = 10;

tree = RedBlackTree();

% Check to see if order is preserved
L = randperm(N);
for i = 1:N
    tree = tree.Insert(L(i), 1);
end

list = zeros(1, N);
for i = 1:N
    list(i) = tree.PopMax();
end
disp(list);

L = randperm(N);
for i = 1:N
    tree = tree.Insert(L(i), 1);
end

list = zeros(1, N);
for i = 1:N
    list(i) = tree.PopMin();
end
disp(list);

% Test removing elements over and over again
rng(305);
L = randperm(N);
for i = 1:N
    tree = tree.Insert(L(i), 1);
end

fprintf('Test Popping Max\n');
% 51 seems to be when it breaks;
for i = 1:50
    [k, d] = tree.PopMax();
    k = k * rand();
    tree = tree.Insert(k, d);
end

%should maybe look at the previous one... because that's where stuff fails
[k, d] = tree.PopMax();
k = k * rand();
tree.Insert(k, d);