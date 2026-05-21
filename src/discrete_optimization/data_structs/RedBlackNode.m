% RedBlackNode
% A node in a RedBlackTree [1, Chapter 13]. The main attributes of a node are
% its data and key. The key is used to determine the location of the node in the
% tree. The other fields are parent, color, left and right. However, these are
% only relevant to structure of the tree. DO NOT MODIFY THESE FIELDS.
%
% Source:
%   [1] T. H. Cormen, Ed., Introduction to algorithms, 3rd ed. Cambridge,
%       Mass: MIT Press, 2009.
%
% Author:
%   - Steven Maio (smaio@sandia.gov)

classdef RedBlackNode < handle

    properties
        color
        key
        data
        parent
        left
        right
    end

    methods (Access = public)
        function this = RedBlackNode(key, data)
            this.color = 1;
            this.key = key;
            this.data = data;
            this.parent = 0;
            this.left = 0;
            this.right = 0;
        end
    end
end
