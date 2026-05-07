classdef RedBlackNode < handle
    % RedBlackNode Node in a RedBlackTree

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
