% RedBlackTree
% Implementation of a red black tree. A red black tree is a self-balancing
% binary search tree. All operations are guaranteed to be O(log n), where n
% is the size of the tree. This implementation is based on [1, Chapter 13].
%
% For convenience, we denote black by 0 and red by 1.
%
% Source:
%   [1] T. H. Cormen, Ed., Introduction to algorithms, 3rd ed. Cambridge,
%       Mass: MIT Press, 2009.

classdef RedBlackTree < handle


    properties
        root
        size
        nil
    end

    methods (Access = public)
        function this = RedBlackTree()
            this.size = 0;
            nil = RedBlackNode(0, 0);
            nil.color = 0;
            nil.left = nil;
            nil.right = nil;
            this.nil = nil;
            this.root = this.nil;
        end

        function this = Insert(this, key, data)
            z = RedBlackNode(key, data);
            y = this.nil;
            x = this.root;
            while x ~= this.nil
                y = x;
                if key < x.key
                    x = x.left;
                else
                    x = x.right;
                end
            end
            z.parent = y;
            if y == this.nil
                this.root = z;
            elseif key < y.key
                y.left = z;
            else
                y.right = z;
            end
            z.left = this.nil;
            z.right = this.nil;
            z.color = 1;
            this = this.Insert_Fixup(z);
            this.size = this.size + 1;
        end

        function u = Min(this)
            u = this.root;
            while u.left ~= this.nil
                u = u.left;
            end
        end

        function u = SubtreeMin(this, z)
            u = z;
            while u.left ~= this.nil
                u = u.left;
            end
        end

        function u = Max(this)
            u = this.root;
            while u.right ~= this.nil
                u = u.right;
            end
        end

        function u = Find(this, key)
            u = this.root;
            while u ~= this.nil
                if u.key == key
                    break
                elseif key < u.key
                    u = u.left;
                else
                    u = u.right;
                end
            end
        end

        function this = Delete(this, key)
            z = this.Find(key);
            if z ~= this.nil
                this.Delete_Node(z);
            end
        end

        function [key, data] = PopMax(this)
            u = this.Max();
            key = u.key;
            data = u.data;
            this.Delete_Node(u);
        end

        function [key, data] = PopMin(this)
            u = this.Min();
            key = u.key;
            data = u.data;
            this.Delete_Node(u);
        end

    end

    methods (Access = private)
        function this = Delete_Node(this, z)
            y = z;
            y_original_color = y.color;
            if z.left == this.nil
                x = z.right;
                this.Transplant(z, z.right);
            elseif z.right == this.nil  
                x = z.left;
                this.Transplant(z, z.left);
            else
                y = this.SubtreeMin(z.right);
                y_original_color = y.color;
                x = y.right;
                if y.parent == z
                    x.parent = y;
                else
                    this.Transplant(y, y.right);
                    y.right = z.right;
                    y.right.parent = y;
                end
                this.Transplant(z, y);
                y.left = z.left;
                y.left.parent = y;
                y.color = z.color;
            end
            if y_original_color == 0
                this.Delete_Fixup(x);
            end
            this.size = this.size - 1;
        end

        function this = Left_Rotate(this, x)
            y = x.right;
            x.right = y.left;
            if y.left ~= this.nil
                y.left.parent = x;
            end
            y.parent = x.parent;
            if x.parent == this.nil
                this.root = y;
            elseif x == x.parent.left
                x.parent.left = y;
            else
                x.parent.right = y;
            end
            y.left = x;
            x.parent = y;
        end

        function this = Right_Rotate(this, y)
            x = y.left;
            y.left = x.right;
            if x.right ~= this.nil
                x.right.parent = y;
            end
            x.parent = y.parent;
            if y.parent == this.nil
                this.root = x;
            elseif y == y.parent.left
                y.parent.left = x;
            else
                y.parent.right = x;
            end
            x.right = y;
            y.parent = x;
        end

        function this = Insert_Fixup(this, z)
            while z.parent.color == 1
                if z.parent == z.parent.parent.left
                    y = z.parent.parent.right;
                    if y.color == 1
                        z.parent.color = 0;
                        y.color = 0;
                        z.parent.parent.color = 1;
                        z = z.parent.parent;
                    else
                        if z == z.parent.right
                            z = z.parent;
                            this = this.Left_Rotate(z);
                        end
                        z.parent.color = 0;
                        z.parent.parent.color = 1;
                        this = this.Right_Rotate(z.parent.parent);
                    end
                else
                    y = z.parent.parent.left;
                    if y.color == 1
                        z.parent.color = 0;
                        y.color = 0;
                        z.parent.parent.color = 1;
                        z = z.parent.parent;
                    else
                        if z == z.parent.left
                            z = z.parent;
                            this = this.Right_Rotate(z);
                        end
                        z.parent.color = 0;
                        z.parent.parent.color = 1;
                        this = this.Left_Rotate(z.parent.parent);
                    end
                end
            end
            this.root.color = 0;
        end

        function this = Transplant(this, u, v)
            if u.parent == this.nil
                this.root = v;
            elseif u == u.parent.left
                u.parent.left = v;
            else
                u.parent.right = v;
            end
            v.parent = u.parent;
        end

        function this = Delete_Fixup(this, x)
            while x ~= this.root && x.color ~= 0
                if x == x.parent.left
                    w = x.parent.right;
                    if w.color == 1
                        w.color = 0;
                        x.parent.color = 1;
                        this.Left_Rotate(x.parent);
                        w = x.parent.right;
                    end
                    if w.left.color == 0 && w.right.color == 0
                        w.color = 1;
                        x = x.parent;
                    else
                        if w.right.color == 0
                            w.left.color = 0;
                            w.color = 1;
                            this.Right_Rotate(w);
                            w = x.parent.right;
                        end
                        w.color = x.parent.color;
                        x.parent.color = 0;
                        w.right.color = 0;
                        this.Left_Rotate(x.parent);
                        x = this.root;
                    end
                else
                    w = x.parent.left;
                    if w.color == 1
                        w.color = 0;
                        x.parent.color = 1;
                        this.Right_Rotate(x.parent);
                        w = x.parent.left;
                    end
                    if w.right.color == 0 && w.left.color == 0
                        w.color = 1;
                        x = x.parent;
                    else
                        if w.left.color == 0
                            w.right.color = 0;
                            w.color = 1;
                            this.Left_Rotate(w);
                            w = x.parent.left;
                        end
                        w.color = x.parent.color;
                        x.parent.color = 0;
                        w.left.color = 0;
                        this.Right_Rotate(x.parent);
                        x = this.root;
                    end
                end
            end
            x.color = 0;
        end
    end
end
