classdef Transient_Adv_Diff_2D < handle
    % Solve the transient advection-diffusion equation
    %
    %     du/dt - \Delta u - v \cdot u = f
    %
    % over a two-dimensional domain with
    % homogeneous Dirichlet boundary conditions.
    % The source term is a sum of Gaussian bubble functions,
    %
    %     f(x,y,t) = \sum_{i=1}^{m} q_i(t)\phi_i(x,y),
    %
    % where
    %
    %     \phi_i(x,y) = exp(-200((x - x_i)^2 + (y - y_i)^2))
    %
    % is a Gaussian function centered at (x_i,y_i).
    % The PDE model and source centers are specified in the constructor.

    properties
        model           % PDE toolbox model, made with createpde.
        control_nodes   % Coordinates of the control nodes (2 x n_q).
        diffusion       % Diffusion coefficient (default = 0.05).
        advection       % Advection coefficient (default = 4.00).
        init_center     % Center of the initial condition blob (2 x 1).
        v_weights       % Weights forcing the velocity to obey no-slip conditions.
    end

    properties (Dependent)
        x               % :math:`x`-coordinates of the spatial nodes.
        y               % :math:`y`-coordinates of the spatial nodes.
        n_q             % Number of control nodes :math:`n_q`.
    end

    methods

        function [x] = get.x(this)
            x = this.model.Mesh.Nodes(1, :)';
        end

        function [y] = get.y(this)
            y = this.model.Mesh.Nodes(2, :)';
        end

        function [n_q] = get.n_q(this)
            n_q = size(this.control_nodes, 2);
        end

    end

    methods (Access = public)

        %% Constructor
        function this = Transient_Adv_Diff_2D(model, init_center, diffusion_coeff, advection_coeff, nodes)
            % Initialize (but do not solve) the problem.
            %
            % Parameters
            % ----------
            % diffusion_coeff : float
            %   Diffusion coefficient. Larger means more diffusion.
            % advection_coeff : float
            %   Advection coefficient. Larger means more advection.
            % nodes : int or str
            %   Number of sources OR a 2xn_q matrix of center coordinates.

            arguments
                model
                init_center
                diffusion_coeff {mustBePositive} = 0.05
                advection_coeff {mustBePositive} = 4.00
                nodes = 25
            end

            % Initialize the control node centers.
            if isscalar(nodes)
                qx = linspace(-1, 1, ceil(sqrt(nodes)) + 2);
                qx = qx(2:end - 1);
                control_nodes = combinations(qx, qx);
                this.control_nodes = control_nodes{:, :}';
            else
                this.control_nodes = nodes;
            end

            % Calculate velocity weights.
            boundaryNodes = findNodes(model.Mesh, 'region', 'edge', 1:model.Geometry.NumEdges);
            coordinates = model.Mesh.Nodes(:, boundaryNodes);
            num_elements = size(model.Mesh.Nodes, 2);
            distances = zeros(num_elements, 1);
            for i = 1:num_elements
                node = model.Mesh.Nodes(:, i);
                nearest_index = dsearchn(coordinates', node');
                nearest_Bnode = coordinates(:, nearest_index);
                distances(i) = sum((nearest_Bnode - node).^2);
            end
            this.v_weights = 1 - exp(-1000 .* distances);

            % Save the model.
            this.diffusion = diffusion_coeff;
            this.advection = advection_coeff;
            this.init_center = reshape(init_center, 2, 1);
            this.model = model;
        end

        function [objective] = Make_Objective(this, center, T, n_t, beta_reg)
            % Construct a Transient_Adv_Diff_2D_Objective corresponding to
            % the model mesh.
            %
            % Parameters
            % ----------
            % center : (2 x 1)
            %   Spatial location where the objective most heavily penalizes
            %   the presence of contaminant.
            % T : float
            %   Final simulation time.
            % n_t : int
            %   Number of time steps.
            % beta_reg : float
            %   Regularization constant for the control magnitudes.
            arguments
                this
                center
                T
                n_t {mustBePositive, mustBeInteger}
                beta_reg = 1e-3
            end

            % Extract the mass matrix from the model.
            M = assembleFEMatrices(this.model, 'M').M;

            % % Assemble the control matrix.
            % xy = [this.x, this.y];
            % Bq = zeros(this.n_y, this.n_q);
            % for j = 1:this.n_q
            %     center = this.control_nodes(:, j)';
            %     Bq(:, n_q) = 50 * exp(-1000 .* sum((xy - this.control_nodes).^2, 1));
            % end

            objective = Transient_Adv_Diff_2D_Objective(center, this.x, this.y, M, T, n_t, this.n_q, beta_reg);
        end

        %% Problem definition
        function [u0] = Initial_Condition(this, loc)
            % Initial condition: a Gaussian blob.
            u0 = 20 * exp(-100 .* sum(([loc.x; loc.y] - this.init_center).^2, 1));
        end

        function [v] = Velocity(this, x, y)
            % Constant velocity field: constant -> in x, sin(x) in y.
            %
            % Parameters
            % ----------
            % x : (n x 1)
            %   x coordinates at which to evaluate the velocity field.
            % y : (n x 1)
            %   y coordinates at which to evaluate the velocity field.
            %
            % Returns
            % -------
            % v : (n x 2)
            %   Velocity in x and y directions at the given coordinates.
            xx = reshape(x, [], 1);
            xmin = min(xx);
            xmax = max(xx);
            xspan = xmax - xmin;
            n = size(xx, 1);
            flow = [ones(n, 1), sin(10 * pi .* xx) / 2];
            for i = 1:n
                dt = (xx(i) - xmin) / xspan;
                angl = pi * (-75 * dt) / 180;
                cosangl = cos(angl);
                sinangl = sin(angl);
                rotation = [cosangl -sinangl; sinangl cosangl];
                flow(i, :) = flow(i, :) * rotation';
                % % Weight velocity field for no-slip condition. EXPENSIVE!
                % idx = dsearchn(this.model.Mesh.Nodes', [x(i), y(i)]);
                % flow(i, :) = flow(i, :) * this.v_weights(idx);
            end
            v = flow;
        end

        function [f] = Source(this, q, x, y)
            % Gaussian source function.
            %
            %     f(x,y,t) = \sum_{i=1}^{n_q} q_i(t)\phi_i(x,y),
            %
            % where
            %
            %     \phi_i(x,y) = exp(-200((x - x_i)^2 + (y - y_i)^2))
            %
            % Parameters
            % ----------
            % q : (n_q x 1)
            %   Coefficients for each Gaussian bubble at the current time,
            %   i.e., q = [q_1(t) ... q_{n_q}(t)]'.
            % x : (n x 1)
            %   x coordinates at which to evaluate the source.
            % y : (n x 1)
            %   y coordinates at which to evaluate the source.
            %   (not used because velocity depends on x only).
            %
            % Returns
            % -------
            % f : (n x 1)
            %   Value of the source function at the given coordinates.
            x = reshape(x, [], 1);
            n = size(x, 1);
            f = zeros(n, 1);
            qflat = reshape(q, 1, []);
            for i = 1:size(x, 1)
                xy = [x(i); y(i)];
                phis = 50 * exp(-1000 .* sum((xy - this.control_nodes).^2, 1));
                f(i) = qflat * phis';
            end
        end

        %% Solver
        function [u] = State_Solve(this, q, t, animate)
            % Solve the PDE over the specified times.
            arguments
                this
                q
                t (:, :) {mustBeNumeric}
                animate = false
            end

            % (re-)Define the PDE coefficients.
            fCoef = @(lc, st) this.AdvectionTerm(lc, st) + this.SourceTerm(q, lc, st);
            specifyCoefficients(this.model, "m", 0, "d", 1, ...
                                "c", this.diffusion, "a", 0, "f", fCoef);

            % Set the initial conditions.
            setInitialConditions(this.model, @this.Initial_Condition);

            u = solvepde(this.model, t);

            if animate
                this.Animate_Solution(u.NodalSolution);
            end
        end

        %% Visualization
        function Plot_Mesh(this)
            % Visualize just the two-dimensional finite element mesh.

            % figure;
            % pdegplot(this.model, "EdgeLabels", "on");
            % axis equal;
            % xlim([-1.1 1.1]);
            % ylim([-1.1 1.1]);
            % title("Geometry with Edge Labels");

            figure;
            pdemesh(this.model);
            % axis equal;
            xlim([min(this.x, [], "all") max(this.x, [], "all")]);
            ylim([min(this.y, [], "all") max(this.y, [], "all")]);
            title("Finite Element Mesh");
        end

        function Plot_Control_Nodes(this)
            % Plot the control node locations on top of the mesh.
            this.Plot_Mesh();
            hold on;
            scatter(this.control_nodes(1, :), this.control_nodes(2, :), 72, "black", "filled", "o");
        end

        function Plot_Field(this, u, name, logscale, animationfig)
            arguments
                this
                u
                name {mustBeText} = ""
                logscale = false
                animationfig = false
            end

            if animationfig
                figure(50);
            else
                figure;
            end

            pdeplot(this.model.Mesh, XYData = u, ColorMap = "sky");
            colormap(redblue);
            if logscale
                set(gca, 'ColorScale', 'log');
            end
            colorbar;
            title(name);
        end

        function [X, Y, VX, VY] = Plot_Vector_Field(this, v, resolution, name)
            % Plot a vector field with a quiver plot over a 2D domain.
            arguments
                this
                v
                resolution {mustBePositive} = 50
                name {mustBeText} = ""
            end

            xmin = min(this.x);
            xmax = max(this.x);
            ymin = min(this.y);
            ymax = max(this.y);
            xx = linspace(xmin, xmax, resolution)';
            yy = linspace(ymin, ymax, resolution)';
            [X, Y] = meshgrid(xx, yy);

            vx = scatteredInterpolant(this.x, this.y, v(:, 1));
            vy = scatteredInterpolant(this.x, this.y, v(:, 2));
            VX = vx(X, Y);
            VY = vy(X, Y);

            quiver(X, Y, VX, VY);
            xlim([xmin, xmax]);
            ylim([ymin, ymax]);
            title(name);
        end

        function Plot_Velocity_Field(this, resolution)
            arguments
                this
                resolution {mustBePositive} = 50
            end

            figure;
            xmin = min(this.x);
            xmax = max(this.x);
            ymin = min(this.y);
            ymax = max(this.y);
            xx = linspace(xmin, xmax, resolution)';
            yy = linspace(ymin, ymax, resolution)';
            [X, Y] = meshgrid(xx, yy);
            VXVY = this.Velocity(reshape(X, [], 1), reshape(Y, [], 1));
            VX = reshape(VXVY(:, 1), size(X));
            VY = reshape(VXVY(:, 2), size(Y));

            quiver(X, Y, VX, VY);
            xlim([xmin, xmax]);
            ylim([ymin, ymax]);
            title('Velocity field (wind)');
        end

        function Animate_Solution(this, u)
            n_t = size(u, 2);
            waittime = 2 / n_t;
            umax = max(abs(u), [], "all") * 2 / 3;
            limits = [-umax * 2 / 3, umax];

            for j = 1:n_t
                this.Plot_Field(u(:, j), ['t = t_{', num2str(j), '}'], false, true);
                clim(limits);
                pause(waittime);
            end
        end

    end

    %% Helper methods
    methods (Access = protected)

        function [out] = AdvectionTerm(this, loc, state)
            % Advection term v \cdot \nabla u.
            %
            % Parameters
            % ----------
            % loc
            %   PDE Toolbox object with ``x`` and ``y`` properties, row
            %   vectors representing ``n`` spatial locations.
            % state
            %   PDE Toolbox object with ``ux`` and ``uy`` properties, row
            %   vectors representing the gradient of the state at the ``n``
            %   spatial locations described by ``loc``.
            %
            % Returns
            % -------
            % out (1 x n)
            %   Advection term evaluated at the spatial locations.
            velocity = this.Velocity(loc.x, loc.y)';
            grad = [state.ux; state.uy];
            out = -this.advection * sum(velocity .* grad, 1);
        end

        function [out] = SourceTerm(this, q, loc, state)
            % Source term f.
            %
            % Parameters
            % ----------
            % q
            %   Function handle mapping time ``state.time`` to the ``n_q``
            %   coefficients :math:`phi_i(t)` defining the input function.
            % loc
            %   PDE Toolbox object with ``x`` and ``y`` properties, row
            %   vectors representing ``n`` spatial locations.
            % state
            %   PDE Toolbox object with ``ux`` and ``uy`` properties, row
            %   vectors representing the gradient of the state at the ``n``
            %   spatial locations described by ``loc``.
            %
            % Returns
            % -------
            % out (1 x n)
            %   Advection term evaluated at the spatial locations.
            qt = q(state.time);
            out = this.Source(qt, loc.x, loc.y)';
        end

    end

    %% Model initializers
    methods (Static)

        function [model] = model_fromfile(loadfile)
            % Initialize (but do not solve) the pde model, loading a mesh
            % from an existing file.
            %
            % Parameters
            % ----------
            % loadfile : str
            %   File to load mesh data from. Must have the following fields.
            %   - 'points'
            %   - 'triangles'

            arguments
                loadfile {mustBeText}
            end

            % Initialize PDE object.
            model = createpde;

            % Set up the geometry.
            load(loadfile, 'points', 'triangles');
            geometryFromMesh(model, points, triangles(1:3, :));

            % Enforce homogeneous Dirichlet boundary conditions everywhere.
            nE = model.Geometry.NumEdges;
            % applyBoundaryCondition(model, "dirichlet", "Edge", 1:nE, u = 0);
            applyBoundaryCondition(model, "neumann", "Edge", 1:nE, g = 0, q = 0);
            % applyBoundaryCondition(model, "dirichlet", "Edge", [1, 44, 47, 74], u = 0);
        end

        function [model] = model_default(Hmax)
            % Initialize (but do not solve) the pde model. This option
            % creates a default box geometry with a given mesh parameter.
            %
            % Parameters
            % ----------
            % Hmax : float
            %   Maximum spacing in the spatial mesh. Smaller Hmax means a
            %   finer mesh (more degrees of freedom).

            arguments
                Hmax {mustBePositive} = 0.05
            end

            % Initialize PDE object.
            model = createpde;

            % Set up the (square) geometry.
            geometryFromEdges(model, @squareg);
            generateMesh(model, "Hmax", Hmax);

            % Enforce homogeneous Dirichlet boundary conditions.
            applyBoundaryCondition(model, "dirichlet", "Edge", 1:4, "u", 0);
        end

    end
end
