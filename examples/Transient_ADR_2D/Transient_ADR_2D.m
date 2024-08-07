classdef Transient_ADR_2D < handle
    % Solve the following two-species transient advection-diffusion-reaction equations
    %
    %     du1/dt - \kappa_1\Delta u_1 = -\alpha_1 v \cdot u - \rho u_1 u_2
    %     du2/dt - \kappa_2\Delta u_2 = -\alpha_2 v \cdot u + f
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
        diffusion       % Diffusion coefficients (default = 0.05).
        advection       % Advection coefficients (default = 4.00).
        reactant        % Reaction coefficient (default = 0.01).
        init_center     % Center of the initial condition blob (2 x 1).
        v_weights       % Weights forcing the velocity to obey no-slip conditions.
    end

    properties (Dependent)
        x               % :math:`x`-coordinates of the spatial nodes.
        y               % :math:`y`-coordinates of the spatial nodes.
        n_x             % Number of spatial nodes.
        n_y             % Dimension :math:`n_y = 2n_x` of the ODE state.
        n_q             % Number of control nodes :math:`n_q`.
    end

    methods

        function [x] = get.x(this)
            x = this.model.Mesh.Nodes(1, :)';
        end

        function [y] = get.y(this)
            y = this.model.Mesh.Nodes(2, :)';
        end

        function [nx] = get.n_x(this)
            nx = length(this.x);
        end

        function [ny] = get.n_y(this)
            ny = 2 * this.n_x;
        end

        function [n_q] = get.n_q(this)
            n_q = size(this.control_nodes, 2);
        end

    end

    methods (Access = public)

        %% Constructor
        function this = Transient_ADR_2D(model, init_center, diffusion_coeffs, advection_coeffs, reaction_coeff, nodes)
            % Initialize (but do not solve) the problem.
            %
            % Parameters
            % ----------
            % model
            %   Initialized pdetoolbox model object (made with ``model = createpde(2);``).
            % init_center : [float, float]
            %   Center of the initial condition blob.
            % diffusion_coeffs : [float, float]
            %   Diffusion coefficients. Larger means more diffusion.
            % advection_coeffs : [float, float]
            %   Advection coefficients. Larger means more advection.
            % reaction_coeff : float
            %   Reaction coefficient. Larger means more reaction.
            % nodes : int or str
            %   Number of sources OR a 2xn_q matrix of center coordinates.

            arguments
                model
                init_center
                diffusion_coeffs {mustBePositive} = [0.05, 0.05]
                advection_coeffs {mustBePositive} = [4.00, 4.00]
                reaction_coeff {mustBePositive} = 10
                nodes = 16
            end

            % Initialize the control node centers.
            if isscalar(nodes)
                qx = linspace(-1, 1, floor(sqrt(nodes)) + 2);
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
            this.diffusion = reshape(diffusion_coeffs, 2, 1);
            this.advection = reshape(advection_coeffs, 2, 1);
            this.reactant = reaction_coeff;
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

            objective = Transient_ADR_2D_Objective(center, this.x, this.y, M, T, n_t, this.n_q, beta_reg);
        end

        %% Problem definition
        function [out] = Initial_Contaminant(this, x, y)
            % Initial condition for the contaminant.
            %
            % Parameters
            % ----------
            % x
            %   x-coordinates at which to evaluate the initial condition.
            out = 50 * exp(-50 .* sum(([x; y] - this.init_center).^2, 1)) + 1;
        end

        function [v] = Velocity(this, x, y)
            % Constant velocity field: constant -> in x, sin(x) in y,
            % with x-dependent rotation.
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

        function [r] = Reaction(~, y1, y2)
            % Reaction term: :math:`-y_1(t) y_2(t)`.
            r = -y1 .* y2;
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
        function [u] = State_Solve(this, q, t)
            % Solve the PDE over the specified times.
            arguments
                this
                q
                t (:, :) {mustBeNumeric}
            end

            % (re-)Define the PDE coefficients.
            fCoef = @(lc, st) this.AdvectionTerm(lc, st) + this.ReactionTerm(lc, st) + this.SourceTerm(q, lc, st);
            specifyCoefficients(this.model, "m", 0, "d", 1, ...
                                "c", this.diffusion, "a", 0, "f", fCoef);

            % Set the initial conditions.
            setInitialConditions(this.model, @this.Initial_Condition);

            u = solvepde(this.model, t);
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

        function Plot_Field(this, y, name, logscale, animationfig)
            % Plot an (n_x x 2) field in two subplots.
            arguments
                this
                y
                name {mustBeText} = ""
                logscale = false
                animationfig = false
            end

            if animationfig
                figure(50);
            else
                fig = figure();
                fig.Position(3:4) = [830, 300];
            end

            if size(y, 2) == 1
                y_new = zeros(this.n_x, 2);
                y_new(:, 1) = y(1:this.n_x, 1);
                y_new(:, 2) = y(this.n_x + 1:end, 1);
                y = y_new;
            end

            for i = 1:2
                subplot(1, 2, i);
                pdeplot(this.model.Mesh, XYData = y(:, i), ColorMap = "parula");
                if logscale
                    set(gca, 'ColorScale', 'log');
                end
                colorbar;
                title(name);
            end
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
            fig = figure(50);
            fig.Position(3:4) = [830, 300];

            if ndims(u) == 2
                unew = zeros(this.n_x, 2, size(u, 2));
                unew(:, 1, :) = u(1:this.n_x, :);
                unew(:, 2, :) = u(this.n_x + 1:end, :);
                u = unew;
            end

            n_t = size(u, 3);
            waittime = 2 / n_t;
            umax = max(abs(u), [], "all");
            umin = min(abs(u), [], "all");
            limits = [umin, umax];

            for j = 1:n_t
                ys = [u(:, 1, j), u(:, 2, j)];
                this.Plot_Field(abs(ys), ['t = t_{', num2str(j), '}'], true, true);
                subplot(1, 2, 1);
                clim(limits);
                subplot(1, 2, 2);
                clim(limits);
                pause(waittime);
            end
        end

    end

    %% Interface methods for PDE toolbox
    methods (Access = protected)

        function [u0] = Initial_Condition(this, loc)
            % Initial condition: a Gaussian blob for y1 and zero for y2.
            M = length(loc.x);
            y1_0 = this.Initial_Contaminant(loc.x, loc.y);
            u0 = [reshape(y1_0, 1, M); zeros(1, M)];
        end

        function [out] = AdvectionTerm(this, loc, state)
            % Advection term, [v \cdot \nabla y_1; v \cdot \nabla y_2].
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
            % out (2 x n)
            %   Advection term evaluated at the spatial locations.
            velocity = this.Velocity(loc.x, loc.y)';
            agrad1 = this.advection(1) .* [state.ux(1, :); state.uy(1, :)];
            agrad2 = this.advection(2) .* [state.ux(2, :); state.uy(2, :)];
            out = [-sum(velocity .* agrad1, 1); -sum(velocity .* agrad2, 1)];
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
            % out (2 x n)
            %   Source term evaluated at the spatial locations.
            n = length(loc.x);
            qt = q(state.time);
            f = this.Source(qt, loc.x, loc.y)';
            out = [zeros(1, n); f];
        end

        function [out] = ReactionTerm(this, loc, state)
            % Reaction term.
            %
            % Parameters
            % ----------
            % loc
            %   PDE Toolbox object with ``x`` and ``y`` properties, row
            %   vectors representing ``n`` spatial locations.
            % state
            %   PDE Toolbox object with a ``u`` property, a row vector
            %   representing the state at the ``n`` spatial locations
            %   described by ``loc``.
            %
            % Returns
            % -------
            % out (2 x n)
            %   Reaction term evaluated at the spatial locations.
            y1y2 = this.Reaction(state.u(1, :), state.u(2, :));
            % n = length(loc.x);
            % out = [this.reactant * y1y2; zeros(1, n)];
            out = this.reactant .* [y1y2; y1y2];
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
            model = createpde(2);

            % Set up the geometry.
            load(loadfile, 'points', 'triangles');
            geometryFromMesh(model, points, triangles(1:3, :));

            % Enforce homogeneous Dirichlet boundary conditions everywhere.
            nE = model.Geometry.NumEdges;
            % applyBoundaryCondition(model, "dirichlet", "Edge", 1:nE, u = [0; 0]);
            applyBoundaryCondition(model, "neumann", "Edge", 1:nE, g = [0; 0], q = [0, 0; 0, 0]);
            % applyBoundaryCondition(model, "dirichlet", "Edge", [1, 44, 47, 74], u = [0; 0]);
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
            model = createpde(2);

            % Set up the (square) geometry.
            geometryFromEdges(model, @squareg);
            generateMesh(model, "Hmax", Hmax);

            % Enforce homogeneous Dirichlet boundary conditions.
            applyBoundaryCondition(model, "neumann", "Edge", 1:4, g = [0; 0], q = [0, 0; 0, 0]);
        end

    end
end
