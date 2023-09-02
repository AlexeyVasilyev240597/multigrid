% Class 'ellipticProblemSolver' implements a solver of elliptic problem,
% include iteration methods: Jacobi, Seidel, SOR and multigrid
classdef ellipticProblemSolver < handle
    properties (Access = public)
        % grid
        grid = []
        
        % finite-difference operator
        A_op = [];
        
        % right side of equation
        f_func = [];               
        % border conditions
        border = [];        
        % exact solution (if it's known)
        u_ex = [];       
        
        % spectral radius
        rho = [];        
        % error
        error = [];   

        % max number of iterations
        k_max = [];
        % accuracy
        eps = []; 
                            
        % params of MG method
        nu_1 = [];
        nu_2 = [];
        level = [];
        smooth_m = []; 
    end
    
    methods (Access = public)
        % constructor of class
        function obj = ellipticProblemSolver(c, e_l, N,...
                a_coeff, b_coeff, q_coeff,...
                f, border, u_ex)
            obj.grid = gridClass(c, e_l, N);
            obj.A_op = finDiffOpClass(a_coeff, b_coeff, q_coeff);
            obj.f_func = f;
            obj.border = border;     
            obj.u_ex = u_ex;
        end       
        
        % func returns vector with values of solution in border and
        % zeros in inner nodes
        function border_vec = fillBorder(obj, grid_m)
            N = grid_m.N;
            border_vec = zeros(N + 2);
            
            %left
            [x, y] = grid_m.leftEdge();
            border_vec(:, 1)   = obj.border(x, y);
            %right
            [x, y] = grid_m.rightEdge();
            border_vec(:, end) = obj.border(x, y);
            %bottom
            [x, y] = grid_m.bottomEdge();
            border_vec(1,   :) = obj.border(x, y);
            %top
            [x, y] = grid_m.topEdge();
            border_vec(end, :) = obj.border(x, y);               
        end
        
        % func returns vector -- right part of finite difference equation
        function g_vec = getRightPart(obj)
            N = obj.grid.N;
            g_vec = obj.fillBorder(obj.grid);
            
            % indices of inner nodes
            i_in = 2:N+1;
            j_in = 2:N+1;            
            
            g_vec(j_in, i_in) = obj.f_func(obj.grid.x(i_in),...
                                           obj.grid.y(j_in));
        end        
        
        % func returns solution of finite difference equation,
        % it is solving by multigrid method
        % PARAM IN:
        %   * nu_1, nu_2 - number of pre- and post-smoothing iterations,
        %   * level_max  - number of grid,
        %   * smooth_m   - method of smoothing.
        function v = solveProblemByMG(obj, nu_1, nu_2, level_max, smooth_m)            
            A_mat = obj.A_op.getMatrix(obj.grid);
            g = obj.getRightPart();
            v_prev = obj.fillBorder(obj.grid);
            v_prev_prev = v_prev;
            
            % if exact solution is known 
            % it will be calculated absolute error
            if ~isempty(obj.u_ex)                 
                u0 = obj.u_ex(obj.grid.x, obj.grid.y);
            end
            
            obj.error = ones(obj.k_max, 1);
            obj.rho = zeros(obj.k_max, 1);
            k = 1;
            while k <= obj.k_max && obj.eps < obj.error(max(k-1 ,1))
                v = MultiGrid(obj, A_mat, v_prev, g, 1, nu_1, nu_2, 1,...
                              level_max, smooth_m);
                if k > 1
                    obj.rho(k) = norm(v - v_prev)/norm(v_prev - v_prev_prev);                    
                end

                if ~isempty(obj.u_ex)
                    obj.error(k) = norm(u0 - v);
                else
                    obj.error(k) = norm(v_prev - v);
                end
                v_prev_prev = v_prev;
                v_prev = v;
                
                sprintf('k = %d, err = %.2e\n', k, obj.error(k))
                k = k + 1;
            end   
            obj.error = obj.error(1:k-1);
            obj.rho = obj.rho(1:k-1);
            %param_str = sprintf('Solver_G_%d_l_%d_s_%s_nu_%d',...
            % obj.grid.N, level_max, smooth_m, nu_1);
			%save(strcat('mats/', param_str, '.mat'), 'obj');
        end
            
        % func returns solution of finite difference equation,
        % it is solving by some iteration method
        % PARAM IN:
        %   * method - name of method,
        %   * omega  - parameter of SOR method
        function v = solveProblemByIter(obj, method, omega) 
            if ~isempty(obj.u_ex)
                u0 = obj.u_ex(obj.grid.x, obj.grid.y);
            end
            
            A_mat = obj.A_op.getMatrix(obj.grid);
            v_prev = obj.fillBorder(obj.grid);
            g = obj.getRightPart();
            
            method_func = [];
            if strcmp(method, 'Jacobi')	
                method_func = @(v_prev) obj.JacobiIter(A_mat, v_prev, g);
            end
            if strcmp(method, 'Seidel')
                method_func = @(v_prev) obj.SeidelIter(A_mat, v_prev, g);
            end
            if strcmp(method, 'SOR')
                method_func = @(v_prev) obj.SORIter(A_mat, v_prev, g, omega);
            end
            if isempty(method_func)
                disp('solveProblemByIter::ERROR: undefined method')
                return
            end
                                    
            v_prev_prev = v_prev;            
            obj.error = ones(obj.k_max, 1);	
            obj.rho = zeros(obj.k_max, 1);
            k = 1;
            while k <= obj.k_max && obj.eps < obj.error(max(k-1 ,1))
                v = method_func(v_prev);
                if k > 1
                    obj.rho(k) = norm(v - v_prev)/norm(v_prev - v_prev_prev);                      
                end
                v_prev_prev = v_prev;
                v_prev = v;     
                if ~isempty(obj.u_ex)                    
                    obj.error(k) = norm(u0 - v);                
                else
                    obj.error(k) = norm(v_prev - v);                
                end
                k = k + 1;
                if ~rem(k, 100)
                    sprintf('k = %d, err = %.2e\n', k, obj.error(k-1))
                end
            end
            obj.error = obj.error(1:k-1);            
            obj.rho = obj.rho(1:k-1);
        end
        
        % func sets params of solver
        % PARAM IN:
        %   * eps   - accuracy of error,
        %   * k_max - max number of iterations
        function setParams(obj, eps, k_max)
            obj.eps = eps;
            obj.k_max = k_max;
        end
        
        % func returns smoothing solution of the problem Av = g
        % PARAM IN:
        %   * A - matrix of problem,
        %   * v - approximate solution,
        %   * g - right part of equation,
        %   * k - number of smoothing iterations,
        %   * method - method of smoothing
        v_s = smoothEll(obj, A, v, g, k, method)
    end
    
    methods (Static)        
        d_h = prolongation(d_2h)
        
        d_2h = restriction(d_h)   
        
        v = JacobiIter(A, v0, g)
        
        v = SeidelIter(A, v0, g)
        
        v = SORIter(A, v0, g, omega)
    end
    
    methods (Access = private)        
        vm_new = MultiGrid(obj, A_m, vm_old, g_m, gamma, nu_1, nu_2,...
                           level, level_max, smooth_m)                  
    end
          
end