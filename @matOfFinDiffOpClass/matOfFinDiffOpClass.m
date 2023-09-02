% Class 'matOfFinDiffOpClass' implements a matrix of finite difference 
% operator that approximates by five-nodes template, i.e.
% a v_{i, j} + d v_{i-1, j} + c v_{i, j-1} + e v_{i, j+1} + b v_{i+1, j}
classdef matOfFinDiffOpClass < handle
    
    properties
        % matrix represents by 5 vectors
        a = [];
        b = [];
        c = [];
        d = [];
        e = [];
        
        N = [];
    end
       
    methods (Access = public)
        function obj = matOfFinDiffOpClass(a, b, c, d, e)
            if numel(a) == numel(b) &&...
                numel(a) == numel(c) &&...
                numel(a) == numel(d) &&...
                numel(a) == numel(e)
            obj.a = a;
            obj.b = b;
            obj.c = c;
            obj.d = d;
            obj.e = e;   
            
            [~, obj.N] = size(a);
            else
                disp('matOfFinDiffOp::ERROR: Sizes of vectors is not equal')
            end
        end
        
        function product = mul(obj, v)
            product = zeros(size(v));
            
            [~, N_v] = size(v);
            if N_v - 2 ~= obj.N
                disp('matOfFinDiffOp::mul::ERROR: Index dimensions mismatch')
                return
            end     
            % finite difference operator don't change values in border nodes
            % for closing system of equation to solve it
            product(1, :)   = v(1,   :);
            product(end, :) = v(end, :);
            product(:, 1)   = v(:,   1);
            product(:, end) = v(:, end);
            % indices of inner nodes
            i_in = 2:obj.N+1;
            j_in = 2:obj.N+1;
            product(j_in, i_in) = obj.d.*v(j_in,   i_in-1)+...
                                  obj.c.*v(j_in-1, i_in)+...
                                  obj.a.*v(j_in,   i_in)+...
                                  obj.e.*v(j_in+1, i_in)+...
                                  obj.b.*v(j_in,   i_in+1);
        end
        
        function v = solveSystem(obj, g)
            if obj.N ~= obj.N
                disp('matOfFinDiffOp class::solveSistemOneInner::ERROR: number of inner nodes do not equals 1')
                v = 0;
                return
            end
            v = obj.TriDiagMatrixAlgorithm(g);
            %A = getFullMat();
            %v = reshape(A\g(:), obj.N + 2, obj.N + 2);
        end        
    end 
    
    methods (Access = private)
        v = TriDiagMatrixAlgorithm(obj, g)
        A = getFullMat(obj)
    end
    
end

