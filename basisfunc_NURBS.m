function R_ik = basisfunc_NURBS( N, K, T, w )
%basisfunc_NURBS computes the basis function for Nonuniform Rational B-Splines
%
%   INPUT:
%   N - number of control points
%       e.g. if there are 5 control points, N = 5.
%   K - order (degree + 1) 
%       e.g. for quadratic B-Splines, degree = 2; order K = 2 + 1 = 3.
%   T - knot vector, a nondecreasing sequence from t_0 to t_(N+K-1)
%       e.g. for N = 4, K = 3, a knot vector can be:
%           [t0 t1 t2 t3 t4 t5 t6 t7]  -> parameters
%           [ 0  0  0  1  2  3  3  3]  -> nonuniform
%   w - weight factors, w_i is the weight(>=0) for each control point
%
%   OUTPUT:
%   basis function R_ik{k,i}, 1 <= k <= K, 0 <= i <= N-1, t_K-1 <= t <= t_N
%       e.g. for N = 5, K = 3, basis function R_ik consist of
%       (K = 1) R_01 R_11 R_21 R_31 R_41
%       (K = 2) R_02 R_12 R_22 R_32 R_42
%       (K = 3) R_03 R_13 R_23 R_33 R_43
%       for each control point, there is a weight function in each order
%       for K = 1, R_i1(t) = 1 when t in [t_i, t_i+1); 0 otherwise
%       for K > 1, R_ik(t) are defined recursively based on equation:
%       R_ik(t) = w_i * N_ik(t)             % numerator 
%               / ( sum(w_i * N_ik(t)) )    % denominator
%       where N_ik(t) can be computed by N_ik = basisfunc_NBS( N, K, T )
%
%       Note that the indicing in matlab starts from 1 instead of 0, so in
%       order to get access to R_01, where i=0 and k=1, one should use
%       R_ik{k,i+1} or R_ik{1,1}; and to see the function at each interval,
%       use R_ik{k,i+1}{#of interval}, e.g. for R_23 over interval t=[1,2)      
%
%   EXAMPLE: (book P303)
%       N = 5; K = 3; T = [0 0 0 1 2 3 3 3]; w = [1 1 .5 1 1];
%       N_ik =  basisfunc_NURBS( N, K, T, w );
%       Result: 3-by-5 cell
%
%   Author: Di Zhu 2016-05-10
%   ref: Curves and Surfaces for Computer Graphics, David Salomon, 2006.

n_seg = N - (K - 1); % number of polynomial segments
% number of polynomial segments also indicate the number of intervals
% for example, when n_seg = 3, there are three intervals for parameter t:
% 1st interval: [0,1)
% 2nd interval: [1,2)
% 3rd interval: [2,3)
% the left bound is closed and right bound is open

% check input variable w
if length(w) ~= N
    error('Length of w should be the same as number of control points N'); 
end

% compute the nonrational basis functions
try
    N_ik = basisfunc_NBS( N, K, T );
catch
    error('Wrong input for N, K, T'); % when unable to compute N_ik, return error message
end

% compute the denominator for rational basis function
denom = num2cell(zeros(n_seg,1)); % preallocate denom with cells of zeros
for i = 1 : N % loop through each control points, sum each w * N_ik
    for j = 1 : n_seg % loop through each interval [0,1),[1,2),...
        denom{j} = denom{j} + w(i) * N_ik{K,i}{j} ;
        % only the Kth order N_ik will be used to compute R_iK
        % denom{:} contains the expressions in each interval
    end
    
end

% preallocate R_ik, each row indicates a weighted function for each control
% point, and each column indicates the basis function in corresponding
% interval (from 1 to n_seg => [0,1),[1,2),...,[n_seg-1,n_seg).)
R_ik = num2cell(zeros(N,n_seg)); 

% compute basis function R_0K,      R_1K, ...,      R_(N-1)K; 
% or in matlab indicing, R_ik{K,1}, R_ik{K,2},...,  R_ik{K,N}.
for i = 1 : N
    for j = 1 : n_seg
        R_ik{i,j} =  w(i) * N_ik{K,i}{j} / denom{j} ;
        if ~isnumeric(R_ik{i,j})
            R_ik{i,j} = simplify(R_ik{i,j});
        end
    end
end

% print the basis functions
Y_N = input('Print R_ik?(1 for yes, anything else for no)\n');
if Y_N == 1
    for i = 1 : N % segment
        fprintf('R_%d%d =\n ',i-1,K);
        var = 0;
        for j = 1 : n_seg
            if R_ik{i,j} ~= '0' % if the expression is not zero
                fprintf('\t\t%s, t=[%d,%d)\n', char(R_ik{i,j}), j-1,j);
                var = var + 1;
            end
        end

        % when the function = 0, print 0 after N_ik
        if var == 0
            fprintf('\t\t0\n');
        else
            fprintf('\n');
        end  
        
    end
end    

