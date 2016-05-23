function N_ik = basisfunc_NBS( N, K, T, t, fileName )
%basisfunc_NBS computes the basis function for Nonrational B-Splines
%
%   INPUT:
%   N - number of control points
%       e.g. if there are 5 control points, N = 5.
%   K - order (degree + 1) 
%       e.g. for quadratic B-Splines, degree = 2; order K = 2 + 1 = 3.
%   T - knot vector, a nondecreasing sequence from t_0 to t_(N+K-1)
%       e.g. for N = 5, K = 3, a knot vector can be:
%           [t0 t1 t2 t3 t4 t5 t6 t7]  -> parameters
%           [ 0  0  0  1  2  3  3  3]  -> nonuniform
%   t - parameter character, has to be of type char
%       e.g. 't'
%   fileName - the name of txt file to write data into
%       e.g. fileName = 'N_ik'
%
%   OUTPUT:
%   basis function N_ik{k,i}, 1 <= k <= K, 0 <= i <= N-1, t_K-1 <= t <= t_N
%       e.g. for N = 5, K = 3, basis function N_ik consist of
%       (K = 1) N_01 N_11 N_21 N_31 N_41
%       (K = 2) N_02 N_12 N_22 N_32 N_42
%       (K = 3) N_03 N_13 N_23 N_33 N_43
%       for each control point, there is a weight function in each order
%       for K = 1, N_i1(t) = 1 when t in [t_i, t_i+1); 0 otherwise
%       for K > 1, N_ik(t) are defined recursively based on equation:
%       N_ik(t) = (t - t_i)/(t_i+k-1 - t_i) * N_i,k-1(t)      % N_ik_left
%               + (t_i+k - t)/(t_i+k - t_i+1) * N_i+1,k-1(t)  % N_ik_right
%
%       Note that the indicing in matlab starts from 1 instead of 0, so in
%       order to get access to N_01, where i=0 and k=1, one should use
%       N_ik{k,i+1} or N_ik{1,1}; and to see the function at each interval,
%       use N_ik{k,i+1}{#of interval}, e.g. for N_23 over interval t=[1,2)
%       in matlab it's actually N_i=2,k=3{k=3,i+1=3} over the second interval
%       so one should put N_ik{3,3}{2} to see the value.
%
%   EXAMPLE: (book P303)
%       N = 5; K = 3; T = [0 0 0 1 2 3 3 3];
%       N_ik =  basisfunc_NBS( N, K, T );
%       Result: 3-by-5 cell
%         order k = 1: 
%         N_01 = 	0
%         N_11 = 	0
%         N_21 = 	1, t=[0,1);
%         N_31 = 	1, t=[1,2);
%         N_41 = 	1, t=[2,3);
%         order k = 2: 
%         N_02 = 	0
%         N_12 = 	1 - t, t=[0,1);
%         N_22 = 	t, t=[0,1);	2 - t, t=[1,2);
%         N_32 = 	t - 1, t=[1,2);	3 - t, t=[2,3);
%         N_42 = 	t - 2, t=[2,3);
%         order k = 3: 
%         N_03 = 	(t - 1)^2, t=[0,1);
%         N_13 = 	- t*(t/2 - 1) - t*(t - 1), t=[0,1);	(t/2 - 1)*(t - 2), t=[1,2);
%         N_23 = 	t^2/2, t=[0,1);	- (t/2 - 3/2)*(t - 1) - (t*(t - 2))/2, t=[1,2);	(t/2 - 3/2)*(t - 3), t=[2,3);
%         N_33 = 	(t/2 - 1/2)*(t - 1), t=[1,2);	- (t - 2)*(t - 3) - (t/2 - 1/2)*(t - 3), t=[2,3);
%         N_43 = 	(t - 2)^2, t=[2,3);
%
%   Author: Di Zhu 2016-05-09
%   ref: Curves and Surfaces for Computer Graphics, David Salomon, 2006.

% check inputs
if length(T) ~= N + K
    error('Wrong number of knots');
end

t = syms(t); % parameter t
N_ik = cell(K, N);      % preallocate N_ik
n_seg = N - (K - 1); % number of polynomial segments
% number of polynomial segments also indicate the number of intervals
% for example, when n_seg = 3, there are three intervals for parameter t:
% 1st interval: [0,1)
% 2nd interval: [1,2)
% 3rd interval: [2,3)
% the left bound is closed and right bound is open

% initiate the basis functions N_i1 for K = 1
ind = 0;
for i = 1 : N
    % preallocate
    N_ik{1,i} = num2cell(zeros(n_seg,1));
    
    if T(i+1) > 0
        if T(i+1) - T(i) ~= 0
            ind = ind + 1;
            N_ik{1,i}{ind} = 1; 
        end
    end
    
end 

% compute basis functions for higher order recursively
for k = 2 : K % for order k in [1,K]  
    for i = 1 : N  % for i in [0,N]
        
        % piecewise function of N_ik based on intervals
        N_ik{k,i} = cell(n_seg,1);
        
        for j = 1 : n_seg
            % determine if N_i,k-1=0
            if N_ik{k-1,i}{j} == 0
                N_ik_left = 0;
            else
                N_ik_left = (t - T(i)) / (T(i + k - 1) - T(i)) * N_ik{k-1,i}{j};
            end
            
            % determine if N_i,k=0
            if i < N 
                if N_ik{k-1,i+1}{j} == 0 
                    N_ik_right = 0;
                else
                    N_ik_right = (T(i + k) - t) / (T(i + k) - T(i + 1)) * N_ik{k-1,i+1}{j};
                end
            else
                N_ik_right = 0;
            end
      
            % combine the two components
            N_ik{k,i}{j} = N_ik_left + N_ik_right;
        end
    end 
end

% write the basis functions to txt file
fid = fopen(strcat(fileName,'.txt'),'wt'); % delete text and rewrite file
fprintf(fid, strcat(fileName,':\n'));
for k = 1 : K % order
    fprintf(fid,'order k = %d: \n',k);
    for i = 1 : N % segment
        fprintf(fid,'N_%d%d = ',i-1,k);
        var1 = 0;
        for j = 1 : n_seg
            
            if isnumeric( N_ik{k,i}{j} ) % if the value is a number
                if N_ik{k,i}{j} ~= 0 % if the function is not zero
                    fprintf(fid,'\t%d, %s=[%d,%d);', N_ik{k,i}{j},char(t),j-1,j);
                    var1 = var1 + 1;
                end
            else
                fprintf(fid,'\t%s, %s=[%d,%d);', char(N_ik{k,i}{j}),char(t),j-1,j);
                var1 = var1 + 1;
            end
        end
        
        % when the function = 0, print 0 after N_ik
        if var1 == 0
            fprintf(fid,'\t0\n');
        else
            fprintf(fid,'\n');
        end
        
    end
end
fclose(fid);
