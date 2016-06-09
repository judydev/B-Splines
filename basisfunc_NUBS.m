function [ N_ik ] = basisfunc_NUBS( N, K, U, u, fileName )
%basisfunc_NUBS computes the basis function for Nonuniform B-Splines
%
%   INPUT:
%   N - number of control points
%       e.g. if there are 5 control points, N = 5.
%   K - order (degree + 1) 
%       e.g. for quadratic B-Splines, degree = 2; order K = 2 + 1 = 3.
%   U - knot vector, a nondecreasing sequence from u_0 to u_(N+K-1)
%       e.g. for N = 5, K = 3, a knot vector can be:
%           [u0 u1 u2 u3 u4 u5 u6 u7]  -> parameters
%           [ 0  0  0  1  2  3  3  3]  -> nonuniform
%   u - parameter character
%       e.g. 'u'
%   fileName - the name of txt file to write data into. or leave blank if
%              do not need to write to file
%       e.g. fileName = 'N_ik'
%
%   OUTPUT:
%   basis function N_ik{k,i}, 1 <= k <= K, 0 <= i <= N-1, u_K-1 <= u <= u_N
%       e.g. for N = 5, K = 3, basis function N_ik consist of
%       (K = 1) N_01 N_11 N_21 N_31 N_41
%       (K = 2) N_02 N_12 N_22 N_32 N_42
%       (K = 3) N_03 N_13 N_23 N_33 N_43
%       for each control point, there is a weight function in each order
%       for K = 1, N_i1(u) = 1 when u in [u_i, u_i+1); 0 otherwise
%       for K > 1, N_ik(u) are defined recursively based on equation:
%       N_ik(u) = (u - u_i)/(u_i+k-1 - u_i) * N_i,k-1(u)      % N_ik_left
%               + (u_i+k - u)/(u_i+k - u_i+1) * N_i+1,k-1(u)  % N_ik_right
%
%       Note that the indicing in matlab starts from 1 instead of 0, so in
%       order to get access to N_01, where i=0 and k=1, one should use
%       N_ik{k,i+1} or N_ik{1,1}; and to see the function at each interval,
%       use N_ik{k,i+1}{#of interval}, e.g. for N_23 over interval u=[1,2)
%       in matlab it's actually N_i=2,k=3{k=3,i+1=3} over the second interval
%       so one should put N_ik{3,3}{2} to see the value.
%
%   EXAMPLE: (book P303)
%       N = 5; K = 3; U = [0 0 0 1 2 3 3 3];
%       N_ik =  basisfunc_NUBS( N, K, U, 'u','N_ik');
%
%   Author: Di Zhu 2016-06-09
%   ref: Curves and Surfaces for Computer Graphics, David Salomon, 2006.

% check inputs
if length(U) ~= N + K
    error('Wrong number of knots');
end

u = sym(u); % parameter u
N_ik = cell(K, N+K-1);      % preallocate N_ik
n_seg = N - (K - 1); % number of polynomial segments
% number of polynomial segments also indicate the number of intervals
% for example, when n_seg = 3, there are three intervals for parameter u:
% 1st interval: [0,1)
% 2nd interval: [1,2)
% 3rd interval: [2,3)
% the left bound is closed and right bound is open

% initiate the basis functions N_i0 for K = 0 (N_i1 for K=1 in Matlab)
ind = 0;
for i = 1 : (N + K - 1)
    % preallocate
    N_ik{1,i} = num2cell(zeros(n_seg,1));
    
    if U(i+1) > 0
        if U(i+1) - U(i) ~= 0
            ind = ind + 1;
            N_ik{1,i}{ind} = 1; 
        else
            N_ik{1,i}{ind} = 0;
        end
    end
end 

% compute basis functions for higher order recursively
for k = 2 : K % for order k in [1,K]  
    for i = 1 : (N + K - k)  % for i in [0,N]
        
        % piecewise function of N_ik based on intervals
        N_ik{k,i} = cell(n_seg,1);
        
        for j = 1 : n_seg
            % determine if N_i,k-1=0
            if N_ik{k-1,i}{j} == 0
                N_ik_left = 0;
            else
                N_ik_left = (u - U(i)) / (U(i + k - 1) - U(i)) * N_ik{k-1,i}{j};
            end
            
            % determine if N_i,k=0
            if i < N 
                if N_ik{k-1,i+1}{j} == 0 
                    N_ik_right = 0;
                else
                    N_ik_right = (U(i + k) - u) / (U(i + k) - U(i + 1)) * N_ik{k-1,i+1}{j};
                end
            else
                N_ik_right = 0;
            end
      
            % combine the two components
            N_ik{k,i}{j} = N_ik_left + N_ik_right;
        end
    end 
end

if nargin == 5
    % write the basis functions to txt file
    fid = fopen(strcat(fileName,'.txt'),'wt'); % delete text and rewrite file
    fprintf(fid, strcat(fileName,':\n'));
    fprintf(fid,'Number of control points = %d\n',N);
    fprintf(fid,'Order K = %d, degree p = %d\n',K,K-1);
    fprintf(fid,'Knot vector T = {');
    for i = 1 : length(U)
        fprintf(fid,' %d',U(i));
    end
    fprintf(fid,' }\n\n');
    
    for k = 1 : K % order
        fprintf(fid,'degree p = %d: \n',k-1);
        for i = 1 : ( N + K - k) % segment
            fprintf(fid,'N_%d%d = ',i-1,k-1);
            var1 = 0;
            for j = 1 : n_seg
                
                if isnumeric( N_ik{k,i}{j} ) % if the value is a number
                    if N_ik{k,i}{j} ~= 0 % if the function is not zero
                        fprintf(fid,'\t%d, %s=[%d,%d);', N_ik{k,i}{j},char(u),j-1,j);
                        var1 = var1 + 1;
                    end
                else
                    fprintf(fid,'\t%s, %s=[%d,%d);', char(N_ik{k,i}{j}),char(u),j-1,j);
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
end
