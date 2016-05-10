function S = closed_cubic_bspline(P, varargin)
%   INPUT:
%   P: points in row vectors (N-by-3 matrix)
%   varargin(optional): 1 for plotting
%
%   Author: Di Zhu 2016-05-06
%   reference book: Curves and Surfaces for Computer Graphics, 
%   David Salomon, 2006

% check the validity of input; at least 4 points needed for degree of 3 
if size(P,1) < 4
    error('Not enough control points');
else
    np = size(P,1);
end

% modify the input of 1D/2D points into 3D points by adding zeros
if size(P,2) == 1
    P = [P,zeros(size(P,1),1),zeros(size(P,1),1)];
elseif size(P,2) == 2
    P = [P,zeros(size(P,1),1)];
else
    error('column dimension of point matrix must be 1, 2 or 3');
end 

if nargin == 2
    plot_switch = varargin{1}; % 1 for plotting
elseif nargin > 2
    error('Too many input arguments');
end
    
% define basis matrix M
M = [-1 3 -3 1 ; 3 -6 3 0 ; -3 0 3 0 ; 1 4 1 0 ]/6; % ref P258
% preallocate for speed
nj = 100;
S = cell(np,1);

for i = 1 : np % form np groups of points
    
    % each group has 4 points (P_i-1, P_i, P_i+1, P_i+2)
    if i + 3 <= np
        Pi = P(i : i + 3, :);
    else
        Pi = [ P(i : end, :); P(1 : 4 - size(P(i : end, :)), :)];
    end
    
    % calculate cubic B-Spline segment
    for j = 1 : nj % able to use parallel for loop
        t = 1/(nj-1) * j;
        B(j,:) = [t.^3, t.^2, t, 1] * M * Pi;
    end
    
    % add Bezier section to spline
    S{i} = B;
end

if plot_switch == 1
    for i = 1 : np
        plot3(S{i}(:,1),S{i}(:,2),S{i}(:,3));
        hold on;
    end
else
end


