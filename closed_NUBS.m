function [C] = closed_NUBS(CP,K,U,u,N_ik)
%closed_NUBS computes and plots a closed nonuniform nonrational B-Spline
%based on given control points, order and knot vector.
%
%   INPUT:
%       CP - control points (N-by-3 matrix), make sure CP_n=CP_0
%       K - order (degree + 1)
%       U - knot vector
%           e.g. U = [zeros(1,K),1:(N-K),(N-K+1)*ones(1,K)]; -> uniform 
%       u - character for parameter, e.g. 'u'
%       fileName - for storing basis functions
%
%   OUTPUT:
%       C - Piecewise Parametric Curve
%       N_ik - basis functions
%
%   SUBFUNCTIONS:
%       [ N_ik ] = basisfunc_NUBS( N, K, U, u, fileName );
%       
%   Author: Di Zhu 2016-06-08

% %let Pn = P0. ignore if P is closed already
% if CP(1,:) == CP(end,:) 
%     CP = CP;
% else
%     CP(end+1,:) = CP(1,:); 
% end

N = length(CP);
Ci = num2cell(zeros(N-K+1,1));

for i = 1 : N-K+1 % 1 2 3 4... n_seg
    for j = 1 : N
        Ci{i} = Ci{i} + N_ik{K,j}{i} * CP(j,:);     
    end
end

% computing the points on curve with step = 0.1
ind = 0;
for i = 1 : N-K+1
    ui = K+i-1;
    for k = 1 : 10 % step = 0.1
        uik = U(ui) + (U(ui+1)-U(ui))*(k-1)/10; 
        ind = ind + 1;
        C(ind,:) = subs(Ci{i},u,uik); % store point coords for plotting
    end
end
C(ind+1,:) = C(1,:); % add the first point as ending point to close the curve

% for plotting
% m = mean(P);
% figure,hold on,daspect([1 1 1])
% plot(C(:,1),C(:,2));
% scatter(m(:,1),m(:,2))
% scatter(P(:,1),P(:,2))
