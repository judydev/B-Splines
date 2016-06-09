function P = subCurvePoint(CP,K,N_ik,u,uu)
% with given u = uu, find the point on clamped B-Spline curve
%
%   INPUT:
%   CP - control points
%   K - order (degree + 1)
%   N_ik - basis functions, from "basisfunc_NUBS"
%   u - parameter character, consistent with char used in N_ik
%   uu - parameter, a real number given in scope of U
%
%   OUTPUT:
%   P - coordinates of the point on curve corresponding to the given uu
%
%   Author: Di Zhu 2016-06-09

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

% compute the point on curve
% find span
span = 0 : N-K+1;
for ii = 1 : N-K+1
    if uu >= span(ii) && uu < span(ii+1)
        break;
    end
end
% compute the point in corresponding interval 
P = subs(Ci{ii},u,uu);
