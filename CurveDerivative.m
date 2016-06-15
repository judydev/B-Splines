function PD = CurveDerivative(CP,U,u,uu)
% clamped B-Spline curve derivative
%
%   INPUT:
%   CP - control points (needs to be closed, or Pn = P0)
%   U - clamped knot vector
%   u - parameter character, consistent with char used in N_ik
%   uu - parameter, a real number given in scope of U
%
%   OUTPUT:
%   PD - coordinates of the point on curve corresponding to the given uu
%
%   Author: Di Zhu 2016-06-09

N = length(CP); % number of control points
K = length(U) - N; 
U_ = U; U_(1) = []; U_(end) = []; % U'
N_ = N - 1;
K_ = K - 1;
N_ik_ = basisfunc_NUBS(N_, K_, U_, u);

Ci_ = num2cell(zeros(N-K+1,1));

% find Qi
Q = zeros(N_,size(CP,2));
for i = 1 : N_
    Q(i,:) = (K-1) *(CP(i+1,:) - CP(i,:)) / (U(i+K) - U(i+1));
end

% find C'(u) in new intervals
for i = 1 : N_-K_+1 % 1 2 3 4... n_seg
    for j = 1 : N_
        Ci_{i} = Ci_{i} + N_ik_{K_,j}{i} * Q(j,:);     
    end
end

% compute the derivative at point on curve
% find span
span = 0 : N_-K_+1;
for ii = 1 : N_-K_+1
    if uu >= span(ii) && uu < span(ii+1)
        break;
    end
end
% compute the point in corresponding interval 
PD = subs(Ci_{ii},u,uu);
