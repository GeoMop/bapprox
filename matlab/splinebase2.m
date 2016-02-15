function [ f ] = splinebase2( T, k, t )

%  T - knot vector
%  k - index of basis function, k = 1,...,length(T)-3 
%  t - parameter t \in [0,1]
%  f - function value

if (t < T(k)) || (t > T(k+3)) % test of support
    f = 0;
    return
end

n = length(T);
N = zeros(n-1, 1); % vector of supports 

for i = 1:n-1
    if (T(i)-T(i+1)) ~= 0
        N(i) = 1;
    end
end

tk = T(k);
tk1 = T(k+1);
tk2 = T(k+2);
tk3 = T(k+3);

if (t >= T(k)) && (t <= T(k+1)) && (N(k) ~=0)
    f = (t-tk)^2 / ((tk2 - tk) * (tk1-tk));
    return
elseif (t >= T(k+1)) && (t <= T(k+2)) && (N(k+1) ~=0)
    f= ((t-tk) * (tk2 -t)) / ((tk2-tk) * (tk2-tk1)) + ((t-tk1) * (tk3 -t)) / ((tk3-tk1) * (tk2-tk1)); 
    return
elseif (t >= T(k+2)) && (t <= T(k+3)) && (N(k+2) ~=0)
    f = (tk3-t)^2 / ((tk3-tk1)*(tk3-tk2));
    return
end


end