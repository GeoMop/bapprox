function [ ] = basetestvec( knots )
n = 100;
n_basf = length(knots)-3;
one = ones(n_basf,1);
y = sparse(zeros(n,n_basf));
for i=1:n
    y(i,:) =  splinebasevec(knots,min(knots) + max(knots)*(i-1)/(n-1));
end
% spy(y)
% pause
x = linspace(min(knots),max(knots),n);
plot(x,y)
figure
if norm(y*one - ones(n,1))< n_basf*eps
    disp('OK')
end
end

