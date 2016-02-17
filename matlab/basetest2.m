function [ ] = basetest2( knots )
n = 100;
n_basf = length(knots)-3;
one = ones(n_basf, 1);
y = zeros(n, n_basf);
for k =1:n_basf
    for i=1:n
        y(i,k) = splinebase2(knots, k, min(knots) + max(knots)*(i-1)/(n-1));
    end
end
x = linspace(min(knots), max(knots), n);

plot(x,y)
figure

if norm(y*one - ones(n,1)) == 0
    disp('OK')
end
