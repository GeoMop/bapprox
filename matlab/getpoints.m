function [ X ] = getpoints()

nu = 60;
nv = 60;
h = 0;
X = zeros(nu*nv,4);

for i=1:nu
    for j = 1:nv
        h = h+1;
        x = 2*(i-1)/(nu-1);
        y = 2*(j-1)/(nv-1);
        X(h,1) = x;
        X(h,2) = y;
        X(h,3) = cos(10*pi*x)*cos(3*pi*y)+exp(2*sin(x*y));
        %X(h,3) = cos(3*pi*exp(-x))*cos(3*pi*y^2)+exp(2*sin(x*y));
        %X(h,3) = cos(3*pi*exp(-x))*cos(3*pi*y^1.2) - 0.5*rand(1,1);%+exp(2*sin(x*y));
        X(h,4) = 1;
    end
end

end

