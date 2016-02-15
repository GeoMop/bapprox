function [ ] = plotresults( u_knots, v_knots, xp, yp, X, z )

[np k] = size(X);

u_n_basf = length(u_knots)-3;
v_n_basf = length(v_knots)-3;


nu = 30;
nv = 30;


Zsurf = zeros(nu,nv);

for j =1:nu
    for k =1:nv
        
        uf = zeros(u_n_basf, 1);
        vf = zeros(v_n_basf, 1);
        u = (j-1)/(nu-1);
        v = (k-1)/(nv-1);
        
        for l =1:u_n_basf
            uf(l) = splinebase2(u_knots, l,u);
        end
        for l =1:v_n_basf
            vf(l) = splinebase2(v_knots, l,v);
        end
        
        Zsurf(k,j) = vf'*kron(eye(v_n_basf),uf')*z;

    end
end

uplot = linspace(xp(1), xp(2), nu);
vplot = linspace(yp(1), yp(2), nv);

Xsurf = kron(uplot, ones(nu,1));
Ysurf = kron(vplot', ones(nv,1)');

surf(Xsurf, Ysurf, Zsurf);
hold on
for k=1:np
    plot3(X(k,1)*(xp(2)-xp(1))+xp(1),X(k,2)*(yp(2)-yp(1))+yp(1),X(k,3),'.','MarkerSize',35);
end
hold off

end

