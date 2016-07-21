function [ ] = plotresultsvec2( u_knots, v_knots,P0,P1,P2,P3,X,z, Err)

[np k] = size(X);

u_n_basf = length(u_knots)-3;
v_n_basf = length(v_knots)-3;


nu = 60;
nv = 60;

Zsurf = zeros(nu,nv);
Xsurf = zeros(nu,nv);
Ysurf = zeros(nu,nv);

X_coor = zeros(u_n_basf,v_n_basf);
Y_coor = zeros(u_n_basf,v_n_basf);


% Compute X & Y control points

for i =1:u_n_basf
    u = (u_knots(i+1)+u_knots(i+2))/2;
    for j =1:v_n_basf
        v = (v_knots(j+1)+v_knots(j+2))/2;
        P = u * (v * P2 + (1-v) * P1) + (1-u) * ( v * P3 + (1-v) * P0) ; 
        X_coor(i,j) = P(1);
        Y_coor(i,j) = P(2);   
    end
end

x = X_coor(:);
y = Y_coor(:);


% Compute fine grid in X & Y & Z for draw

uf = spalloc(u_n_basf,nu,3*nu); 
vf = spalloc(v_n_basf,nv,3*nv); 

for j =1:nu
    u = (j-1)/(nu-1);
    uf(:,j) = splinebasevec(u_knots,u,0);
end

for k =1:nv
    v = (k-1)/(nv-1);
    vf(:,k) = splinebasevec(v_knots,v,0);
end

for k =1:nv
        Zsurf(:,k) = kron(vf',uf(:,k)') * z; 
        Xsurf(:,k) = kron(vf',uf(:,k)') * x;
        Ysurf(:,k) = kron(vf',uf(:,k)') * y;
end



surf(Xsurf,Ysurf,Zsurf);
hold on
for k=1:np
    if X(k,4) ~=0
        plot3(X(k,1),X(k,2),X(k,3),'.','MarkerSize',35);
    end
end
% hold off
% 
%  figure
%  
%  Xsurferr = kron(xv(2:u_n_basf-1)',ones(v_n_basf-2,1));
%  Ysurferr = kron(yv(2:v_n_basf-1),ones(u_n_basf-2,1)');
% 
%  
% surf(Xsurferr,Ysurferr,Err);

end

