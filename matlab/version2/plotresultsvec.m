function [ ] = plotresultsvec( u_knots, v_knots,xp,yp,X,z , Err)

[np k] = size(X);

u_n_basf = length(u_knots)-3;
v_n_basf = length(v_knots)-3;


xv = (u_knots(2:u_n_basf+1)+u_knots(3:u_n_basf +2))/2*(xp(2)-xp(1))+xp(1);
yv = (v_knots(2:v_n_basf+1)+v_knots(3:v_n_basf +2))/2*(yp(2)-yp(1))+yp(1);


nu = 60;
nv = 60;

Zsurf = zeros(nu,nv);
Xsurf = zeros(nu,nv);
Ysurf = zeros(nu,nv);


uf = spalloc(u_n_basf,nu,3*nu); 
vf = spalloc(v_n_basf,nv,3*nv); 

for j =1:nu
    u = (j-1)/(nu-1);
    uf(:,j) = splinebasevec(u_knots,u);
end

for k =1:nv
    v = (k-1)/(nv-1);
    vf(:,k) = splinebasevec(v_knots,v);
end


for k=1:nu
    Zsurf(:,k) =  kron(vf',uf(:,k)')*z; % I_vn
end


Xsurf = kron(xv'*uf,ones(nv,1));
Ysurf = kron(ones(1,nu),vf'*yv);


surf(Xsurf,Ysurf,Zsurf);
hold on
for k=1:np
    if X(k,4) ~=0
        plot3(X(k,1),X(k,2),X(k,3),'.','MarkerSize',35);
    end
end
hold off

 figure
 
 Xsurferr = kron(xv(2:u_n_basf-1)',ones(v_n_basf-2,1));
 Ysurferr = kron(yv(2:v_n_basf-1),ones(u_n_basf-2,1)');

 
surf(Xsurferr,Ysurferr,Err);

end

