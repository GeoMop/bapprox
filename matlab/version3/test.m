
clc
clear all
close all

P0 = [0.2; 0];
P1 = [2; 1.5];
P2 = [1.5; 2];
P3 = [0.5; 1.7];

P  = [P0 P1 P2 P3 P0];

[ X ] = getpoints();

[ Xp X] = transform2( X,P0,P1,P2,P3 );

[n m] = size(Xp);

Xp;

Xp= solve2( Xp,P0,P1,P2,P3 );


plot(P(1,:),P(2,:));
hold on
for j =1:n
   plot(Xp(j,5),Xp(j,6),'ro');
end

%compute center

 u_n_basf = 10;
 v_n_basf = 10;
 
u_knots  = get_knot_vector( u_n_basf );

v_knots  = get_knot_vector( v_n_basf );


X_coor = zeros(u_n_basf,v_n_basf);
Y_coor = zeros(u_n_basf,v_n_basf);

for i =1:u_n_basf
    u = (u_knots(i+1)+u_knots(i+2))/2;
    for j =1:v_n_basf
        v = (v_knots(j+1)+v_knots(j+2))/2;
        %P = v * ( u * P2 + (1-u) * P3 ) +  (1-v)*( (1-u) * P0 + u * P1); 
        P = u * (v * P2 + (1-v) * P1) + (1-u) * ( v * P3 + (1-v) * P0) ; 
        X_coor(i,j) = P(1);
        Y_coor(i,j) = P(2);
          
    end
end


% A = [P0' P0(1)*P0(2) 1;P1' P1(1)*P1(2) 1;P2' P2(1)*P2(2) 1;P3' P3(1)*P3(2) 1]
% B = [0 0 0 1; 1 0 0 1; 0 1 0 1; 1 1 0 1];
% 
% inv(A)*B

xc = X_coor(:);
yc = Y_coor(:);

plot(xc',yc','xk')