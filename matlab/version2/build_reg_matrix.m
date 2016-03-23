function [ A ] = build_reg_matrix( u_knots,v_knots, Xp )

u_n_basf = length(u_knots)-3;
v_n_basf = length(v_knots)-3;
u_n_inter =  length(u_knots)-5;
v_n_inter =  length(v_knots)-5;

A = zeros(u_n_basf*v_n_basf);
for j =1:u_n_basf
    for k =1:v_n_basf
        for i= 1:u_n_inter
                Tu  = (u_knots(i+2)+u_knots(i+3))/2;
                Du = u_knots(i+3)- u_knots(i+2);
                u1 = Tu - sqrt(3)/6 * Du;
                u2 = Tu + sqrt(3)/6 * Du;
                fu1 = splinebase2(u_knots,j,u1);
                fu2 = splinebase2(u_knots,j,u2);
                fdu1 = splinebasediff2(u_knots,j,u1);
                fdu2 = splinebasediff2(u_knots,j,u2);    
                
            for l =v_n_inter
                Tv  = (v_knots(l+2)+v_knots(l+3))/2;
                Dv = v_knots(l+3)- v_knots(l+2);
                v1 = Tv - sqrt(3)/6 * Dv;
                v2 = Tv + sqrt(3)/6 * Dv;
                fv1 = splinebase2(v_knots,k,u1);
                fv2 = splinebase2(v_knots,k,u2);
                fdv1 = splinebasediff2(v_knots,k,u1);
                fdv2 = splinebasediff2(v_knots,k,u2);
                
               
                A(j+(k-1)*u_n_basf,k+(j-1)*v_n_basf) = A(j,k)+ fdu1*fv1+fdu1*fv2+fdu2*fv1+fdu2*fv2   +  fdv1*fu1+fdv1*fu2+fdv2*fu1+fdv2*fu2;
                %
            end
        end
    end
end

% A = A + A'
% 
% spy(A)
% pause
end

