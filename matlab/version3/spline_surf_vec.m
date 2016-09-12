clc
clear all
close all

% knots = [0 0 0 1/3 2/3  1 1 1]';
% %knots = [0 0 0 1/2  1 1 1]';
% 
% u_knots = knots;
% v_knots = knots;
% u_n_basf = size(knots-3);
% v_n_basf = size(knots-3);

u_n_basf = 15;
v_n_basf = 19;
u_knots = get_knot_vector(u_n_basf);
v_knots = get_knot_vector(v_n_basf);

%%% base test

%basetestvec(v_knots);

%%% patch boundary

P0 = [0.2; 0.3];
P1 = [1.2; 0.5];
P2 = [1.5; 2];
P3 = [0.9; 1.7];

%%% Reduce points to be approximatex (i.e., points which are in convex hull of the points P0,P1,P2,P3)

X = getpoints();
[Xp X] = transform2( X,P0,P1,P2,P3 );
[np k] = size(Xp);
[ Xp ] = solve2( Xp,P0,P1,P2,P3 )

%%% Construction of the matrix 

%interv ??
[B, Interv ] = build_LS_matrix( u_knots,v_knots, Xp);

W = sparse(diag(Xp(:,4)));

g = Xp(:,3);
b = B'*W*g;

C = B'*W*B;
nnzC = nnz(C);


 A = build_reg_matrix( u_knots,v_knots, P0,P1,P2,P3,nnzC);
 
 nC = norm(full(C))
 nA = norm(full(A))
  r = norm(full(C))/  norm(full(A))
 %r = 0.0;
 
 S = C+0.3*r*A;

z = pcg(A,b,1e-12,500);

%%% Solution

% Direct Solver
% C = B'*W*B;%+A;
% b = B'*W*g;
% z = C\b;


% Errors
Err = get_errors(abs(W*B*z-g),Interv,u_n_basf,v_n_basf);


%%% PLOT results

plotresultsvec2(u_knots, v_knots,P0,P1,P2,P3,X,z,Err)



