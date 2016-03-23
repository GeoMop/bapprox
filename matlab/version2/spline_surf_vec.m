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

xp = [0.0;1.0];
yp = [0.0;2.0];


%%% Interpolating points

X = getpoints();
[Xp X] = transform(X,xp,yp);
[np k] = size(Xp);

%%% Construction of the matrix 

%interv ??
[B, Interv ] = build_LS_matrix( u_knots,v_knots, Xp);
%A = build_reg_matrix( u_knots,v_knots, Xp);

W = sparse(diag(Xp(:,4)));

g = Xp(:,3);

%%% Solution

% Unstable
% C = B'*W*B;%+A;
% b = B'*W*g;
% z = C\b;

% Stable
% [q r] = qr(B);
% z = r\q'*g;

% SVD
%[U S V] = svd(full(B),'econ');
[U S V] = svd(full(B));
[a b] = size(S);
Sm = S > 1e-3;
S =S.*Sm;
r = rank(S);
Sdi = diag(diag(S(1:r,1:r)).^(-1));
Si = sparse(zeros(a,b));
Si(1:r,1:r) = Sdi;
z = V*Si'*U'*g;
% norm(full(B)- U*S*V')
% pause


% Errors
Err = get_errors(abs(W*B*z-g),Interv,u_n_basf,v_n_basf);


%%% PLOT results

plotresultsvec(u_knots, v_knots,xp,yp,X,z,Err)



