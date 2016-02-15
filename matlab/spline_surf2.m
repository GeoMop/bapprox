clc
clear all
close all

% knots = [0 0 0 1/3 2/3  1 1 1]';
% knots = [0 0 0 1/2  1 1 1]';

n_basf = 15;
knots = zeros(n_basf+5, 1);
knots(3:n_basf+3) = linspace(0, 1, n_basf+1);
knots(n_basf+4:n_basf+5) = 1;

u_knots = knots;
v_knots = knots;

%%% base test

%basetest2(knots);

%%% patch boundary

xp = [0.1; 0.9];
yp = [0.1; 1.9];

%%% Interpolating points

X = getpoints();
Xp = transform(X, xp, yp);
[np k] = size(Xp);

%%% Construction of the matrix 

u_n_basf = length(u_knots) - 3;
v_n_basf = length(v_knots) - 3;

B = zeros(np, u_n_basf * v_n_basf);


uf = zeros(u_n_basf, 1);
vf = zeros(v_n_basf, 1);

for j = 1:np
    for k =1:u_n_basf
        uf(k) = splinebase2(u_knots, k, Xp(j,1));
    end
    for k =1:v_n_basf
        vf(k) = splinebase2(v_knots, k, Xp(j,2));
    end
    B(j,:) = vf'*kron(eye(v_n_basf),uf');    
end

g = Xp(:,3)

%%% Solution

% Unstable
% W = diag(Xp(:,4));
% A = B'*W*B;
% b = B'*W*g;
% z = A\b

% Stable
[q r] = qr(B);
z = r\q'*g;
 
%%% PLOT results

plotresults(u_knots, v_knots, xp, yp, Xp, z)

%TODO
% chyby