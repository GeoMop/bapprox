
clc
clear all
close all

P0 = [0.2; 0];
P1 = [2; 0.5];
P2 = [1.5; 2];
P3 = [0.5; 1.7];

P  = [P0 P1 P2 P3 P0];

[ X ] = getpoints();

[ Xp X] = transform2( X,P0,P1,P2,P3 );

[n m] = size(Xp);

Xp


plot(P(1,:),P(2,:));
hold on
for j =1:n
   plot(Xp(j,5),Xp(j,6),'ro');
end


