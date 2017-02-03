function [ Xp X] = transform( X,xp,yp )

[np k] = size(X);

h = 0;
for j=1:np
    if ((X(j,1) >= xp(1)) && (X(j,1) <= xp(2)) && (X(j,2) >= yp(1)) && (X(j,2) <= yp(2)))
      h = h+1;
      Xp(h,1) = (X(j,1)-xp(1))/(xp(2) -xp(1)) ;
      Xp(h,2) = (X(j,2)-yp(1))/(yp(2) -yp(1)) ;
      Xp(h,3:4) = X(j,3:4);
    else
        X(j,3:4) = 0;
    end
end

end

