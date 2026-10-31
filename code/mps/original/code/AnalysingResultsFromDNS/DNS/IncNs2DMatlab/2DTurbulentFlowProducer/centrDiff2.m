function dv = centrDiff2(v,h, dir)
%computes 9-stencil 2nd order diff

coeffs = [-205/72 8/5 -1/5 8/315 -1/560]/h^2;
%coeffs = [-2 1 0 0 0]/h^2;


dv = coeffs(1)*v;
dv = dv + coeffs(2)*( +circshift(v,1,dir) + circshift(v,-1,dir));
dv = dv + coeffs(3)*( +circshift(v,2,dir) + circshift(v,-2,dir));
dv = dv + coeffs(4)*( +circshift(v,3,dir) + circshift(v,-3,dir));
dv = dv + coeffs(5)*( +circshift(v,4,dir) + circshift(v,-4,dir));

end

