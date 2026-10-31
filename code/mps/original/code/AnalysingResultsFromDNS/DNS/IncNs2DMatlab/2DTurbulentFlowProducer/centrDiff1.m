function dv = centrDiff1(v,h, dir)
%computes 9-stencil 1st order diff

coeffs = [0 4/5 -1/5 4/105 -1/280]/h;
%coeffs = [0 1/2 0 0 0]/h;

dv = coeffs(1)*v;
dv = dv + coeffs(2)*( -circshift(v,1,dir) + circshift(v,-1,dir));
dv = dv + coeffs(3)*( -circshift(v,2,dir) + circshift(v,-2,dir));
dv = dv + coeffs(4)*( -circshift(v,3,dir) + circshift(v,-3,dir));
dv = dv + coeffs(5)*( -circshift(v,4,dir) + circshift(v,-4,dir));

end

