Ekin(t) = sum(sum((U.^2 + V.^2))) %*0.5*dx*dy/L^2;

disp(sum(sum((U.^2 + V.^2))))

dudx = (real(ifft2( i*KX.*fft2(U)))); dudy = (real(ifft2( i*KY.*fft2(U)))); 
dvdx = (real(ifft2( i*KX.*fft2(V)))); dvdy = (real(ifft2( i*KY.*fft2(V)))); 

Diss(t) = nu*sum(sum((dudx.^2 + dudy.^2 + dvdx.^2 + dvdy.^2)*dx*dx))/L^2;
