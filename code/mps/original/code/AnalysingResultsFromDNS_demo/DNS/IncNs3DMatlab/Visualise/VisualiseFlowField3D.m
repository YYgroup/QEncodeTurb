%% Visualises flow field
% 
%Assume spacing is identical in both x and y direction:
[Nx,Ny,Nz] = size(U);
x = (0:(Nx-1))/Nx;
y = (0:(Ny-1))/Ny;
z = (0:(Nz-1))/Nz;
[X,Y,Z] = meshgrid(x,y,z);


xmax = max(X(:));
xmin = min(X(:));
ymax = max(Y(:));
zmax = max(Z(:));

num = min(4,Nx);
range = 1:num:Nx;

cx = X(range,range,range)/xmax;
cy = Y(range,range,range)/ymax;
cz = Z(range,range,range)/zmax;

cu = U(range,range,range);
cv = V(range,range,range);
cw = W(range,range,range);

%%%%%
s = round(Nz/2);
vorticityZ = curl(U(:,:,s),V(:,:,s));

ekin = (U.^2+V.^2+W.^2)*0.5*dx^3;%/L^3;
%vortmean = sum(sum(sum(vorticity)))./(N^3);

title1 = sprintf('Velocity Field and kinetic energy\nRe = %i    t = %0.2g   <E_t/E_0> = %0.3g',Re,t, sum(sum(sum( ekin )))/Ekin(1) );
title2 = sprintf('Vorticity_z at z = 0.5\nRe = %i   t = %0.2g',Re,t);

clf;
subplot(3,1,1)
view(-8,37)
quiver3(cx,cy,cz,cu,cv,cw,'r')
hold on;
hsurfaces = slice(x,y,z,ekin,[0,max(x(:))],max(y(:)),0);
set(hsurfaces,'FaceColor','interp','EdgeColor','none')
title(title1);
xlabel('x'); ylabel('y'); zlabel('z');
hold off;

subplot(3,1,2)
imagesc(x,y,vorticityZ);
title(title2);
xlabel('x'); ylabel('y');

subplot(3,1,3)
FU=fftn(U);FV=fftn(V);FW=fftn(W);
[K,EK] = energyspec3D(FU,FV,FW,KX,KY,KZ);
loglog(K,EK,'*-');
axis([ 0 max(K) 0 1])

%figure(3)
%isosurface(X,Y,Z,vorticity3D,vortmean);
%title(sprintf('Curl of velocity\nRe = %0.1g   t = %0.2g',Re,t));

% 
% 
% 
% 
% 
% 
% box
% 
%         %clf, contourf(x,y,curl(U,V),20,'w-'), hold on;
%         imagesc(x,y,vorticity); set(gca,'YDir','normal');
%         
%         
%         colorbar('Ticks',[.75*min(vorticity(:) ),.75*max(vorticity(:))],'TickLabels',{'Low', 'High'})
%         
%      
%         Len = 1;%sqrt(U.^2+V.^2+eps);
%        % quiver(x,y,(U./Len),(V./Len),.4,'k-')
%         hold off, axis equal, axis([x(1) x(end) y(1) y(end)])
%         %title(sprintf('Curl of velocity\nRe = %0.1g   t = %0.2g',Re,t))
%         title('b                                    ');%sprintf('Vorticity;    Re = %0.0f, t = %0.0f',Re,t))
%         %xlabel('Normalised X'); ylabel('Normalised Y');
%         xlabel('x_1'); ylabel('x_2');
%         drawnow
%         
%         set(gca,'FontSize',30,'fontWeight','bold')
% 
%         set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')