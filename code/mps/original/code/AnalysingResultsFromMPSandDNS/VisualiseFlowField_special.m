%% Visualises flow field

function [fieldCurlMax, fieldCurlMin] = VisualiseFlowField_special(U,V,sigma,showAxis,fontsz,linWidth,showQuiver,showContour)

% %Only look at the fluctuations
%U = U-mean(U,2);
%V = V-mean(V,2);

%Assume spacing is identical in both x and y direction:
[Nx,Ny] = size(U);
x = (0:(Nx-1))/Nx;
y = (0:(Ny-1))/Ny;
dx = 1/Nx;

%div = sqrt(sum(sum((divergence(U,V).^2))))/(Nx*Ny)
%Ek = sqrt(sum(sum(U.^2+V.^2)))/(Nx*Ny)
if(sigma > 0)
    U = imgaussfilt(U,sigma);
    V = imgaussfilt(V,sigma);
end

dVdx = centrDiff1(V,dx,2);
dUdy = centrDiff1(U,dx,1);
Ek = (U.^2+V.^2);
Urms = Ek.^0.5;

y0 = 0.198;%.197;%.25;
y1 = .802;%.803;%.75;

%figure(1), box
% clf, contourf(x,y,(curl(U,V)).^1,10,'white-'), hold on;
fieldCurl = -dUdy+dVdx;

fieldCurl(abs(fieldCurl)<5)=NaN;

if(showContour==1)
    img= imagesc(x, y,fieldCurl); set(gca,'YDir','normal'); hold on;
end
% colorbar
Len = 2;%sqrt(U.^2+V.^2+eps);
if(showQuiver==1)
    arrowSpacing = 24; %48 is high res
    d = round(Nx/arrowSpacing);
    startingSpace_x = 1;
    endingSpace_x = 1;
    startingSpace_y = 1;
    endingSpace_y = 1;
    
    quiver(x(startingSpace_x:d:end-endingSpace_x),y(startingSpace_y:d:end-endingSpace_y) ...
    ,U(startingSpace_y:d:end-endingSpace_y,startingSpace_x:d:end-endingSpace_x)./Len./ ...
    Urms(startingSpace_y:d:end-endingSpace_y,startingSpace_x:d:end-endingSpace_x), ... 
    V(startingSpace_y:d:end-endingSpace_y,startingSpace_x:d:end-endingSpace_x)./Len./ ... 
    Urms(startingSpace_y:d:end-endingSpace_y,startingSpace_x:d:end-endingSpace_x),.5,'color',[0.5 0.5 0.5]);

end
hold off, axis equal, axis([x(1) x(end) y(round(1+Ny*y0)) y(round(end*y1))])
set(img,'AlphaData',double(~isnan(fieldCurl)))
%drawnow

% Axes
xticks([0.2 0.5 .8])
yticks([0.3 .5 .7])

showAxis_x=showAxis(1);
showAxis_y=showAxis(2);
showLabels=showAxis(3);

if(showAxis_x==1 && showAxis_y==1)
    xticklabels({'0.2','0.5','0.8'});
    yticklabels({'0.3','0.5','0.7'});
elseif(showAxis_x==1 && showAxis_y==0)
    xticklabels({'0.2','0.5','0.8'});
    yticklabels({'','',''});ylabel('')
elseif(showAxis_x==0 && showAxis_y==1)
    yticklabels({'0.3','0.5','0.7'});
    xticklabels({'','',''});xlabel('')
else
    xticklabels({'','',''});
    yticklabels({'','',''});
    xlabel(''); ylabel('');
end

if(showLabels==1)
       xlabel('$x/L_{\mathrm{box}}$','interpreter','latex'); ylabel('$y/L_{\mathrm{box}}$','interpreter','latex');
end
set(gca,'FontSize',fontsz)
set(gca,'TickLabelInterpreter', 'latex');

set(gca,'LineWidth', linWidth);

fieldCurlMax = max(fieldCurl(:));
fieldCurlMin = min(fieldCurl(:));

end