function subplotEntanSpec(SingVals,chi_crit,contourLevels,fontsz,linWid,xmajticks,xmajticklabels,xminticks,yTickLabs,showaxes_x,showaxes_y,showlabel_x,showlabel_y)


%% Prepare for plotting
% Extract parameters
ChiMax = length(SingVals(5,:));
Biparts = length(SingVals(:,1));
width = 100;

% Transform Data
dummy1 = SingVals;
dummy2= chi_crit;
SingVals=zeros(Biparts*width,ChiMax);
chi_crit = zeros(Biparts*width,1);
for n = 1:Biparts
    for m = 1:width
        SingVals((n-1)*width+m,:) = dummy1(n,:);
        chi_crit((n-1)*width+m) = dummy2(n);
        %(n-1)*width+m
    end
end

%Define axes
x=log2(1:ChiMax);
y=1:Biparts*width;


%% Plot
contourf(x,y,log10(SingVals),contourLevels,'edgecolor','none');hold on;
plot(log2(chi_crit), y,'color','black','linewidth',linWid,'linestyle','--');
%set(h, 'EdgeColor', 'none');
%
set(gca, 'ydir', 'reverse' )
set(gca, 'XTick', xmajticks);
set(gca, 'YTick', (1:width:width*Biparts)+ceil(width/2) );
set(gca,'XAxisLocation','top','YAxisLocation','left');
ax=gca;
ax.LineWidth=linWid;set(gca,'layer','top');
ax.FontSize=fontsz;

if(showaxes_x ==1)
    xticklabels(xmajticklabels);
else
    set(gca,'xticklabel',[]);
end
if(showaxes_y==1)
    yticklabels(yTickLabs);
else
    set(gca,'yticklabel',[]);
end
if(showlabel_x ==1)
    xlabel('$\alpha$','interpreter','latex' ,'FontSize',fontsz);
end
if(showlabel_y==1)
    %ylabel('$\Delta x_{\mathrm{rough}}/\Delta x_{\mathrm{fine}}$','interpreter','latex');
    ylabel('$n$-th Bipartition','interpreter','latex');
end

set(gca,'TickLabelInterpreter', 'latex');
%
% %Minor Ticks x
axesH = gca;
%set(axesH, 'XMinorTick','on', 'XMinorGrid','on')
%axesH.XAxis.MinorTickValues = xminticks;
axesH.XGrid = 'on';
%
%Grid of Ticks y
%axesH = gca;
%set(axesH, 'YMinorTick','on', 'YMinorGrid','on')
%axesH.YAxis.MinorTickValues = (1:width:width*Biparts+ceil(width/2));


end