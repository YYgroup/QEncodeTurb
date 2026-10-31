%
%
rng(0);
clf;

subplot_tight(2,2,[1,2],[0.02 0.02],'Position',[0.1,0.7,0.3,0.2])
plot(randn(10,3)); 
subplot_tight(2,2,3,[0.02 0.02])
plot(randn(10,3)); 
subplot_tight(2,2,4,[0.02 0.02])
plot(randn(10,3)); 

