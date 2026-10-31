function [Smax,bond,errSmax_vec,errT_vec] = TNanalyseFlow(SingMats, entanEntros)
%Performs an entanglement-spectrum analysis of the inputed MPS. 
%First the maxmium entanglement entropy is identified along with the 
%associated bond number. Then the total truncation error as a function of
%chi, along with the truncation error specific to the Singular Values
%Matrix with the greatest entanglement entropy also as a function of chi,
%are calculated and outputed.

%% Identify maximum Entanglement Entropy along with associated bond
[Smax,bond]= max(entanEntros(1,:));

%% Calculate Truncation errors
numBonds = length(SingMats);
D_max = size(SingMats{round(numBonds/2)}); D_max = D_max(1);
Sdiags = zeros(D_max,numBonds);
for i = 1:numBonds
    D_i = size(SingMats{i}); D_i = D_i(1);
    Sdiags(1:D_i,i) = diag(SingMats{i});
end
Sdiags = Sdiags.^2./sum(Sdiags.^2);

errSmax_vec = zeros(D_max,1); 
errT_vec = errSmax_vec;
errSmax = 0;
errT = 0;
for chi = D_max:(-1):1
    errSmax = errSmax + Sdiags(chi,bond);
    errT = errT + sum(Sdiags(chi,:));
    %
    errSmax_vec(chi) = errSmax;
    errT_vec(chi) = errT;
end

% %% Test which distributions the errors follow more closely
% n = 1:length(errT_vec);
% 
% expFit = fit(n',errT_vec,'exp1');
% plFit = fit(n',errT_vec,'power1');
% 
% res_expFit = errT_vec-expFit(n);
% res_plFit =  errT_vec-plFit(n);
% 
% [~,expFit_adp] = adtest(res_expFit);
% [~, plFit_adp] = adtest(res_plFit);
% 
% expFit_lastR = res_expFit(end);
% plFit_lastR = res_plFit(end);
% 
% %expFit_adp = 0; expFit_lastR = 0; plFit_adp = 0; plFit_lastR = 0;

end