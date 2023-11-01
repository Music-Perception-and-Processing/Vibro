function [rval, pval, cil, ciu] = corrplot2(x,y, varargin)
% corrplot with cool line ! 
% if flag == 1 then NO print! 

b = regress(y, [ones(size(x,1),1), x]); % model
yhat = b(2)*x + b(1);  % prediction
%figure('units','normalized','outerposition',[.35 .3 .35 .55]);
plot(x, yhat, 'color', [.75 .75 .75], 'linewidth', 4) % line
hold on; % cool down
plot(x, y, 'o', 'color', [.25 .25 1]) % scatter
[rval, pval, cil, ciu] = corrcoef(x,y);
rval = rval(1,2);
pval = pval(1,2);
cil = cil(1,2); 
ciu = ciu(1,2);
if nargin == 2
    text(min(x)+abs(max(x)-min(x))/10, max(y)-abs(max(y)-min(y))/10, strcat('r=', num2str(round(rval*100)/100)), 'fontsize', 20) 
end
hold off
grid on
end