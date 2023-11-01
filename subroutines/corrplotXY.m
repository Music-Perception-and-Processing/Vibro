function [rval, pval] = corrplot3(x,y, varargin)
% corrplot with cool line ! 
% if flag == 1 then NO print! 

b = regress(y, [ones(size(x,1),1), x]); % model
yhat = b(2)*x + b(1);  % prediction
%figure('units','normalized','outerposition',[.35 .3 .35 .55]);
hold on; % cool down
%plot(x, y, 'o', 'color', [0 .447 .7410]) % scatter
[rval, pval] = corr(x,y, 'rows', 'complete');
if pval < 0.05
    plot(x, yhat, 'color', [.75 .75 .75], 'linewidth', 4) % line
    if nargin == 2
        %text(min(x)+abs(max(x)-min(x))/10, max(y)-abs(max(y)-min(y))/10, strcat('r=', num2str(round(rval*100)/100)), 'fontsize', 20) 
        text(12, 3, strcat('r=', num2str(round(rval*100)/100)), 'fontsize', 16) 
    end
else 
    %text(min(x)+abs(max(x)-min(x))/10, max(y)-abs(max(y)-min(y))/10, strcat('r=', num2str(round(rval*100)/100), '(ns)'), 'fontsize', 20) 
    %text(min(x)+abs(max(x)-min(x))/10, max(y), strcat('r=', num2str(round(rval*100)/100), '(ns)'), 'fontsize', 20) 
    text(12, 3, strcat('r=', num2str(round(rval*100)/100), '(ns)'), 'fontsize', 16) 
end
hold off
grid on
set(gca, 'box', 'off')
grid off
end