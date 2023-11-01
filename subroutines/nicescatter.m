function [] = nicescatter(x,y, symbsty, col)

x = squeeze(x);
y = squeeze(y); 

plot(x, y, ...
    'LineStyle', 'none', 'Marker', symbsty, ...
           'MarkerFaceColor', col,'MarkerEdgeColor', [.5 .5 .5], ...
           'MarkerSize', 8); 
end
