function [ layer ] = Question1c( p )
%QUESTION1C Generates a plot of the mean firing rate 
%in each of the eight modules for a 1000ms run
%
%   INPUTS
%       p     - rewiring probability
%
%   OUTPUTS
%       layer - same object as in Tutorial 2
%             - layer{1} contains all the excitatory neurons
%             - layer{2} contains all the inhibatory neurons


[ layer ] = Question1(p, 1000);


firings = layer{1}.firings;

x = zeros(1000, 8);

% COLLECT FIRINGS IN EACH RANGE
for t=1:1000
 
    x(t, 1) = sum(firings(:, 1) == t & firings(:, 2) <= 100); 
    x(t, 2) = sum(firings(:, 1) == t & firings(:, 2) > 100 & firings(:, 2) <= 200); 
    x(t, 3) = sum(firings(:, 1) == t & firings(:, 2) > 200 & firings(:, 2) <= 300);
    x(t, 4) = sum(firings(:, 1) == t & firings(:, 2) > 300 & firings(:, 2) <= 400);
    x(t, 5) = sum(firings(:, 1) == t & firings(:, 2) > 400 & firings(:, 2) <= 500);
    x(t, 6) = sum(firings(:, 1) == t & firings(:, 2) > 500 & firings(:, 2) <= 600);
    x(t, 7) = sum(firings(:, 1) == t & firings(:, 2) > 600 & firings(:, 2) <= 700);
    x(t, 8) = sum(firings(:, 1) == t & firings(:, 2) > 700 & firings(:, 2) <= 800);
    
end

% CALCULATE MEANS
means = zeros(50, 8);

for d=1:50
    
    bIndex = (d-1)*20 + 1;
    eIndex = mod(bIndex + 49, size(x, 1));
    
    means(d, 1) = mean(x(bIndex:eIndex, 1));
    means(d, 2) = mean(x(bIndex:eIndex, 2));
    means(d, 3) = mean(x(bIndex:eIndex, 3));
    means(d, 4) = mean(x(bIndex:eIndex, 4));
    means(d, 5) = mean(x(bIndex:eIndex, 5));
    means(d, 6) = mean(x(bIndex:eIndex, 6));
    means(d, 7) = mean(x(bIndex:eIndex, 7));
    means(d, 8) = mean(x(bIndex:eIndex, 8));
    
end

% GENERATE PLOT

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',{'0','100','200','300','400','500','600','700','800','900','1000'});
%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 50]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to plot
plot(1:50,means,'Parent',axes1);

% Create xlabel
xlabel('Time (ms)');

% Create ylabel
ylabel('Mean firing rate');

% Create title
title(sprintf('p = %.2f', p));

end

