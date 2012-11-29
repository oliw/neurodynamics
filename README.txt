QUESTION 1a:
To generate the matrix connectiviy plot run Question1.m
Then type "spy([layer{1}.S{1} layer{1}.S{2}; layer{2}.S{1} layer{2}.S{2}])"

Plots - 
results/1a/connectivity00.fig
results/1a/connectivity01.fig
results/1a/connectivity02.fig
results/1a/connectivity03.fig
results/1a/connectivity04.fig
results/1a/connectivity05.fig


QUESTION 1b:
To generate the raster plot run Question1.m
Then type "scatter(layer{1}.firings(:,1), layer{1}.firings(:,2))"

Plots - 
results/1b/firings00.fig
results/1b/firings01.fig
results/1b/firings02.fig
results/1b/firings03.fig
results/1b/firings04.fig
results/1b/firings05.fig


QUESTION 1c:
To generate the mean firings rate plot run Question1c.m

Plots - 
results/1c/meanrate00.fig
results/1c/meanrate01.fig
results/1c/meanrate02.fig
results/1c/meanrate03.fig
results/1c/meanrate04.fig
results/1c/meanrate05.fig

Matlab Data
results/1c/meansfirings60seconds.mat/

QUESTION 2:
Question2.m runs multiple simulations over 60 seconds each and ouputs neural complexity.
complexity.mat contains our simulation results over 161 runs
The graph can be plotted by typing "scatter(c(:,2), c(:,1))"  
 
Plots - 
results/2/neuralcomplexity.fig/

Matlab figures -
results/2/complexity.mat
