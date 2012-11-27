function [ complexity ] = Question2RunSimComplexity(numberOfRuns)
%RUNSIM2 Summary of this function goes here
%   Detailed explanation goes here

P_MIN = 0.05;   % instead choose random probabilities
P_MAX = 0.5;
STEP = (P_MAX-P_MIN)/numberOfRuns;

complexity = zeros(numberOfRuns+1, 1);

for r=1:numberOfRuns+1
   
   fprintf('\nStarting run %d...\n', r); 
    
   p = P_MIN+(r-1)*STEP;
    
   [ means, layer ] = Question2RunMinute( p );
   
   complexity(r) = Question2CalcComplexity(means);
   
   complexity(r)
    
end

end

