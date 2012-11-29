function [ complexity ] = Question2RunSimComplexity(numberOfRuns)
%RUNSIM2 Summary of this function goes here
%   Detailed explanation goes here

complexity = zeros(numberOfRuns+1, 2);

for r=1:numberOfRuns+1
   
    
   p = 0.05 + (rand(1) * 0.45);
   fprintf('\nStarting run %d with probability %.2f...\n',r , p); 
   
   [ means, layer ] = Question2RunMinute( p );
   
   complexity(r, 1) = Question2CalcComplexity(means);
   complexity(r, 2) = p;
   
   complexity(r, :)
   
   name = sprintf('run%d.mat', r);
   save(name, 'complexity');
    
end

end

