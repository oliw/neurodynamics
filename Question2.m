function [ complexity ] = Question2(numberOfRuns)
%RUNSIM2 Top-level function for Question 2

complexity = zeros(numberOfRuns+1, 2);

for r=1:numberOfRuns+1
    
   p = 0.05 + (rand(1) * 0.45);
   fprintf('\nStarting run %d with probability %.2f...\n',r , p); 
   
   [ means, layer ] = Question2RunMinute( p );
   
   complexity(r, 1) = Question2CalcComplexity(means);
   complexity(r, 2) = p;
   
   %complexity(r, :)    % print status for long simulations
   
   %Save intermediate results for long simulations
   %name = sprintf('run%d.mat', r);
   %save(name, 'complexity');
    
end

end

