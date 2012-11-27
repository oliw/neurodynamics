function [ output ] = Question2CalcComplexity( time_series )
%QUESTION2 Calculates neural complexitity, for a set of time series (nvar x
%   nobs data matrix)
%   There is a variable for each module

% Apply differencing twice to the set time series
time_series = aks_diff(time_series');
time_series = aks_diff(time_series);

output = complexity(time_series');

end

function output = complexity(S)
output = 0;

for i=1:size(S,2)
   output = output + (mutualInformation(i, S));
end

output = output - integration(S);

end

function out_i = integration(S)
out_i = 0;
for i=1:size(S,2)
    out_i = out_i + entropy(S(:,i)); 
end

out_i = out_i - entropy(S);

end

function out_mi = mutualInformation(X, S)
Xs = S(:,X);    % H(X)
Ss = S;         % H(S)
Ss(:,X) = [];   % H(S-{X})

out_mi = entropy(Xs) + entropy(Ss) - entropy(S);
end

function out_e = entropy(variables)
% variables has M columns (each column a variable)
% variables has N rows for each observation of each variable
N = size(variables,2);
out_e = 0.5*log((2*pi*exp(1))^N*det(cov(variables)));

end