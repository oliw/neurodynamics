function [ output ] = Question2( time_series )
%QUESTION2 Calculates neural complexitity, for a set of time series (nvar x
%nobs data matrix)
% There is a variable for each module

% Apply differencing twice to the set time series
time_series = aks_diff(time_series);
time_series = aks_diff(time_series);

output = complexity(time_series);

end

function output = complexity(S)
i = 0;
for n=1:size(S,2)
   i = i + (mutualInformation(n, S) - integration(S));
end
end

function output = integration(S)
i = 0;
for n=1:size(S,2)
    i = i + (entropy(S(:,n)) - entropy(S)); 
end
end

function output = mutualInformation(X, S)
% Separate X from S
Xs = S(:,X);
Ss = S;
Ss(:,X) = [];
output = entropy(Xs) + entropy(Ss) + entropy(S);
end

function output = entropy(variables)
% variables has M columns (each column a variable)
% variables has N rows for each observation of each variable
n = size(variables,2);
output = 0.5*log((2*pi*e)^n*abs(cov(variables)));
end