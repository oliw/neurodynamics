function [C ] = Question1( p )
%QUESTION1 Summary of this function goes here
%   Detailed explanation goes here

% Build small-world modular networks of Izhikevich neurons

modules = 8;
excitoryNeuronsPerModule = 100;
inhibatoryNeurons = 200;
randomWiringsPerModule = 1000;

% Initialise connectivity matrix
C = zeros(modules*excitoryNeuronsPerModule+inhibatoryNeurons);

% Each module contains 1000 randomly assigned directed inner connections (before
% rewiring)
% TODO Some neurons may be completely cut off, is this ok?
for module=1:modules
    firstNeuron = (module-1)*excitoryNeuronsPerModule+1;
    lastNeuron = module*excitoryNeuronsPerModule;
    C(firstNeuron:lastNeuron,firstNeuron:lastNeuron) = randomWiring(randomWiringsPerModule, excitoryNeuronsPerModule);
end

% Rewiring with probability p


% Each inhibatoryNeuron projects to every neuron in the whole network.
% TODO At the moment an inhibatory neuron inhibits itself
C(modules*excitoryNeuronsPerModule+1:size(C,2),:) = 1;

% Each inhibatory neuron recieves input from FOUR excitatory neurons
% which must all come from the same module
for module=1:modules
    randomList = randperm(excitoryNeuronsPerModule);
    randomList = randomList(1:4);
    randomList = randomList + (module-1)*excitoryNeuronsPerModule;
    C(randomList,801:1000) = 1;
end

C(1:800,1:800) = rewire(C(1:800,1:800), p, modules, excitoryNeuronsPerModule);

% Each neuron has a background chance of firing with lambda = 0.01

% Each connection has a weight, scaling facotr, projecton pattern(?),
% conduction delay

% Build Weight, Scaling Factor and ConductionDelay
Weight = zeros(size(C));
ScalingFactor = zeros(size(C));
ConductionDelay = zeros(size(C));
for i=1:800
   for j=1:800
       if C(i,j)==1
          Weight(i,j) =1;
          ScalingFactor(i,j) = 17;
          ConductionDelay(i,j) = randi([0,20]);
       end
   end
end
for i=1:800
   for j=801:1000
       if C(i,j)==1
          Weight(i,j) = rand;
          ScalingFactor(i,j) = 50;
          ConductionDelay(i,j) = 1;
       end
   end
end    
for i=801:1000
   for j=1:800
      if C(i,j)==1
          Weight(i,j) =-rand;
          ScalingFactor(i,j) = 2;
          ConductionDelay(i,j) = 1;
      end
   end
end
for i=801:1000
   for j=801:1000
      if C(i,j)==1
          Weight(i,j) =-rand;
          ScalingFactor(i,j) = 1;
          ConductionDelay(i,j) = 1;
      end
   end
end

% NEXT TODO: Replicate the firing




end

function module = randomWiring(numOfWirings, size)
module = zeros(size);
wiringsMade = 0;
while wiringsMade < numOfWirings
    randomI = randi([1, size]);
    randomJ = randi([1, size]);
    if randomI ~= randomJ
        if module(randomJ, randomI) == 0 && module(randomI, randomJ) == 0
           module(randomI, randomJ) = 1;
           wiringsMade = wiringsMade + 1;
        end
    end
end
end

% Pre: oldmatrix contains just the excitatory modules 
function matrix = rewire(oldmatrix, p, modules, excitoryNeuronsPerModule)

matrix = zeros(size(oldmatrix));
for module=1:modules
    firstNeuron = (module-1)*excitoryNeuronsPerModule+1;
    lastNeuron = module*excitoryNeuronsPerModule;
    for i=firstNeuron:lastNeuron
       for j=firstNeuron:lastNeuron
           if oldmatrix(i,j)==1
              if rand <= p
                 % Rewire edge
                 edges = 1:800;
                 edges(firstNeuron:lastNeuron) = [];
                 while true
                     newDest = edges(randi([1, size(edges,2)]));
                     if matrix(newDest, i) == 0
                        matrix(i,newDest) = 1;
                        break;
                     end
                 end
              else
                 % Keep unchanged
                 matrix(i,j) = 1;
              end
           end
       end
    end
end

end
