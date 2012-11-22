function [layer] = Question1( p )
%QUESTION1 Summary of this function goes here
%   Detailed explanation goes here
% CONSTANTS

% Build small-world modular networks of Izhikevich neurons

MODULES = 8;
EXCITATORY_NEURONS_PER_MODULE = 100;
RANDOM_WIRING_PER_MODULE = 1000;

EXCITATORY_NEURONS = MODULES*EXCITATORY_NEURONS_PER_MODULE;
INHIBITORY_NEURONS = 200;
NEURONS = EXCITATORY_NEURONS + INHIBITORY_NEURONS;

% Build two layer structure

% Layer 1 (Exitatory neurons) Regular Spiking
r = rand(EXCITATORY_NEURONS,1);
layer{1}.rows = EXCITATORY_NEURONS;
layer{1}.columns = 1;
layer{1}.a = 0.02*ones(EXCITATORY_NEURONS,1); % Set every a of each neuron to 0.02
layer{1}.b = 0.2*ones(EXCITATORY_NEURONS,1); % b is the same for everyone
layer{1}.c = -65+15*r.^2; %  c's are all slightly different thanks to r
layer{1}.d = 8-6*r.^2; % d's are all different thanks to r

% Layer 2 (Inhibatory neurons) Fast Spiking
r = rand(INHIBITORY_NEURONS,1);
layer{2}.rows = INHIBITORY_NEURONS;
layer{2}.columns = 1;
layer{2}.a = 0.02+0.08*r;
layer{2}.b = 0.25-0.05*r; 
layer{2}.c = -65*ones(INHIBITORY_NEURONS,1); 
layer{2}.d = 2*ones(INHIBITORY_NEURONS,1); 

% Clear connectivity matrices
L = length(layer);
for i=1:L
   for j=1:L
      layer{i}.S{j} = [];
      layer{i}.factor{j} = [];
      layer{i}.delay{j} = [];
   end
end

% Initialize connections from excitatory neurons
layer{1}.S{1} = zeros(EXCITATORY_NEURONS);
layer{2}.S{1} = zeros(INHIBITORY_NEURONS, EXCITATORY_NEURONS);

% Each module contains 1000 randomly assigned directed inner connections (before
% rewiring)
% TODO Some neurons may be completely cut off, is this ok?
for module=1:MODULES
    firstNeuron = (module-1)*EXCITATORY_NEURONS_PER_MODULE+1;
    lastNeuron = module*EXCITATORY_NEURONS_PER_MODULE;
    % Connections between excitatory neurons have weight 1
    layer{1}.S{1}(firstNeuron:lastNeuron,firstNeuron:lastNeuron) = randomWiring(RANDOM_WIRING_PER_MODULE, EXCITATORY_NEURONS_PER_MODULE);
end
% Rewiring with probability p
layer{1}.S{1} = rewire(layer{1}.S{1}, p, MODULES, EXCITATORY_NEURONS_PER_MODULE);

% Each inhibitory neuron projects to every neuron in the whole network.
% Connections from inhibitory neurons all have a weight between -1 and 0.
layer{2}.S{2} = -rand(INHIBITORY_NEURONS);
layer{1}.S{2} = -rand(EXCITATORY_NEURONS,INHIBITORY_NEURONS);

% Each inhibitory neuron has connections from FOUR excitory neurons
% All must come from the same module
for ineuron=1:INHIBITORY_NEURONS
    module = randi([1,MODULES]);
    randomList = randperm(EXCITATORY_NEURONS_PER_MODULE);
    chosenExNeurons = randomList(1:4) + (module-1)*EXCITATORY_NEURONS_PER_MODULE;
    layer{2}.S{1}(ineuron, chosenExNeurons) = rand;
end

% Set scaling factors
layer{1}.factor{1} = 17;
layer{1}.factor{2} = 2;
layer{2}.factor{1} = 50;
layer{2}.factor{2} = 1;

% Set conduction delay
layer{1}.delay{1} = 20*rand(EXCITATORY_NEURONS,EXCITATORY_NEURONS);
layer{1}.delay{2} = ones(EXCITATORY_NEURONS,INHIBITORY_NEURONS);
layer{2}.delay{1} = ones(INHIBITORY_NEURONS,EXCITATORY_NEURONS);
layer{2}.delay{2} = ones(INHIBITORY_NEURONS, INHIBITORY_NEURONS);


% SIMULATE

Tmax = 1000;
Ib = 15;
Dmax = 25; % maximum propagation delay in milliseconds. The time it takes to go from one neuron to another.

% Generate a random time for each neuron to fire
% A neuron will only background fire once in a second
layer{1}.backgroundFire = randi([1, 1000], EXCITATORY_NEURONS,1);
layer{2}.backgroundFire = randi([1,1000], INHIBITORY_NEURONS, 1);

% Initialise layers with intial values of v and u and notes of which are
% firing
for lr=1:length(layer)
   layer{lr}.v = -65*ones(layer{lr}.rows,layer{lr}.columns); % Every neuron in the layer lr is given an intial voltage of -65
   layer{lr}.u = layer{lr}.b.*layer{lr}.v; % Initial value of u
   layer{lr}.firings = []; % None of the neurons are firing?
end

for t=1:Tmax
    
%    layer{1}.I = zeros(EXCITATORY_NEURONS,1);
%    layer{2}.I = zeros(INHIBITORY_NEURONS,1);
%    layer{1}.I(layer{1}.backgroundFire(:)==t) = Ib; 
%    layer{2}.I(layer{2}.backgroundFire(:)==t) = Ib; 
%    
%    if t==1
%     layer{1}.I = Ib*ones(EXCITATORY_NEURONS,1); % Layer 1 always has a constant input of current 
%     layer{2}.I = Ib*ones(INHIBITORY_NEURONS,1); % Layer 2 always has a constant input of current
%    end
   % Deliver a constant base current to layer 1
   
   
   
   
%    layer{1}.I = Ib*ones(EXCITATORY_NEURONS,1); % Layer 1 always has a constant input of current
%    layer{2}.I = Ib*ones(INHIBITORY_NEURONS,1); % Layer 2 always has a constant input of current
   
   lambda = 0.01;
   layer{1}.I = 15*(poissrnd(lambda, EXCITATORY_NEURONS, 1) > 0);
   layer{2}.I = 15*(poissrnd(lambda, INHIBITORY_NEURONS, 1) > 0);   
   
   % Update all the neurons
   for lr=1:length(layer)
      layer = IzNeuronUpdate(layer,lr,t,Dmax); % Takes the neurons, the layer index we want to update, the time, and the propagation delay. It calculates the new v and u values 
   end
   
    v1(t,1:800) = layer{1}.v;
    v2(t,1:200) = layer{2}.v;

    u1(t,1:800) = layer{1}.u;
    u2(t,1:200) = layer{2}.u;
    
end


firings1 = layer{1}.firings;
firings2 = layer{2}.firings;

% Add Dirac pulses (mainly for presentation)
if ~isempty(firings1)
   v1(sub2ind(size(v1),firings1(:,1),firings1(:,2))) = 30;
end
if ~isempty(firings2)
   v2(sub2ind(size(v2),firings2(:,1),firings2(:,2))) = 30;
end

figure(3)
clf

subplot(2,1,1)
if ~isempty(firings1)
   plot(firings1(:,1),firings1(:,2),'.')
end
% xlabel('Time (ms)')
xlim([0 Tmax])
ylabel('Neuron number')
ylim([0 800+1])
set(gca,'YDir','reverse')
title('Population 1 firings')

subplot(2,1,2)
if ~isempty(firings2)
   plot(firings2(:,1),firings2(:,2),'.')
end
xlabel('Time (ms)')
xlim([0 Tmax])
ylabel('Neuron number')
ylim([0 200+1])
set(gca,'YDir','reverse')
title('Population 2 firings')

drawnow


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
