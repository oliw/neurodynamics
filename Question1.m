function [layer ] = Question1( p )
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
layer{2}.columns = INHIBITORY_NEURONS;
layer{2}.a = 0.02*ones(INHIBITORY_NEURONS,1); % Set every a of each neuron to 0.02
layer{2}.b = 0.25*ones(INHIBITORY_NEURONS,1); % b is the same for everyone
layer{2}.c = -65+15*r.^2; %  c's are all slightly different thanks to r
layer{2}.d = 2-6*r.^2; % d's are all different thanks to r

% Clear connectivity matrices
L = length(layer);
for i=1:L
   for j=1:L
      layer{i}.S{j} = [];
      layer{i}.factor{j} = [];
      layer{i}.delay{j} = [];
   end
end

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

% Each inhibatoryNeuron projects to every neuron in the whole network.
% Connections from inhibitory neurons all have a weight between -1 and 0.
layer{2}.S{2} = -rand(INHIBITORY_NEURONS);
layer{1}.S{2} = -rand(EXCITATORY_NEURONS,INHIBITORY_NEURONS);

for module=1:MODULES
    randomList = randperm(EXCITATORY_NEURONS_PER_MODULE);
    chosenExNeurons = randomList(1:4) + (module-1)*EXCITATORY_NEURONS_PER_MODULE;
    % Any Excitatory to Inhibatory connections have a random weight between 0
    % and 1
    layer{2}.S{1}(chosenExNeurons,:) = rand;
end

% Set scaling factors
layer{1}.factor{1} = 17;
layer{1}.factor{2} = 2;
layer{2}.factor{1} = 50;
layer{2}.factor{2} = 1;

% Set conduction delay
layer{1}.delay{1} = 20.*rand(EXCITATORY_NEURONS,EXCITATORY_NEURONS);
layer{1}.delay{2} = ones(EXCITATORY_NEURONS,INHIBITORY_NEURONS);
layer{2}.delay{1} = ones(INHIBITORY_NEURONS,EXCITATORY_NEURONS);
layer{2}.delay{2} = ones(INHIBITORY_NEURONS, INHIBITORY_NEURONS);


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

function layer = IzNeuronUpdate(layer,t,Dmax)
% Updates membrane potential v and reset rate u for neurons in layer 1
% using Izhikevich's neuron model and the Euler method. Dmax is the maximum
% conduction delay
i = 1;
dt = 0.2; % Euler method step size
% Calculate current from incoming spikes
for j=1:length(layer)
   S = layer{i}.S{j};
   if ~isempty(S)
      firings = layer{j}.firings;
      if ~isempty(firings)
         % Find incoming spikes (taking account of propagation delays)
         delay = layer{i}.delay{j};
         F = layer{i}.factor{j};
         % Sum current from incoming spikes
         k = size(firings,1);
         while (k>0 && firings(k,1)>t-(Dmax+1))
            spikes = (delay(:,firings(k,2))==t-firings(k,1));
            if ~isempty(layer{i}.I(spikes))
               layer{i}.I(spikes) = layer{i}.I(spikes)+S(spikes,firings(k,2))*F;
            end
            k = k-1;
         end;
         % Don't let I go below zero (shunting inhibition)
         % layer{i}.I = layer{i}.I.*(layer{i}.I > 0);
      end
   end
end

% Update v and u using Izhikevich's model in increments of dt
for k=1:1/dt
   v = layer{i}.v;
   u = layer{i}.u;
   layer{i}.v = v+(dt*(0.04*v.^2+5*v+140-u+layer{i}.I));
   layer{i}.u = u+(dt*(layer{i}.a.*(layer{i}.b.*v-u)));
   % Reset neurons that have spiked
   fired = find(layer{i}.v>=30); % indices of spikes
   if ~isempty(fired)
      layer{i}.firings = [layer{i}.firings ; t+0*fired, fired];
      layer{i}.v(fired) = layer{i}.c(fired);
      layer{i}.u(fired) = layer{i}.u(fired)+layer{i}.d(fired);
   end
end
end
