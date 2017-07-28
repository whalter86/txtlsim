%% TXTL Tutorial
% txtl_tutorial.m - basic usage of the TXTL modeling toolbox
% Vipul Singhal, 28 July 2017
%
% This file contains a simple tutorial of the TXTL modeling toolbox. You
% will learn about setting up a negative autoregulation circuit, simulating it, 
% plotting the results, creating variations of the circuit, and
% understanding the object structure of the models. 

%% Initializing the toolbox
% Use this command to add the subdirectories needed to your matlab path. To
% be run each time you begin a new TXTL toolbox session. 
txtl_init;

%% Negative Autoregulation - A simple example
% Here we demonstrate the setup of a genetic circuit where a transcription
% factor represses its own expression. 
% 
% Set up the standard TXTL tubes
% These load up the RNAP, Ribosome and degradation enzyme concentrations
% ``E30VNPRL'' refers to a configuration file 
tube1 = txtl_extract('E30VNPRL');
tube2 = txtl_buffer('E30VNPRL');

% Now set up a tube that will contain our DNA
tube3 = txtl_newtube('gene_expression');

% Define the DNA strands (defines TX-TL species + reactions)
dna_deGFP = txtl_add_dna(tube3, ...
  'p70(50)', 'rbs(20)', 'deGFP(1000)', ...	% promoter, rbs, gene
   30, ...					% concentration (nM)
  'plasmid');					% type

% Mix the contents of the individual tubes
Mobj = txtl_combine([tube1, tube2, tube3]);

%
% Run a simulaton
%
% At this point, the entire experiment is set up and loaded into 'Mobj'.
% So now we just use standard Simbiology and MATLAB commands to run
% and plot our results!
%



cs = getconfigset(Mobj);
set(cs.RuntimeOptions, 'StatesToLog', 'all');
tic
[simData] = txtl_runsim(Mobj,14*60*60);
toc
t_ode = simData.Time;
x_ode = simData.Data;



%% plot the result

txtl_plot(simData,Mobj);



% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:
