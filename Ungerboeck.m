%% Ungerboeck Observation Model
%
% Ik -> X(Z) -> + ->yk
%               ^
%               |
%               vk
% X(Z): x_{-L} ... x_0 ... x_L
% X(Z) = F(Z) F*(1/z*) 
sim_nSymbols = 10;

clear('-regexp', '^(?!sim_.*$).');

h = load('channel4ecen776.mat');
Ik = randi([0 1], sim_nSymbols, 1);
