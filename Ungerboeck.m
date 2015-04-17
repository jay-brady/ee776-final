%% Ungerboeck Observation Model
%
% Ik -> X(Z) -> + ->yk
%               ^
%               |
%               vk
% X(Z): x_{-L} ... x_0 ... x_L
% X(Z) = F(Z) F*(1/z*) 

sim_nSymbols = 10;
sim_rolloff = .5;
sim_sampPerSym = 16;
sim_trimEnergy = 100; % percent
sim_SNR = 1; %
t_tmp = 1;
another_thing = 1;
clear('-regexp', '^(?!sim_.*$).');

% get channel & pulse shape
c = load('channel4ecen776.mat');
g = rcosdesign(sim_rolloff, 1, sim_sampPerSym);

% for Ungerboekc, need X(z)
h = conv(g, c.h_m4a);
h = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 .1]; % easy channel
X_us = xcorr(h); % high sample rate autocorr

% find center tap offset before downsampleing. Need a starting point from
% 1-16 (one idexed)
[~, t_ctrIdx]  = max(X_us);
ds_begin = mod(t_ctrIdx-1, 16) + 1; % mod point before ctr, then inc by 1 
X_unpruned = X_us(ds_begin:16:end); % downsample

% prune down to given amount of energy
totalEnergy = sum(abs(X_unpruned).^2);
prunedEnergy = sim_trimEnergy/100 * totalEnergy;
t_pruned = false;
X = nan(size(X_unpruned));
[~, t_ctrIdx] = max(X_unpruned);
for i=t_ctrIdx:length(X_unpruned)
    t_offset = t_ctrIdx - i;
    X(t_ctrIdx + t_offset) = X_unpruned(t_ctrIdx + t_offset);
    X(t_ctrIdx - t_offset) = X_unpruned(t_ctrIdx - t_offset);
    t_energy = nansum(abs(X).^2);
    % go until you get over the pruned energy, then remove the excess
    if t_energy > prunedEnergy
%        X(t_ctrIdx + t_offset) = NaN; % remove to keep under 85%?
%        X(t_ctrIdx - t_offset) = NaN;
       X(isnan(X)) = [];
       break;
    end    
end

% Get F(z) from X
roots_X = roots(X);
F = roots_X(abs(roots_X) < 1);
Fconj = roots_X(abs(roots_X) > 1);
if any(abs(roots_X) == 1)
   error('root of X == 1'); 
end

% find scale factor adjustment
% Fscale =  max(X) / max(abs(poly(X)));
Fscale =  sqrt(max(X) / max(conv(poly(F), poly(Fconj))));
F = Fscale * F;
Fconj = Fscale * Fconj;

% generate symbol and noise sequences
I = randi([0 1], sim_nSymbols, 1); % random 0s and 1s
I(I==0) = -1; % turn 0s into -1s
n = randn(sim_nSymbols*2, 1) * 1/sim_SNR; % arbitrarily generate 2x noise samples to account for difference in X, F length

% pass symbols thru channel, add noise
Ifilt = conv(X, I);
nu = conv(Fconj, n);
y = Ifilt + nu(1:length(Ifilt));

figure(1);
subplot(3,1,1);
stem(I);
% set(gca, 'xlim', [0 length(y) + 1], 'ylim', [-2 2]);
subplot(3,1,2);
stem(Ifilt);
% set(gca, 'xlim', [0 length(y) + 1], 'ylim', [-2 2]);
subplot(3,1,3);
stem(y)
% set(gca, 'xlim', [0 length(y) + 1], 'ylim', [-2 2]);


clear -regexp ^t_