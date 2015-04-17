function [ I_hat ] = viterbi_BPSK( I, x, model_type)
%VITERBI 

% First, get channel response
y = conv(I, x);
L = (length(x) - 1) / 2;
y = y(L+1:end-L);

state = nan(2, length(I));
survivor = nan(2, length(I));


end

