I = [1 -1 1 1 -1 1];
x = [1/2 1 1/2];

% viterbi(I, x, 'unger');

% First, get channel response
y = conv(I, x);
L = (length(x) - 1) / 2;
y = y(L+1:end-L);

Mk = nan(2, length(I));
survivor = nan(2, length(I));

states = dec2bin(0:2^L-1);
states = flipud(states);
% t_states = dec2bin(0:2^L-1);
% for i=1:2^L
%     for j=1:L
%         states(i,j) = str2double(t_states(i,j));
%     end
% end
% states(states == 0) = -1;
alphabet = [-1 1];
% Metric: Mk-1(Ik-1) + Re{ I*k[2yk - ikxo -2sumiL(xiIk-i)]}
% k: index on time
% xo: center tap of channel autocorr
% xi: causal portion of channel autocorr
xo = x(L+1);
xi = x(L+2:end);
M = nan(length(states), length(I));
S = nan(length(states), length(I));
Ikm1 = '';
% iterate over time (left->right)
for k=1:length(I)
    % iterate over states(top->bottom)
    for i=1:length(states)
        Ik = bin2dec(states(i,:));
        if Ik == 0, Ikbpsk = -1; else Ikbpsk = Ik; end;
        % iterate over prev. states
        Mb = zeros(2,1);
        if k == 1
            M(i,k) = real(Ikbpsk' * (2*y(k) - Ikbpsk*xo));
            fprintf('k: %i Ikm1: %2i Ik: %2i yk: %.1f Mb: %2i\n', k, 0, Ikbpsk, y(k), M(i,k));
        else
            Ikm1(1) = states(floor(Ik/2)+1,:); % get prev state
            Ikm1bpsk = str2double(Ikm1(1,:)); % conv to double
            Ikm1bpsk(Ikm1bpsk==0) = -1; % change 0 to -1 (for BPSK)
            Mb(1) = M(1,k-1) + real(Ikbpsk' * (2*y(k) - Ikbpsk*xo - 2*sum(xi.*Ikm1bpsk)));
            fprintf('k: %i Ikm1: %2i Ik: %2i yk: %.1f Mb: %2i\n', k, Ikm1bpsk, Ikbpsk, y(k), Mb(1));
            
            Ikm1(2,:) = (states(floor(Ik/2)+length(states)/2+1,:));
            Ikm1bpsk = str2double(Ikm1(2,:));
            Ikm1bpsk(Ikm1bpsk==0) = -1;
            Mb(2) = M(2,k-1) + real(Ikbpsk' * (2*y(k) - Ikbpsk*xo - 2*sum(xi.*Ikm1bpsk)));
            fprintf('k: %i Ikm1: %2i Ik: %2i yk: %.1f Mb: %2i\n', k, Ikm1bpsk, Ikbpsk, y(k), Mb(2));
            
            [maxMb, maxIdx] = max(Mb);
            M(i,k) = maxMb;
            S(i,k) = str2double(Ikm1(maxIdx,:));
            
            
        end
    end
end
