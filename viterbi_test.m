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
        % iterate over prev. states
        Mb = zeros(2,1);
        if k == 1
            Mb(1) = real(Ik' * (2*y(k) - Ik*xo));
        else
            Ikm1(1) = states(floor(Ik/2)+1,:);
            Mb(1) = real(Ik' * (2*y(k) - Ik*xo - 2*sum(xi.*str2double(Ikm1(1,:)))));
            Ikm1(2,:) = (states(floor(Ik/2)+length(states)/2+1,:));
            Mb(2) = real(Ik' * (2*y(k) - Ik*xo - 2*sum(xi.*str2double(Ikm1(2,:)))));
            [maxMb, maxIdx] = max(Mb);
            M(i,k) = maxMb;
            S(i,k) = str2double(Ikm1(maxIdx,:));
        end
    end
end
