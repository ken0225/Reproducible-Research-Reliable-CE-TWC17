function [ H ] = beamspace_channel(n, K, L)
% n: number of transmit beams (transmit antennas)
% K: number of users
% L: number of paths, L=1 for LoS, L>1 for multipath
lamada = 1; % wavelength
d = lamada / 2; % distance between two antennas
H = zeros(n, K);
theta=zeros(L, K);
for k=1 : K
	beta=zeros(1,L); 
	% complex gain
	beta(1) = 1; % gain of the LoS
	beta(2 : L) = 10 ^ (-0.5) * (randn(1, L-1) + 1i * randn(1, L-1)) / sqrt(2);
	% NLoS
	% DoA
	theta(1, k) = pi * rand(1) - pi / 2;
	theta(2 : L, k) = pi * rand(1, L-1) - pi / 2;
	% channel for the k th user
for j = 1 : L
		H(:, k) = H(:, k) + beta(j) * array_respones(theta(j, k), n, d, lamada);
end
end
