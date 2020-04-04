clear;
clc;
SNR_dB = (0 : 5 : 40);  
SNR_linear = 10 .^ (SNR_dB / 10);
N_iter = 200; 
n = 256; % number of beams (transmit antennas)
K = 16; % number of users
Q = K * 6;
L = 3; % number of paths per user
V = 8; % Retained number of elements for each component
lamada = 1; % wavelength
N = K; % number of retained RF chains
d = lamada / 2;
NMSE = zeros(1, length(SNR_dB));
for i_snr = 1 : length(SNR_dB) 
	disp(['the SNR is :', num2str(i_snr)]);
	SNR = SNR_linear(i_snr);
	sigma2 = 1/SNR_linear(i_snr);
	temp = 0; temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0; temp5 = 0; temp6 = 0;
	error0 = 0; error1 = 0; error2 = 0; error3 = 0; error4 = 0; error5 = 0; error6 = 0; 
	for iter = 1 : N_iter
		H = beamspace_channel(n, K, L); % generate the signalspace channel
		U = zeros(n, n); 
		deta = 1 / n;
		for i = -(n-1) / 2 : 1 : (n-1) / 2
			U(:,i+(n+1)/2) = sqrt(1/n)*exp(1i*[0:n-1]*2*pi*deta*i).';
		end
        H_beam = U.'*H; % beamspace channel
		for j = 1 : K
			h = H_beam(:, j);
			Phi = rand(6 * K, n); Phi = Phi > 0.5; Phi = 2 * Phi - 1; 
			noise = sqrt(sigma2) * (randn(6 * K, 1) + 1i * randn(6 * K, 1)) / sqrt(2);
			x = Phi * h + noise;

			%%%SMD%%%
			%%%OMP%%%
			h_hat0 = OMP_new(x, Phi, 16, 16);
			error0 = error0 + norm(h - h_hat0, 2)^2 / norm(h, 2)^2;
                        %%%OMP end%%%

			h_hat1 = second_LS(x, Phi, h_hat0, h, L);
			error1 = error1 + norm(h - h_hat1, 2)^2 / norm(h, 2)^2;
			%%%SMD end%%%

			%%%SD%%%
			[h_hat2, support2] = SD(x, Phi, L, V);
			error2 = error2 + norm(h - h_hat2, 2)^2 / norm(h, 2)^2;
			%%%SD end%%%
		end
		F = H_beam * inv(H_beam' * H_beam);
		beta = sqrt(K / trace(F * F'));
		H_eq = H_beam' * F;
		for k=1 : K
			sum_inf = sum(abs(H_eq(:, k)).^ 2) - abs(H_eq(k, k))^2;
			temp = temp + log2(1 + abs(H_eq(k, k))^2 / (sum_inf + K/(SNR * beta^2)));
		end
	end
    NMSE0(i_snr) = error0/K/N_iter;
    NMSE1(i_snr) = error1/K/N_iter;
    NMSE2(i_snr) = error2/K/N_iter;  

end
 
figure(1)
semilogy(SNR_dB, NMSE0, '-.g', 'Marker', 's');
hold on 
semilogy(SNR_dB, NMSE1, '--b', 'Marker', 'o');
hold on 
semilogy(SNR_dB, NMSE2, '-r', 'Marker', 'd');
hold on
grid on
xlabel('SNR (dB)');
ylabel('NMSE (dB)');
axis([0, 40, 1e-2,1e1]);


leg = legend('Conventional OMP-based channel estimation',...
    'Conventional SMD-based channel estimation',...
    'Proposed SD-based channel estimation');
set(leg, 'Location', 'NorthEast') 
