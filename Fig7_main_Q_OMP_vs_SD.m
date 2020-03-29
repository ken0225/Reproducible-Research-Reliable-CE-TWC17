clear;
clc;
SNR_dB=10;  
SNR_linear=10 .^ (SNR_dB/10);
N_iter=200; 
n = 256; % number of beams (transmit antennas)
K = 16; % number of users
Q = [64, 80, 100, 120, 140, 160, 180, 200, 220, 240, 256]; % number of pilots, Q=M*K, M is number of blocks
L = 3; % number of paths per user
V = 8; % Retained number of elements for each component
lamada = 1; % wavelength
N = K; % number of retained RF chains
d = lamada / 2;
NMSE = zeros(1, length(SNR_dB));
% sigma2=1/SNR_linear(5);
for i_Q=1 : length(Q) 
	disp(i_Q);
	SNR=SNR_linear;
	sigma2=1 / SNR_linear;
	temp = 0; temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0; temp5 = 0; temp6 = 0;
	error0 = 0; error1 = 0; error2 = 0; error3 = 0; error4 = 0; error5 = 0; error6 = 0; 
	for iter = 1 : N_iter
		H = beamspace_channel(n, K, L); % generate the signalspace channel
		U = zeros(n, n); 
		deta = 1 / n;
		for i = -(n-1)/2 : 1 : (n-1)/2
			U( : , i+(n+1)/2 ) = sqrt(1/n) * exp(1i * [0 : n-1] * 2 * pi * deta * i).';
		end
		H_beam = U.' * H; % beamspace channel
		for j = 1 : K
			h = H_beam( : , j );
			Phi = rand(Q(i_Q) , n); Phi = Phi>0.5; Phi=2*Phi-1; 
			noise = sqrt(sigma2) * (randn(Q(i_Q), 1) + 1i * randn(Q(i_Q), 1)) / sqrt(2);
			z = Phi * h + noise;
			h_hat1 = OMP_new(z, Phi, 16, 16);
			error1 = error1 + norm(h - h_hat1, 2)^2 / norm(h, 2)^2;
			[h_hat2, support2] = SD(z, Phi, L, V);
			error2 = error2 + norm(h - h_hat2, 2)^2 / norm(h, 2)^2;
		end
        % ZF precoding
        % Full digital
		F = H_beam*inv(H_beam'*H_beam);
		beta = sqrt(K/trace( F * F' ));
		H_eq=H_beam' * F;
		for k=1 : K
			sum_inf=sum(abs(H_eq( : , k )).^2)-abs(H_eq( k , k ))^2;
			temp=temp+log2(1+abs(H_eq( k , k ))^2/(sum_inf+K/(SNR*beta^2)));
		end
	end
	NMSE1(i_Q) = error1/K/N_iter;
	NMSE2(i_Q) = error2/K/N_iter;  
   %NMSE3(i_snr) = error3/K/N_iter;
end
 
figure(1)
semilogy(Q, NMSE1, '--b', 'Marker', 'o');
hold on 
semilogy(Q, NMSE2, '-r', 'Marker', 'd');
hold on
grid on
xlabel('Total number of instants Q');
ylabel('NMSE (dB)');
axis([64, 256, 1e-2, 1e0]);

leg = legend('Conventional OMP-based channel estimation',...
	'Proposed SD-based channel estimation');
set(leg, 'Location', 'NorthEast') 
