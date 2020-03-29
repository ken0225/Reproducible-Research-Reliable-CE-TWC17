function [a, support, used_iter] = OMP_new(z, Phi, s,  remaind_max_iter)

err = 10 ^ (-3); % residual error

m = size(Phi, 2); 
n = size(z, 2); % m x n sparse matrix
a = zeros(m, n); % initial sparse matrix
v = z - Phi * a;
it = 0;
stop = 0;

	while ~ stop
		y = Phi' * v;
		[~, temp2] = sort(sum(abs(y).^2, 2),'descend');
		Omega = temp2(1);

		if it == 0
			T = Omega;
		else
			T = union(Omega, T_last);
		end
		T_last = T;
		%Signal Estimation
		b = zeros(m, n);
		T_index = T;
		b(T_index, :) = Phi(:, T_index) \ z;
		%Prune
		a = b; 
		%Sample Update
		v = z - Phi * a;
		%Iteration counter
		it = it + 1;
		%Check Halting Condition
		if (it >= max(remaind_max_iter) ||  norm(v)<=err*norm(z)) % norm(r_n)<=err  % Normally: max_iter = 6*(s+1)
			stop = 1;
		end
	end

% transmitter antenna

used_iter = it;
end
