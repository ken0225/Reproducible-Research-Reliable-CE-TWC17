function [h_out, support] = SD(x, Phi, L, V)
n = size(Phi, 2);
x_temp = x;
	for l = 1 : L

		y = Phi' * x_temp;
		[~, temp2] = sort(sum(abs(y) .^ 2, 2), 'descend');
		order = temp2(1);

		if mod(V, 2) == 0
			select(l, :) = order - (V)/2 : 1 : order + (V-2)/2; % V=8
		elseif mod(V,2) == 1
			select(l, :) = order - (V-1)/2 : 1 : order + (V-1)/2; % V=9
		end
		for i = 1 : length(select(l, :))
			if select(l, i) > n
				select(l, i) = select(l, i) - n;
			elseif select(l, i) < 1
				select(l, i) = select(l, i) + n;
			end
		end
		Phi2 = Phi(:, select(l, :));
		h_hat2 = inv(Phi2' * Phi2) * Phi2' * x_temp; % LS
		temp = zeros(n, 1);
		temp(select(l, :)) = h_hat2;
		if l >= 1
			x_temp = x_temp - Phi * temp;
		end
	end
select_final = unique(reshape(select, L * V, 1));
Phi_final = Phi(:, select_final);
est = inv(Phi_final' * Phi_final) * Phi_final' * x;
support = select_final;
h_out = zeros(n, 1);
h_out(select_final) = est;
end