function [h_out] = second_LS(x, Phi, h_hat, h, L)
n = length(h_hat);
x_temp = x;
for l = 1 : L
	[~, order] = max(abs(h_hat));
	select(l, :) = order - (9-1)/2 : 1 :order + (9-1)/2;
	for i = 1 : length(select(l, :))
		if select(l, i) > n
			select(l, i) = select(l, i) - n;
		elseif select(l, i) < 1
			select(l, i) = select(l, i) + n;
		end
	end
	Phi2 = Phi(:, select(l, :));
	h_hat2 = inv(Phi2' * Phi2) * Phi2' * x_temp;
	temp = zeros(n, 1);
	temp(select(l, :)) = h_hat2;
	if l >= 1
		x_temp = x_temp - Phi * temp;
		h_hat = OMP_new(x_temp, Phi, 9, 9);
	end
end
select_final = unique(reshape(select,L * 9,1));
Phi_final = Phi(:, select_final);
est =  inv(Phi_final' * Phi_final) * Phi_final' * x;
h_out = zeros(n, 1);
h_out(select_final) = est;
end
