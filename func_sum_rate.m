function  SR = func_sum_rate(W,H,s2)

%[K,~] = size(H);
%SINR = zeros(K,1);
%alpha = zeros(K,1);

%for k = 1 : K
%    for j = 1 : K
%        alpha(k) = alpha(k) + abs(H(k,:) * W(:,j)) ^ 2;
%    end
%    SINR(k) = abs(H(k,:) * W(:,k)) ^ 2 / (alpha(k) - abs(H(k,:) * W(:,k)) ^ 2 + s2);
%end
%SR = real(sum(log2(1 + SINR)));



alpha = sum(abs(H * W).^2,2);
pow = abs(diag(H * W)).^2;
SINR = pow ./ (alpha - pow + s2);
SR = real(sum(log2(1 + SINR)));

end