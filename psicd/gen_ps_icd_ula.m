function v = gen_ps_icd_ula(N, alpha, u, icdIterNum, hybridIterNum, rfNum, digitNum)
K = N*4;
Omega = -1+(2*(1:K)-1)/K;
left = u-alpha;
right = u+alpha;
pattern = zeros(K, 1);
pattern(Omega>left & Omega<right) = 1;
v = gen_ps_icd(pattern, N, icdIterNum);
[v, ~] = ps_icd_optimize(v, rfNum, digitNum, hybridIterNum);
end