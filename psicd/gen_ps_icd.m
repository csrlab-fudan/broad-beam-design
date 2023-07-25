function v = gen_ps_icd(pattern, N, iterNum)
K = length(pattern);
Omega = -1+(2*(1:K)-1)/K;
A = exp(1j*pi*(0:N-1).'*Omega);
g = pattern .* gen_phase(A, pattern, iterNum);
v = 1/K*A*g;
v = v/norm(v);
end

function f = gen_phase(A, pattern, iterNum)
[~, K] = size(A);
f = exp(1j*2*pi*rand(K, 1));
g= f.*pattern;
AHA = A'*A;
R = [real(AHA), -imag(AHA); imag(AHA), real(AHA)];
for i = 1:iterNum
    k = mod(i-1, K)+1;
    f(k) = phase_opt(R, g, k);
    g(k) = abs(g(k))*f(k);
end
end

function phase = phase_opt(R, g, k)
K = length(g);
t = [real(g); imag(g)];
p = R(k, :)*t - R(k,[k,k+K])*t([k,k+K]);
q = R(k+K, :)*t - R(k+K,[k,k+K])*t([k,k+K]);
dp = t.'*R(:, k) - t([k,k+K]).'*R([k,k+K], k);
dq = t.'*R(:, k+K) - t([k,k+K]).'*R([k,k+K], k+K);
alpha = p+dp;
beta = q+dq;
phase = exp(1j*angle(alpha+1j*beta));
end