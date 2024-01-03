function w = gen_gsc_sequence(N, alpha, u, P)
f0 = mod(u, 2);
b = N*(f0-alpha)/(2*P*alpha) + 1/2;
w = zeros(N, 1);
for n = 0:N-1
    r = mod(n, P);
    s = (n-r)/P;
    phi = 2*pi/(N)*P*alpha*(n*(s+b)-P*s*(s+1)/2);
    w(n+1) = exp(1j*phi);
end
end