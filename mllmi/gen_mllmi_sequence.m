function s = gen_mllmi_sequence(N, alpha, u)
up = u+alpha;
low = u-alpha;
rho = 20;
epsilon = 1e-3;

e = zeros(N, 1);
e(1) = 1;
[w, NA] = gen_bmwss_sequence(N, alpha, u);
if NA~=N
    w = exp(1j*2*pi*rand(1, N));
end
% close all
% figure
% spectrum = abs(fft(w, 40*N)).^2;
% spectrum = fftshift(spectrum);
% plot(linspace(-1, 1, length(spectrum)), spectrum)

Wk = w'*w;
loss_old = -inf;
while true
    [V, D] = eig(Wk);
    d = diag(D);
    [~, idx] = max(d);
    s = V(:, idx);
    cvx_begin sdp
    variable delta
    variable W(N, N) complex semidefinite
    variable X1(N, N) complex semidefinite
    variable X2(N, N) complex semidefinite
    variable Z(N-1, N-1) complex semidefinite
    variable tau
    expression r(N)
    maximize(delta - rho*(norm_nuc(W) - norm(Wk, 2) - real(sum(s.*(s'*(W-Wk)).'))))
    subject to
    W>=0;
    X1>=0;
    X2>=0;
    Z>=0;
    diag(W)==1
    for k = 1:N
        r(k) = Omega(W, k);
    end
    r + (1j*tau-delta)*e == Phi(X1) - Gamma(Z, low*pi, up*pi)
    r == Phi(X2)
    cvx_end

    if cvx_optval - loss_old <= epsilon*loss_old
        break
    end
    loss_old = cvx_optval;
    Wk = W;
end

[V, D] = eig(Wk);
[~, idx] = max(diag(D));
s = V(:, idx);
s = s*sqrt(N);
end
% close all
% figure
% spectrum = abs(fft(s, 40*N)).^2;
% spectrum = fftshift(spectrum);
% plot(linspace(-1, 1, length(spectrum)), spectrum)
% figure
% spectrum = abs(fft(r, 40*N)).^2;
% spectrum = fftshift(spectrum);
% plot(linspace(-1, 1, length(spectrum)), spectrum)



function y = Phi(X)
N = size(X, 1);
cvx_begin
expression y(N);
y(1) = Omega(X, 1);
for k = 2:N
    y(k) = Omega(X, k);
end
cvx_end
end

function y = Gamma(Z, low, up)
d0 = cos(low)+cos(up)-cos(up-low)-1;
d1 = (1-exp(1j*low))*(exp(1j*up)-1);

N = size(Z, 1)+1;

cvx_begin
expression y(N)
y(1) = d0*Omega(Z, 1) + conj(d1)*Omega(Z, 2);
for k = 2:N-2
    y(k) = 2*d0*Omega(Z, k) + d1*Omega(Z, k-1) + conj(d1)*Omega(Z, k+1);
end
y(N-1) = d0*Omega(Z, N-1) + d1*Omega(Z, N-2);
y(N) = d1*Omega(Z, N-1);
cvx_end

end

function e = Omega(X, k)
e = 0;
N = size(X, 1);
for l = 0:N-k
    e = e + X(l+k, l+1);
end
end