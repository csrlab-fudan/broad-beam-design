function w= HessMultFcn(x,lambda, v)
load par.mat M N alpha beta factor;
theta = reshape(x(1:end/2), M, N);
phi = reshape(x(end/2+1:end), M, N);
A = exp(1j*theta);
B = exp(1j*phi);
M1 = M*factor; N1 = N*factor;
M2 = round(M1*alpha); % pass band point number
N2= round(N1*beta); % pass band point number
Fa = 1/sqrt(M*N) * fft2(A, M1, N1);
Fb = 1/sqrt(M*N) * fft2(B, M1, N1);
P = zeros(M1, N1);
passband = abs(Fa(1:M2, 1:N2)).^2 + abs(Fb(1:M2, 1:N2)).^2;
P(1:M2, 1:N2) =  passband-mean(passband,"all");
Sa = Fa; Sb = Fb;
Sa(1:M2, 1:N2) = 0;
Sb(1:M2, 1:N2) = 0;
X = reshape(v(1:end/2), M, N);
Y = reshape(v(end/2+1:end), M, N);

% Hessian multiply
[wCons, fftjAX, fftjBY] = cons_Hess_mul(A, B, Sa, Sb, X, Y, M2, N2);
wLoss = loss_Hess_mul(A, B, P, Fa, Fb, X, Y, fftjAX, fftjBY, M2, N2);
w = wLoss + lambda.ineqnonlin*wCons;
end

function [wCons, fftjAX, fftjBY] = cons_Hess_mul(A, B, Sa, Sb, X, Y, M2, N2)
[wTheta, fftjAX] = cons_Hess_mul_sub(A, Sa, X, M2, N2);
[wPhi, fftjBY] = cons_Hess_mul_sub(B, Sb, Y, M2, N2);
wCons = [wTheta; wPhi];
end

function [w,fftjAX] = cons_Hess_mul_sub(A, Sa, X, M2, N2)
[M, N] = size(A);
[M1, N1] = size(Sa);
tmp1 = (M1*N1)/sqrt(M*N) * ifft2(Sa, M1, N1);
tmp1 = tmp1(1:M, 1:N);
fftTmp = 1/sqrt(M*N) * fft2(1j*A.*X, M1, N1);
fftjAX = fftTmp;
fftTmp(1:M2, 1:N2) = 0;
tmp2 = (M1*N1)/sqrt(M*N) * ifft2(fftTmp, M1, N1);
tmp2 = tmp2(1:M, 1:N);
w = 2*imag(-1j*conj(A).*X.*tmp1 + conj(A).*tmp2);
w = w(:);
end

function wLoss = loss_Hess_mul(A, B, P, Fa, Fb, X, Y, fftjAX, fftjBY, M2, N2)
Pax = comp_dpx(Fa, fftjAX, M2, N2);
Pby = comp_dpx(Fb, fftjBY, M2, N2);
w1= loss_Hess_mul_sub(A, P, Fa, X, fftjAX, M2, N2, Pax, Pby);
w2 = loss_Hess_mul_sub(B, P, Fb, Y, fftjBY, M2, N2, Pax, Pby);
wLoss = [w1; w2];
end

function Pax = comp_dpx(Fa, fftjAX, M2, N2)
tmp = 2*real(conj(Fa(1:M2, 1:N2)).*fftjAX(1:M2, 1:N2));
Pax = zeros(size(Fa));
Pax(1:M2, 1:N2) = tmp - mean(tmp,"all");
end

function w = loss_Hess_mul_sub(A, P, Fa, X, fftjAX, M2, N2, Pax, Pby)
[M, N] = size(A);
[M1, N1] = size(Fa);
Pa = zeros(size(Fa));
Pa(1:M2, 1:N2) = Fa(1:M2, 1:N2);
Xa = P.*Fa-sum(P, 'all')/(M2*N2)*Pa;
grad = (M1*N1)/sqrt(M*N) * ifft2(Xa, M1, N1);
grad = grad(1:M, 1:N);
tmp1 = -1j*conj(A).*X.*grad;
tmp21 = (Pax+Pby).*Fa + P.*fftjAX;
fftjAXTmp = zeros(size(fftjAX));
fftjAXTmp(1:M2, 1:N2) = fftjAX(1:M2, 1:N2);
tmp22 = 1/(M2*N2) * (sum(Pax+Pby,"all")*Pa+sum(P,"all")*fftjAXTmp);
tmp2 = (M1*N1)/sqrt(M*N) * ifft2(tmp21-tmp22, M1, N1);
tmp2 = tmp2(1:M, 1:N);
tmp2 = conj(A).*tmp2;
w = 4*imag(tmp1+tmp2);
w = w(:);
end


