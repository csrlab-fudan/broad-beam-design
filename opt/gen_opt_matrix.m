function [W1, W2, loss] = gen_opt_matrix(M, N, alpha, beta, u, v, k, factor, init)
[x0, initStop] = initialize(M, N, alpha, beta, factor, init); % initial value and stopband leakage
k = k*initStop;
fun = @(x) cost_fun(x, M, N, alpha, beta, factor);
nonlcon = @(x) constraint_fun(x, M, N, alpha, beta, factor, k);

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
checkGrad = false;
givenGrad = true;
maxIter = inf;
maxFunc = inf;
scale = true;
% Hessian = "finite-difference";
display = 'iter';
plotfcn = 'optimplotfval';
% "HessianApproximation", Hessian, 
% "HessianMultiplyFcn",@HessMultFcn
% plotfcn = [];
options = optimoptions('fmincon','Display', display, 'SpecifyObjectiveGradient',givenGrad,...
    'SpecifyConstraintGradient',givenGrad, 'SubproblemAlgorithm','cg', 'UseParallel', true,...
    'PlotFcn',plotfcn, 'CheckGradients',checkGrad, "MaxIterations",maxIter,...
    "MaxFunctionEvaluations", maxFunc, "ScaleProblem", scale,  ...
    "HessianMultiplyFcn",@HessMultFcn);
save("par.mat", "M", "N", "alpha", "beta", "factor");
[x,fval,exitflag,output,lambda,grad,hessian]  = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
h = get(gca, 'Children');
loss = get(h, 'Ydata');
theta = reshape(x(1:end/2), M, N);
phi = reshape(x(end/2+1:end), M, N);
W1= exp(1j*theta);
W2 = exp(1j*phi);
W1 = matrix_shift(W1, alpha, beta, u, v);
W2 = matrix_shift(W2, alpha, beta, u, v);
end

function [x0, initStop] = initialize(M, N, alpha, beta, factor, init)
switch M
    case 96
        P = 2/alpha;
    case 120
        P = 2/alpha;
    otherwise
        P = 1/alpha;
end
b = 1/2;
x1 = gen_general_partial_sequence(M, P, alpha, b);
x2 = gen_general_partial_sequence(N, P, beta, b);
X0 = x1*x2.';

M1 = M*factor; N1 = N*factor;
M2 = round(M1*alpha); % pass band point number
N2 = round(N1*beta); % pass band point number
F = 1/sqrt(M*N) * fft2(X0, M1, N1); 
F(1:M2, 1:N2) = 0;
initStop = 2*(norm(F,"fro")^2);

if isequal(init, 'random')
    x0 = 2*pi*rand(2*M*N, 1);
else
    x0 = mod(angle(X0(:)), 2*pi);
    x0 = [x0; x0];
end
end

function [f, g] = cost_fun(x, M, N, alpha, beta, factor)
theta = reshape(x(1:end/2), M, N);
phi = reshape(x(end/2+1:end), M, N);
A = exp(1j*theta);
B = exp(1j*phi);
M1 = M*factor; N1 = N*factor;
M2 = round(M1*alpha); % pass band point number
N2 = round(N1*beta); % pass band point number
Fa = 1/sqrt(M*N) * fft2(A, M1, N1);
Fb = 1/sqrt(M*N) * fft2(B, M1, N1);
P = zeros(M1, N1);
passband = abs(Fa(1:M2, 1:N2)).^2 + abs(Fb(1:M2, 1:N2)).^2;
P(1:M2, 1:N2) =  passband-mean(passband,"all");
f = norm(P,"fro")^2;

% gradient
if nargout>1
    Pa = zeros(size(Fa)); Pa(1:M2, 1:N2)=Fa(1:M2, 1:N2);
    Pb = zeros(size(Fa)); Pb(1:M2, 1:N2)=Fb(1:M2, 1:N2);
    Xa = P.*Fa-sum(P, 'all')/(M2*N2)*Pa;
    Xb = P.*Fb-sum(P, 'all')/(M2*N2)*Pb;
    ga = 2*(M1*N1)/sqrt(M*N) * ifft2(Xa, M1, N1);
    ga = 2*imag(ga(1:M, 1:N).*conj(A));
    gb = 2*(M1*N1)/sqrt(M*N) * ifft2(Xb, M1, N1);
    gb = 2*imag(gb(1:M, 1:N).*conj(B));
    g = [ga(:); gb(:)];
end
end

function [c, ceq, GC, GCeq] = constraint_fun(x, M, N, alpha, beta, factor, k)
theta = reshape(x(1:end/2), M, N);
phi = reshape(x(end/2+1:end), M, N);
A = exp(1j*theta);
B = exp(1j*phi);
M1 = M*factor; N1 = N*factor;
M2 = round(M1*alpha); % pass band point number
N2 = round(N1*beta); % pass band point number
Fa = 1/sqrt(M*N) * fft2(A, M1, N1); Fa(1:M2, 1:N2) = 0;
Fb = 1/sqrt(M*N) * fft2(B, M1, N1); Fb(1:M2, 1:N2) = 0;
c = norm(Fa,"fro")^2+norm(Fb,"fro")^2-k;
ceq = [];

% gradient
if nargout>2
    ga = (M1*N1)/sqrt(M*N) * ifft2(Fa, M1, N1);
    ga = 2*imag(ga(1:M, 1:N).*conj(A));
    gb = (M1*N1)/sqrt(M*N) * ifft2(Fb, M1, N1);
    gb = 2*imag(gb(1:M, 1:N).*conj(B));
    GC = [ga(:); gb(:)];
    GC = GC(:);
    GCeq = [];
end
end

function A = matrix_shift(A, alpha, beta, u, v)
[N1, N2] = size(A);
shift1 = u-alpha;
shift2 = v-beta;
a1 = exp(1j*pi*(0:N1-1).'*shift1);
a2 = exp(1j*pi*(0:N2-1)*shift2);
A = (a1*a2) .* A;
end