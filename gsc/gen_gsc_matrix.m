function [W, w1, w2] = gen_gsc_matrix(N1, N2, alpha, beta, u, v, flag)
if isequal(flag, "GC")
    P1 = 1;
    P2 = 1;
else
    P1 = find_p(N1, alpha);
    P2 = find_p(N2, beta);
end

w1 = gen_gsc_sequence(N1, alpha, u, P1);
w2 = gen_gsc_sequence(N2, beta, v, P2);
W = w1 * w2.';
W = W/norm(W, "fro");
end

function P = find_p(M, alpha)
switch M
    case 24
        P = 12;
    case 36
        P = 12;
    case 48
        P = 16;
    case 60
        P = 15;
    case 64
        P = 16;
    otherwise
        P = 1/alpha;
end
end