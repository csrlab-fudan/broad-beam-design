function W = gen_mllmi_matrix(M, N, alpha, beta, u, v)
s1 = gen_mllmi_sequence(M, alpha, u);
s2 = gen_mllmi_sequence(N, beta, v);
W = s1 * s2.';
W = W/norm(W, "fro");
end