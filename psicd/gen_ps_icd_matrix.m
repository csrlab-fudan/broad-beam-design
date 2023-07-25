function W = gen_ps_icd_matrix(N1, N2, alpha, beta, u, v, rfNum, digitNum)
icdIterNum = 4000;
hybridIterNum = 200;
w1 = gen_ps_icd_ula(N1, alpha, u, icdIterNum, hybridIterNum, rfNum, digitNum);
w2 = gen_ps_icd_ula(N2, beta, v, icdIterNum, hybridIterNum, rfNum, digitNum);
W = w1 * w2.';
end


