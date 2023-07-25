function pair = golay_glue(x, y, a, b, c, d)
u = 1/4 * (x+flip(x)+y-flip(y));
v = 1/4 * (x-flip(x)+y+flip(y));
s = kron(a, u) + kron(b, v);
t = kron(a, flip(conj(v))) - kron(b, flip(conj(u)));
e = kron(c, s) + kron(flip(conj(d)), t);
f = kron(flip(conj(c)), t) - kron(d, s);
pair = [e; f];