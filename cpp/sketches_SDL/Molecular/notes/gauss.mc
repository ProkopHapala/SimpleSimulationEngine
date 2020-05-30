
g(x) := exp(-a*x^2 + -b*(x-s)^2 );

G(x) := integrate(g(x),x);

G(x);

g_(x) := diff(G(x),x,2);

g_(x) - g(x);