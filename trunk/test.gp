Kw=1e-14
h0=1e-2;
w0=Kw/h0;
s0=h0+w0;
d0=h0-w0;

h1=1e-10;
w1=Kw/h1;
s1=h1+w1;
d1=h1-w1;

D=1;

d(t,x) = d1 + (d0-d1) * erfc( x/sqrt(4*D*t) );
s(t,x) = sqrt(s1**2+d(t,x)**2-d1**2);

h(t,x) = 0.5 * (d(t,x)+s(t,x));

pH(t,x) = -log10(h(t,x));
