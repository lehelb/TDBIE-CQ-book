%%
[A,b,c] = RKdata(0); %2stage RK
m = length(c);
z = sym('z');
ons = ones(m,1);

R = 1+z*b'*((eye(m)-z*A)\ons);
C = coeffs(taylor(R-exp(z),'Order',6));
C(end-1)
C(end)
%%
a = sym(sqrt(6));
A = [(88-7*a)/360 (296-169*a)/1800 (-2+3*a)/225;
         (296+169*a)/1800 (88+7*a)/360 (-2-3*a)/225;
        (16-a)/36 (16+a)/36 1/9];
    b = [(16-a)/36 (16+a)/36 1/9].';
    c = [(4-a)/10 ; (4+a)/10 ; 1];
%3stage RK
m = length(c);
z = sym('z');
ons = ones(m,1);

R = 1+z*b'*((eye(m)-z*A)\ons);
C = coeffs(taylor(R-exp(z),'Order',8));
simplify(C(end-1:end))
%%
R = 1/(1-z);
