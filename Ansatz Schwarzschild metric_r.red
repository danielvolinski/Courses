
% Lecture 4 The Cartan formalism
% https://www.youtube.com/watch?v=ql1P63zXnFc
% General Rlativity
% Lecture by Ruth Gregory

load_package excalc$

pform {A,B,C}=0$
fdomain A=A(r),B=B(r),C=C(r)$

write "Define the Ansatz Schwarzschild coframe"$
coframe o(t)      = A              * d t, 
        o(r)      = B              * d r,
        o(theta)  = C              * d theta,
        o(phi)    = C * sin(theta) * d phi
  with metric g = + o(t) * o(t) - o(r) * o(r) 
                  - o(theta) * o(theta) - o(phi) * o(phi)$

frame e$
displayframe;
DETM!*;

on fancy;
on nero$
factor o,^$

write "Verify"$
e(-k) _| o(l);

clear omega$
riemannconx omega$
write "Display the connection form"$
omega(k,-l);

write "Display the connection form in Matrix"$
coords := {t, r, theta, phi}$
matrix Momega(4, 4)$
for k := 1:4 do for l := 1:4 do
Momega(k, l) := omega(part(coords, k), part(coords, l))$
Momega;

clear curv,riemann,ricci,riccisc$
pform curv(k,l)=2,{riemann(a,b,c,d),ricci(a,b),riccisc}=0$

index_symmetries curv(k,l): antisymmetric,
                 riemann(k,l,m,n): antisymmetric in {k,l},{m,n}
                                   symmetric in {{k,l},{m,n}},
                 ricci(k,l): symmetric;

write "Display the curvature form"$
curv(k,-l) := d omega(k,-l) + omega(k,-m) ^ omega(m,-l);

write "Display the curvature form in Matrix"$
linelength 200$
matrix Mcurv(4, 4)$
for k := 1:4 do for l := 1:4 do
Mcurv(k, l) := curv(part(coords, k), part(coords, l))$
Mcurv;

write "Display the Riemann Tensor all up"$
riemann(a,b,c,d) := e(d) _| ( e(c) _| curv(a,b) );
write "Display the Riemann Tensor"$
riemann(a,-b,-c,-d);
write "Display the Riemann Tensor all down"$
riemann(-a,-b,-c,-d);

write "Display the Ricci Tensor"$
ricci(-a,-b) := riemann(c,-a,-d,-b) * g(-c,d);

write "Display the Ricci Scalar"$
riccisc := ricci(-a,-b) * g(a,b);

write "Display the Einstein Tensor"$
clear einstein$
pform einstein(a)=3$
einstein(-a) := (1/2) * curv(b,c) ^ #( o(-b) ^ o(-c) ^ o(-a) );

showtime;
end;

