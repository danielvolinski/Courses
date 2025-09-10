
% Based on Narcos Alpha Playlist
% https://www.youtube.com/playlist?list=PLhS8YPbZkfgIxU3dwrqsj30E9FJqNAMOj
% PSI 18/19 - Gravitational Physics Review
% Lectures by Ruth Gregory

load_package excalc$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write "Define Polar coordinates coframe"$
coframe o(r)     = 1 * d r,
        o(theta) = r * d theta
  with metric g = o(r) * o(r) + o(theta) * o(theta)$

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
coords := {r, theta}$
matrix Momega(2, 2)$
for k := 1:2 do for l := 1:2 do
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
matrix Mcurv(2, 2)$
for k := 1:2 do for l := 1:2 do
Mcurv(k, l) := curv(part(coords, k), part(coords, l))$
Mcurv;

write "Display the Riemann Tensor"$
riemann(a,b,c,d) := e(d) _| ( e(c) _| curv(a,b) );

write "Display the Ricci Tensor"$
ricci(-a,-b) := riemann(c,-a,-d,-b) * g(-c,d);

write "Display the Ricci Scalar"$
riccisc := ricci(-a,-b) * g(a,b);

write "Display the Einstein Tensor"$
clear einstein$
pform einstein(a)=3$
einstein(-a) := (1/2) * curv(b,c) ^ #( o(-b) ^ o(-c) ^ o(-a) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write "Define the 2-Sphere coframe"$
coframe o(theta)  =          1 * d theta,
        o(phi)    = sin(theta) * d phi
  with metric g = o(theta) * o(theta) + o(phi) * o(phi)$

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
coords := {theta, phi}$
matrix Momega(2, 2)$
for k := 1:2 do for l := 1:2 do
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
matrix Mcurv(2, 2)$
for k := 1:2 do for l := 1:2 do
Mcurv(k, l) := curv(part(coords, k), part(coords, l))$
Mcurv;

write "Display the Riemann Tensor"$
riemann(a,b,c,d) := e(d) _| ( e(c) _| curv(a,b) );

write "Display the Ricci Tensor"$
ricci(-a,-b) := riemann(c,-a,-d,-b) * g(-c,d);

write "Display the Ricci Scalar"$
riccisc := ricci(-a,-b) * g(a,b);

write "Display the Einstein Tensor"$
clear einstein$
pform einstein(a)=3$
einstein(-a) := (1/2) * curv(b,c) ^ #( o(-b) ^ o(-c) ^ o(-a) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write "Define Spherical coordinates coframe"$
coframe o(r)     = 1              * d r,
        o(theta) = r              * d theta,
        o(phi)   = r * sin(theta) * d phi
  with metric g = o(r) * o(r) + o(theta) * o(theta) + o(phi) * o(phi)$

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
coords := {r, theta, phi}$
matrix Momega(3, 3)$
for k := 1:3 do for l := 1:3 do
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
matrix Mcurv(3, 3)$
for k := 1:3 do for l := 1:3 do
Mcurv(k, l) := curv(part(coords, k), part(coords, l))$
Mcurv;

write "Display the Riemann Tensor"$
riemann(a,b,c,d) := e(d) _| ( e(c) _| curv(a,b) );

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

