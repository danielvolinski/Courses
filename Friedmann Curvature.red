
%% GRAVITATION (MTW)
%% Charles W. MISNER, Kip S. THORNE, John Archibald WHEELER
%% Box 14.5 CURVATURE COMPUTED USING EXTERIOR DIFFERENTIAL
%% FORMS (METRIC FOR FRIEDMANN COSMOLOGY), Page 355

load_package excalc$

pform a=0$
fdomain a=a(t);

write "Define the Friedmann coframe"$
coframe o(t)      = 1                         * d t,
        o(chi)    = a                         * d chi,
        o(theta)  = a * sin(chi)              * d theta,
        o(phi)    = a * sin(chi) * sin(theta) * d phi
  with metric g = - o(t) * o(t) + o(chi) * o(chi) 
                  + o(theta) * o(theta) 
                  + o(phi) * o(phi)$
frame e$
displayframe;
DETM!*;

on fancy$
on nero$
factor o,^$

clear omega$
riemannconx omega$
write "Display the connection form"$
omega(k,-l);

write "Display the connection form in Matrix"$
linelength 100$
coords := {t, chi, theta, phi}$
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

write "Display the Riemann Tensor"$
riemann(a,b,c,d) := e(d) _| ( e(c) _| curv(a,b) );
riemann(a,-b,-c,-d);
riemann(-a,-b,-c,-d);

write "Display the Ricci Tensor"$
ricci(-a,-b) := riemann(c,-a,-d,-b) * g(-c,d);

write "Display the Ricci Scalar"$
riccisc := ricci(-a,-b) * g(a,b);

write "Display the Einstein Tensor"$
clear einstein$
pform einstein(a)=3$
einstein(-a) := (1/2) * curv(b,c) ^ #( o(-b) ^ o(-c) ^ o(-a) );

% Define the velocity and acceleration components.
% Those should bear contravariant indices.
pform {u(k),a(k)} = 0;
remflag('(u a),'covariant);
factor u,a;

write "Display the Geodesic equations"$
a(k) := e(-m) _| omega(-l, k) * u(l) * u(m);
a(-k);

b := killing_vector vec$
write "Killing Vectors Formal expression"$
part(b,1);
write "PDEs for the coefficients"$
for k := 1 step 1 until length(part(b,2)) do write part(part(b,2),k);

showtime;
end;

