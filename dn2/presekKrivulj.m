function [P, Q] = presekKrivulj(p, pdot, intp, q, qdot, intq, h)
% function [P, Q] = presekKrivulj(p, pdot, intp, q, qdot, intq, h)
% PARAMETRI:
% p in q sta kazalec na funkciji, ki opisujeta dve ravninsko krivulji,
% pdot in qdot sta kazalca na odvode funkcij p in q
% intp in intq sta intervala, na katerem sta parametrizirani krivulju 
% h je dolzina podintervalov, na katere sta razdeljena intp in intq
% VRNJENE VREDNOSTI:
% P je seznam presecisc K in L
% Q je seznam presecisc lomljenk K' in L'

% Razdelimo interval intp na podintervale dolzine h
interval_p = [intp(1):h:intp(2)];
interval_q = [intq(1):h:intq(2)];

% Izracunamo tocke za lomljenki
K_lomljenka_tocke = p(interval_p);
L_lomljenka_tocke = q(interval_q);

% Izracunamo tocke za krivulji K in L na podanem intervalu
K = p(linspace(intp(1), intp(2)));
L = q(linspace(intq(1), intq(2)));

% Q so presecisca lomljenk, T so parametri pri katerih se lomljenke sekajo
[Q, T] = presecisca(K_lomljenka_tocke, L_lomljenka_tocke, h, intp(1), intq(1));

% Narisemo graf lomljenke K' in L'
subplot(2, 1, 1)
plot(K_lomljenka_tocke(1, :), K_lomljenka_tocke(2, :), 'r')
hold on
plot(L_lomljenka_tocke(1, :), L_lomljenka_tocke(2, :), 'b')
hold on
plot(Q(1, :), Q(2, :), 'x', 'color', 'k')
xlabel("t")
ylabel("p(t)")
legend("K'", "L'")
title("Lomljenki K' in L'")
%T = T(1, :)

% pripravimo matrike za Newtonovo metodo
F = @(x) [p(x(1))(1) - q(x(2))(1); p(x(2))(2) - q(x(1))(2)];

% pripravimo Jacobijevo matriko
JF = @(x) [pdot(x(1))(1) , -qdot(x(2))(1); pdot(x(1))(2), -qdot(x(2))(2)];

st_parametrov = length(T);
P = zeros(2, st_parametrov);

% Vsak parameter v matriki T uporabimo kakor zacetni priblizek za
% Newtonovo iteracijo in rezultat vnesemo v eno od funkcij (p ali q) ter
% dobljeno tocko shranimo v matriko P
for i=1:st_parametrov
  P(:, i) = p(newton(F, JF, T(:, i))(1));
endfor

% Narisemo graf za krivulji K in L
hold on
subplot(2, 1, 2)
plot(K(1, :), K(2, :), 'r')
hold on
plot(L(1, :), L(2, :), 'b')
hold on
plot(P(1, :), P(2, :), 'x', 'color', 'k')
xlabel("t")
ylabel("p(t)")
legend("K", "L")
title("Krivulji K in L")

% testi

%!test
% p = @(t) [t; cos(t)];
% pdot = @(t) [1; -sin(t)];
% intp = [-2*pi, 2*pi];
% q = @(u) [u; sin(u)];
% qdot = @(u) [1; cos(u)];
% intq = [-2*pi, 2*pi];
% h = 0.1;
% P = [-5.4978, -2.3562, 0.7854, 3.9270; 0.7071, -0.7071, 0.7071, -0.7071]
% Q = [-5.4978, -2.3562, 0.7854, 3.9270; 0.7067, -0.7064, 0.7063, -0.7068]
% assert(presekKrivulj(p, pdot, intp, q, qdot, intq, h), [P, Q], eps)

%!test
% p = @(t) [t; cos(t)+sin(t)];
% pdot = @(t) [1; cos(t)-sin(t)];
% intp = [-2*pi, 2*pi];
% q = @(u) [u;u.^2];
% qdot = @(u) [1; 2*u];
% intq = [-2*pi, 2*pi];
% h = 0.1;
% P = [-0.561, 1.15; 0.315, 1.321]
% Q = [-0.561, 1.15; 0.315, 1.321]
% assert(presekKrivulj(p, pdot, intp, q, qdot, intq, h), [P, Q], eps)

%!test
% p = @(t) [t; t.*cos(t)];
% pdot = @(t) [1; cos(t)-t.*sin(t)];
% intp = [-2*pi, 2*pi];
% q = @(u) [u; sin(u)];
% qdot = @(u) [1; cos(u)];
% intq = [-2*pi, 2*pi];
% h = 0.1;
% P = [-4.493, 0.000, 7.725; 0.976, 0.000, 0.992];
% Q = [-4.493, 0.000, 7.725; 0.976, 0.000, 0.992];
% assert(presekKrivulj(p, pdot, intp, q, qdot, intq, h), [P, Q], eps)

%!test
% p = @(t) [t; e.^t];
% pdot = @(t) [1; e.^t];
% intp = [-2, 2];
% q = @(u) [u; u];
% qdot = @(u) [1; 1];
% intq = [-2, 2];
% h = 0.1;
% P = [];
% Q = [];
% assert(presekKrivulj(p, pdot, intp, q, qdot, intq, h), [P, Q], eps)























