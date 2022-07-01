function p = parameters(x, y, model, k)
% function p = parameters(x, y, model, k) za dani linearni model
% izracuna optimalne parametre po metodi najmanjsih kvadratov,
% x in y sta stolpca vrednosti neodvisne in odvisne spremenljivke
% model je funkcija, ki opisuje model (y = model(p, t) vrne 
% vrednosti y, ki jo predvideva model s parametri p v tocki t),
% k je stevilo parametrov v modelu

st_podatkov = size(x)(1);

% zgradimo matriko A
A = zeros(st_podatkov, k);

for i = 1:k
    % pripravimo vektor, ki ima enico samo na i-tem mestu
    % s tem izracunamo potrebne komponente matrike A
    v = zeros(k, 1);
    v(i) = 1;
    
    for j = 1:st_podatkov
        A(j, i) = model(v, x(j));
    endfor
endfor

p = A\y;

%!test
%! model = @(p, x) p(1) + p(2)*x;
%! assert(parameters([0; 1], [1; 1], model, 2), [1; 0], eps)

%!test
%! model = @(p, x) p(1) + p(2)*x;
%! assert(parameters([1;2;3;4], [1;2;1;2], model, 2), [1; 0.2], eps)

%!test
%! model = @(p, x) p(1) + p(2)*x + p(3)*(x^2);
%! assert(parameters([0;1;2;3], [1;2;5;10], model, 3), [1; 0; 1], 2*eps)

%!test
%! model = @(p, x) p(1) + p(2)*x + p(3)*(x^2);
%! assert(parameters([0;1;2], [1;2;5], model, 3), [1; 0; 1], eps)

%!test
%! model = @(p, x) p(1)*sin(x) + p(2)*cos(2*x); 
%! assert(parameters([pi/2; 3*pi/2], [0; -4], model, 2), [2; 2], eps)
