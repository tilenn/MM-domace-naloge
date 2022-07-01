function [P, T] = presecisce(A, B)
%[P, T] = presecisce(A, B) vrne presecisce daljic (oz. prazen stolpec, 
%ce se daljici ne sekata) A1A2 in B1B2 v ravnini.
%A = [A1, A2], A1 in A2 sta krajevna vektorja tock na prvi daljici.
%B = [B1, B2], B1 in B2 sta krajevna vektorja tock na drugi daljici.
%P... krajevni vektor (stolpec) presecisca daljic.
%T = [t; u]... razmerji med razdaljo presecisca 
%do zacetne tocke daljice in celotno dolzino daljice.

%Resimo izpeljani sistem enacb za parametra presecisca obeh premic nosilk.
T = [A(:, 2) - A(:, 1), -(B(:, 2) - B(:, 1))]\(B(:, 1) - A(:, 1));
%Preverimo, ce sta parametra znotraj intervala [0, 1] in ...
if(T(1) <= 1 && T(2) <= 1 && T(1) >= 0 && T(2) >= 0)
	%... vrnemo presecisce, ce sta, oziroma ...
	P = A(:, 1) + T(1)*(A(:, 2) - A(:, 1));
else
	%... vrnemo prazen stolpec, ce nista.
	P = zeros(2, 0);
	T = zeros(2, 0);
end