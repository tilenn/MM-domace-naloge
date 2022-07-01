function [P, T] = presecisca(A, B, h = 1, a = 0, c = 0)
%[P, T] = presecisca(A, B) poisce presecisca 
%P = [P1, P2,... , Pm] lomljenk 
%A = [A1, A2, ..., Ak] in B = [B1, B2, ..., Bl].
%(Privzetek: Daljice so v genericni legi - presecisca so transverzalna.) 
%Opcijski parametri so korak h ter zacetni tocki intervalov a in c
%za funkcijo presekKrivulj.

%shranimo st. tock na obeh lomljenkah
k = length(A);
l = length(B);
%pripravimo prazna stolpca presecisc P in parametrov T
P = zeros(2, 0);
T = zeros(2, 0);

%poiscemo presecisca vseh parov daljic
for i = 1:(k - 1) % zanka po daljicah iz prve lomljenke (A)
	for j = 1:(l - 1) % zanka po daljicah iz druge lomljenka (B)
		[Q, U] = presecisce(A(:,[i i+1]), B(:, [j j+1]));
		%Ce se daljici sekata, presecisce dodamo v nabor presecisc P ...
		P = [P, Q];
		%... in dodamo se vrednost parametra, pri katerem se sekata, v T.
		if(length(U) != 0)
			T = [T, U*h + [i*h + a; j*h + c]];
		end
	end
end
%P = unique(P.', 'rows').'
%T = unique(T.', 'rows').'