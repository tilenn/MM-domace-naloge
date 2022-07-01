function X = presekPloskev(f1, gradf1, f2, gradf2, X0, h, n, tol, maxit)
    % function X = presekPloskev(f1, gradf1, f2, gradf2, X0, h, n, tol, maxit)
    % Funkcija, s pomocjo Newtonove metode, vrne n + 1 tock na preseku nivojnic f1 in f2
    % VRNJENA VREDNOST: 
    % X - matrika velikosti 3 x (n + 1), kjer stolpci predstavljajo zaporedne dobljene tocke na preseku nivojnic
    % PARAMETRI:
    % f1, f2 - funkcije ploskev oz. nivojnic
    % gradf1, gradf2 - gradienta funkcij f1 in f2
    % X0 - zacetna aproksimacija tocke na krivulji
    % h - dolzina koraka
    % n - stevilo zaporednih tock na krivulji, ki jih iscemo
    % tol - toleranca za napako pri Newtonovi metodi
    % maxit - maksimalno stevilo iteracij za Newtonovo metodo
    
    % Definiramo vektor tangenten na presek nivojnic
    F = @(x) [cross(gradf1(x), gradf2(x)) / norm(cross(gradf1(x), gradf2(x)))];

    % Definiramo funkcijo f, ki izracuna aproksimacijo za naslednjo tocko
    % Premik za h v smeri vektorja F (dolzina F je enaka 1)
    y = @(x) x + h * F(x);

    % Definiramo sistem nelinearnih enacb, ki ga bomo resili z Newtonovo metodo
    f = @(x) [f1(x); f2(x); F(y(X0))' * x - F(y(X0))' * y(X0)];

    % Definiramo Jacobijevo matriko matrike f
    Df = @(x) [gradf1(x) gradf2(x) F(y(X0))]';

    % Podana tocka X0 ne nujno lezi na preseku nivojnic
    % Problem resimo tako, da resimo sistem in s tem dobimo prvo tocko na preseku 
    % Tocko shranimo v matriko X
    tocka_na_krivulji = newton(f, Df, X0, tol, maxit);
    X = [tocka_na_krivulji];


    % V vsaki iteraciji dobimo novo tocko (skupno 1 + n tock)
    for i = 1:n
        % Iz tocke, izracunane iz prejsnega koraka, se pomaknemo v smeri vektorja, ki je tangenten na isto tocko
        aproksimacija_tocke = y(X(:,i));

        % Ponovno definiramo sistem enacb in Jacobijevo matriko z novimi parametri
        f = @(x) [f1(x); f2(x); F(y(aproksimacija_tocke))' * x - F(y(aproksimacija_tocke))' * y(aproksimacija_tocke)];
        Df = @(x) [gradf1(x) gradf2(x) F(y(aproksimacija_tocke))]';

        % Resimo sistem enacb s pomocjo Newtonove metode
        tocka_na_krivulji = newton(f, Df, aproksimacija_tocke, tol, maxit);

        % Doblejno tocko shranimo
        X = [X, tocka_na_krivulji];
    endfor
    
endfunction

% Newtonova metoda (source: https://ucilnica.fri.uni-lj.si/pluginfile.php/31698/mod_resource/content/2/newton.m)
function [X, n] = newton(F, JF, X0, tol = 1e-10, maxit = 100)
    %X = newton(F, JF, X0, tol, maxit) solves the (nonlinear) 
    %system F(X) = 0 using the Newton's iteration with initial
    %guess X0. (JF is the Jacobi matrix of F.)

    for n = 1:maxit
        %Execute one step of Newton's iteration...
        X = X0 - feval(JF, X0)\feval(F, X0);
        %... and check if the new approximation is within prescribed tolerance.
        if(norm(X - X0) < tol)
            break;
        end
        X0 = X;
    end

    %A warning in case the last approximation is not within specified tolerance.
    if(n == maxit)
        warning("no convergence after maxit iterations")
    end
endfunction

%!test1
%! f1 = @(x)[x(1) + 2*x(2) - 3*x(3) + 1];
%! gradf1 = @(x) [1; 2; -3];
%! f2 = @(x)[3*x(1) + 2*x(2) - 3*x(3) + 1];
%! gradf2 = @(x)[3; 2; -3];
%! X0 = [0; 0; 0.5];
%! K = presekPloskev(f1, gradf1, f2, gradf2, X0, 0.1, 10, 1e-12, 50);
%! for i = 1:columns(K)
%!     assert(f1(K(:,i)), 0, 1e-12);
%!     assert(f2(K(:,i)), 0, 1e-12);
%! endfor

%!test2
%! f1 = @(x)[x(1).^2 + x(2).^2 + x(3).^2 - 4];
%! gradf1 = @(x) [2*x(1); 2*x(2); 2*x(3)];
%! f2 = @(x)[x(1).^2 + x(2).^2 - x(3).^2 - 1];
%! gradf2 = @(x)[2*x(1); 2*x(2); -2*x(3)];
%! X0 = [1.5; 0; 1.3];
%! K = presekPloskev(f1, gradf1, f2, gradf2, X0, 0.1, 10, 1e-12, 50);
%! for i = 1:columns(K)
%!     assert(f1(K(:,i)), 0, 1e-12);
%!     assert(f2(K(:,i)), 0, 1e-12);
%! endfor

%!test3
%! f1 = @(x)[x(1) + 5*x(2) - x(3).^2];
%! gradf1 = @(x) [1; 5; -2*x(3)];
%! f2 = @(x)[x(1)];
%! gradf2 = @(x)[1;0;0];
%! X0 = [0; 0; 0];
%! K = presekPloskev(f1, gradf1, f2, gradf2, X0, 0.5, 7, 1e-12, 50);
%! for i = 1:columns(K)
%!     assert(f1(K(:,i)), 0, 1e-12);
%!     assert(f2(K(:,i)), 0, 1e-12);
%! endfor

