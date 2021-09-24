%% #1. Testing (for the challenge). Squares of side 2, segments of odd length (here I use squares of side 1, segments of length .5, 1.5, ...)
result_all = [];
for L = ((1:15)*2-1)/2
    result = 0;
    for theta = linspace(0,pi/4,30e3)
        points = linspace(0,L,3e3*L)*exp(1j*theta);
        num_squares = numel(unique(floor(points)))+2;
        result = max(result, num_squares);
    end
    result_all(end+1) = result
end

%3,5,6,7,9,10,12,13,15,16,17,19,20,22,23


%% #2. Computing  (for the challenge). Squares of side 2, segments of odd length (here I use squares of side 1, segments of length .5, 1.5, ...)
clear
result_all = [];
for L = ((1:15)*2-1)/2
    t = (0:L).^2;
    [ii, jj] = find(t + t.' < L.^2);
    result_all(end+1) = max(ii+jj)+1; % -2 por 0-based indexing; -1 por los cuadrados hasta (ii,jj);
                                      % +4 por desplazamiento para coger trocitos de cuadrados en los extremos
end


%% #3. Computing. Squares of size 1, segments of length 1, 2, 3, ...
clear
result_all = [];
for L = 1:50
    t = (0:L).^2;
    [ii, jj] = find((t + t.' <= L.^2) & (logical(t)==logical(t.'))); % con igualdad. No contamos las líneas verticales u horizontales, pero sí el origen
    M = max(ii+jj)-1; % -2 por 0-based indexing; -1 por los cuadrados hasta (ii,jj); +2 por desplazamientos para coger trocitos
    [ii, jj] = find(t + t.' < L.^2);
    m = max(ii+jj)+1; % -2 por 0-based indexing; -1 por los cuadrados hasta (ii,jj); +4 por desplazamientos para coger trocitos
    result_all(end+1) = max(M,m); % M: -2 por 0-based indexing; -1 por los cuadrados hasta (ii,jj)
                                            % +4 ó +2 por desplazamientos para coger trocitos. Si M y m son iguales: +4; si no (m==M-1 ó M es vacío): +2
end
% La parte de M nunca aporta nada (he comprobado hasta L=1000). Es porque realmente sólo suma 1: el punto que se añade con "<=",
% el cual no estaba con "<". Pero los de al lado sí estaban: sólo estamos sumando 1. Y luego se
% suman 2 en vez de 4: siempre es peor.


%% #4. Calculando. Pruebo fórmula. ¡Funciona!
L = 1:1e8;
ii = floor(L/sqrt(2));
jj = ceil(sqrt(L.^2-ii.^2))-1;
result_all_form = ii+jj+3;
%plot(ii, jj, '.-')


%% #5. Testing randomly. Squares of side 1, segments of length 1, 2, 3, ...
L = 50;
N = 1e5;
result = NaN(1,N);
for n = 1:N
    z0 = rand + 1j*rand;
    theta = rand*pi/4;
    points = linspace(0,L,.1e3*L)*exp(1j*theta)+z0;
    num_squares = numel(unique(floor(points)));
    result(n) = num_squares;
end
%histogram(result)
max(result)


%% #6. Calculando. Pruebo (pág. 5). Funciona
result_all = [];
for L = 1:1e5
    M = 1:2*L+1;
    result_all(end+1) = find(M.^2-6*M+10-mod(M,2) < 2*L^2, 1, 'last');
end


%% #7. Calculando. Fórmula de @xnor. ¡Funciona! He probado hasta 1e8
result_all = floor(sqrt(2*L.^2-2))+3; %


%% #8. Problema inverso
M = 3:1e3;
L_1 = floor(sqrt((M.^2-6*M+10-mod(M,2))/2))+1; % mi referencia (pág. 3)
L_2 = ceil(sqrt((M-3).^2/2+1)); % "invirtiendo" fórmula para M de @xnor. ¡Funciona!


%% #9. Figura para OEIS

close all
L_all = [1 2 3];
z0_all = [.8+.55j .78+2.9j 1.92+3.02j];
theta_all = [35 48 -45.3]/180*pi;
for k = 1:numel(z0_all)
    z0 = z0_all(k);
    theta = theta_all(k);
    L = L_all(k);
    plot([z0 z0+L*exp(1j*theta)], 'linewidth', 1.5);
    hold on
end
grid on, axis([0 5 0 5]), axis square
set(gca, 'GridAlpha', .4, 'XTick', 1:5, 'YTick', 1:5)
box on


%% #10. Figures with examples, for the paper
clf
hold on
sq = [1+1j 1+2j 2+2j   4+4j 4+3j 4+2j 5+2j 5+1j   7+4j 8+4j 8+3j 9+3j 9+2j 10+2j 10+1j];
for k = 1:numel(sq)
    patch(real([sq(k) sq(k)+1 sq(k)+1 sq(k) sq(k)]), ...
          imag([sq(k) sq(k) sq(k)+1j sq(k)+1j sq(k)]), ...
          .88*[1 1 1], 'edgecolor', 'none')
end
set(gca, 'layer', 'top')
xlabel x
ylabel y

z1 = 1.3+1.8j;
r = 1;
%z1 = 1.55+1.9j;
%r = .8;
theta = 30/180*pi;
plot(z1+[0 r*exp(1j*theta)], '-', 'linewidth', 1)

z1 = 4.7+4.08j;
r = 2.2;
theta = -78.5/180*pi;
set(gca, 'colororderindex', 1)
plot(z1+[0 r*exp(1j*theta)], '-', 'linewidth', 1)

z1 = 7.9+4.25j;
r = 3.5;
theta = -42.5/180*pi;
set(gca, 'colororderindex', 1)
plot(z1+[0 r*exp(1j*theta)], '-', 'linewidth', 1)

axis equal
grid on, axis([0 12 0 6])
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4, 'XTick', xl(1):xl(2), 'YTick', yl(1):yl(2))


%% #11. Figure: touched squares, no grid points
clf
hold on
sq = [0+0j 1+0j 1+1j 2+1j 3+1j 3+2j 4+2j 4+3j];
for k = 1:numel(sq)
    patch(real([sq(k) sq(k)+1 sq(k)+1 sq(k) sq(k)]), ...
          imag([sq(k) sq(k) sq(k)+1j sq(k)+1j sq(k)]), ...
          .88*[1 1 1], 'edgecolor', 'none')
end
set(gca, 'layer', 'top')
patch([0 5 5 0 0], [0 0 4 4 0], 'k', 'linestyle', '--', 'facecolor', 'none', 'linewidth', 1);
z1 = .6 +.2j ;
z2 = 4.8 + 3.1i;
plot([z1 z2], '-', 'linewidth', 1)
axis equal
grid on, axis([-1 6 -1 5])
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4, 'XTick', xl(1):xl(2), 'YTick', yl(1):yl(2))
xlabel x
ylabel y


%% #12. Figure: touched squares, one grid point
clf
hold on
sq = [0+0j 1+0j 1+1j 2+1j 3+2j 4+2j 4+3j];
for k = 1:numel(sq)
    patch(real([sq(k) sq(k)+1 sq(k)+1 sq(k) sq(k)]), ...
          imag([sq(k) sq(k) sq(k)+1j sq(k)+1j sq(k)]), ...
          .88*[1 1 1], 'edgecolor', 'none')
end
set(gca, 'layer', 'top')
patch([0 5 5 0 0], [0 0 4 4 0], 'k', 'linestyle', '--', 'facecolor', 'none', 'linewidth', 1);
z1 = .6 +.2j;
z2 = 4.8 + 3.35i;
plot([z1 z2], '-', 'linewidth', 1)
axis equal
grid on, axis([-1 6 -1 5])
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4, 'XTick', xl(1):xl(2), 'YTick', yl(1):yl(2))
xlabel x
ylabel y


%% #13. Figure: L ineq sqrt, i=4, j=3
clf
hold on
sq = [0+0j  3+2j];
for k = 1:numel(sq)
    patch(real([sq(k) sq(k)+1 sq(k)+1 sq(k) sq(k)]), ...
          imag([sq(k) sq(k) sq(k)+1j sq(k)+1j sq(k)]), ...
          .88*[1 1 1], 'edgecolor', 'none')
end
plot([0+1j 1+1j], 'k:', 'linewidth', .75)
plot([1+0j 1+1j], 'k:', 'linewidth', .75)
plot([3+2j 4+2j], 'k:', 'linewidth', .75)
plot([3+2j 3+3j], 'k:', 'linewidth', .75)
set(gca, 'layer', 'top')
patch([0 4 4 0 0], [0 0 3 3 0], 'k', 'linestyle', '--', 'facecolor', 'none', 'linewidth', 1);
z1 = .8 +.95j ;
z2 = 3.2 + 2.05i;
plot([z1 z2], '-', 'linewidth', 1)
z1 = .05+.1j;
z2 = 3.9+2.95j;
set(gca, 'colororderindex', 1)
plot([z1 z2], '-', 'linewidth', 1)
axis equal
grid on, axis([-1 5 -1 4])
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4, 'XTick', xl(1):xl(2), 'YTick', yl(1):yl(2))
xlabel x
ylabel y


%% #14. Figure: L ineq sqrt, i=4, j=2
clf
hold on
sq = [0+0j  3+1j];
for k = 1:numel(sq)
    patch(real([sq(k) sq(k)+1 sq(k)+1 sq(k) sq(k)]), ...
          imag([sq(k) sq(k) sq(k)+1j sq(k)+1j sq(k)]), ...
          .88*[1 1 1], 'edgecolor', 'none')
end
plot([0+1j 1+1j], 'k:', 'linewidth', .75)
plot([1+0j 1+1j], 'k:', 'linewidth', .75)
plot([3+1j 4+1j], 'k:', 'linewidth', .75)
plot([3+1j 3+2j], 'k:', 'linewidth', .75)
set(gca, 'layer', 'top')
patch([0 4 4 0 0], [0 0 2 2 0], 'k', 'linestyle', '--', 'facecolor', 'none', 'linewidth', 1);
z1 = .85 +.95j ;
z2 = 3.1 + 1.15i;
plot([z1 z2], '-', 'linewidth', 1)
z1 = .04+.06j;
z2 = 3.93+1.95j;
set(gca, 'colororderindex', 1)
plot([z1 z2], '-', 'linewidth', 1)
axis equal
grid on, axis([-1 5 -1 4])
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4, 'XTick', xl(1):xl(2), 'YTick', yl(1):yl(2))
xlabel x
ylabel y


%% #15. Figure: i,j,L,S
clf
hold on
set(gca, 'layer', 'top')
axis equal
grid on
plot([2+2j 2+7.6j], 'k')
plot([2+2j 7.6+2j], 'k')
plot([2+2j 2+3j 3+3j 3+4j 4+4j 4+5j], 'ko', 'markerfacecolor', 'k', 'markersize', 5)
plot([4+2j 2+4j 5+2j 2+5j 3+5j 5+3j 2+6j 6+2j 3+6j 6+3j 7+2j 2+7j 3+2j 4+3j 5+4j], 'ko', 'markerfacecolor', 'none', 'markersize', 5)
theta = linspace(0,pi/2,200);
for r = [1 sqrt(2) sqrt(5) sqrt(8) sqrt(13)] 
    set(gca, 'colororderindex', 1)
    plot(2+2j+r*exp(j*theta))
end
plot([1.8 2.4], [2.2 1.6], 'k--')
plot([1.8 3.4], [3.2 1.6], 'k--')
plot([1.8 4.4], [4.2 1.6], 'k--')
plot([1.8 5.4], [5.2 1.6], 'k--')
plot([1.8 6.4], [6.2 1.6], 'k--')
plot([1.8 7.4], [7.2 1.6], 'k--')
axis([0 7 0 7])
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4, 'XTick', xl(1):xl(2), 'YTick', yl(1):yl(2))
xlabel i, ylabel j
text(1.4333, 1.3597, ['i+j' char(8211) '1 = 3         4          5         6          7         8'])
axis([0 8 0 8])


%% #16. Direct formula for funt, for L not necessarily integer: checked
clear, clf
L_step = .001;
L_max = 14;
L = L_step:L_step:L_max;
i = ceil(L/sqrt(2)) + 1;
j = ceil(sqrt(L.^2-(i-2).^2)) + 1;
funt = i+j-1; 
plot(L, funt, 'linewidth', .6)
L = 1:L_max;
i = ceil(L/sqrt(2)) + 1;
j = ceil(sqrt(L.^2-(i-2).^2)) + 1;
funt = i+j-1; 
hold on
set(gca, 'colororderindex', 1)
plot(L, funt, '.', 'markersize', 7, 'linewidth', .6)
grid on
set(gcf, 'position', [318 263 618 333])
xlabel(char(8467))
ylabel([char(963) '(' char(8467) ')'])


%% #17. Formula for funt using funl, for L not necessarily integer: checked
L_step = .001;
L_max = 14;
L = L_step:L_step:L_max;
n = 4:2*L_max;
funl = [0 0 0 sqrt((n.^2-6*n+10-mod(n,2))/2)];
for k = 1:numel(L)
    funt(k) = find(funl < L(k), 1, 'last');    
end


%% #18. @xnor's formula for L integer: checked
L = 1:14;
funt = floor(sqrt(2*L.^2-2))+3;
plot(L, funt, '.')


%% #19. Formula for funl:
M = 1:20;
funl = [0 0 0 sqrt((M(4:end).^2-6*M(4:end)+10-mod(M(4:end),2))/2)];
stem(M, funl, '.', 'markersize', 9, 'linewidth', .6)
grid on
set(gcf, 'position', [318 263 618 333])
xlabel S
ylabel([char(955) '(S)'])


%% #20. Formula for funli:
S = 1:20;
%funl = [0 0 0 sqrt((S(4:end).^2-6*S(4:end)+10-mod(S(4:end),2))/2)];
%funli = floor(funl)+1;
funli = [1 1 1 ceil(sqrt((S(4:end)-3).^2/2 + 1))];
stem(S, funli, '.', 'markersize', 9, 'linewidth', .6)
grid on
set(gcf, 'position', [318 263 618 333])
xlabel S
ylabel([char(923) '(S)'])


%% #21. Rectangular case, initial tests
clear
clf
a = 1.35; b = 1;
d_max = 20;
[ii, jj] = ndgrid(0:d_max/min(a,b));
ii = ii(:)+2; jj = jj(:)+2;
d2 = ((ii-2)*a).^2 + ((jj-2)*b).^2;
ind = d2<=d_max^2;
ii = ii(ind); jj = jj(ind); d2 = d2(ind);
ij1 = ii+jj-1;
data = [ii jj ij1 d2];
data = sortrows(data, [3 4]);
ind = find(diff([0; data(:, 3)])~=0); % first (minimum) d2 for each ij1. If there are coincidences only one point is taken
data_chosen = data(ind, :);
%data_chosen
plot(data_chosen(:,1)-2, data_chosen(:,2)-2, 'o'), %axis equal
hold on
xl = [0 max(data_chosen(:,1))-2];
%plot(xl, (2*a^2*xl+a^2-b^2)/2/b^2)
%plot(xl, (2*a^2*xl)/2/b^2)
%plot(xl, (2*a^2*(xl-1)+a^2-b^2)/2/b^2)
k = a^2/b^2; if mod(k,1)<1e-9, k = fix(k); end % fix floating-point inaccuracy
plot(xl, k*xl-ceil(k/2)) % checked. Valid only for k (=a^2/b^2) integer
grid on

    
%% #22. For a >= b: values of jj that are repeated because ii increases instead of jj:
a = 1.35; b = 1; ii=3:6; jj_rep = ceil((2*a^2*(ii-2)+a^2-b^2)/2/b^2)+2


%% #23. Page 11. k (=a^2/b^2) even. Checked
k = 4;
b = 1;
a = sqrt(k*b^2);
L = 0.01:.01:20;
iip = ceil(real( (k^2 + sqrt(-k^3 + 4*k*(k+1)*L.^2/b^2)) /2/k/(k+1) ))-1;
jjp = ceil( sqrt(L.^2-(iip*a).^2)/b )-1;
plot(iip, jjp, 'o')


%% #24. Page 12. k (=a^2/b^2) odd. Checked
k = 5;
b = 1;
a = sqrt(k*b^2);
L = 0.01:.01:30;
iip = ceil(real( (k+k^2 + sqrt(-k^3-2*k^2-k + 4*k*(k+1)*L.^2/b^2)) /2/k/(k+1) ))-1;
jjp = ceil( sqrt(L.^2-(iip*a).^2)/b )-1;
plot(iip, jjp, 'o')


%% #25. Computation, unifying odd and even:
clear
k = 7;
L = .001:.001:14; %4.7;
b = 1;
a = sqrt(k*b^2);
if mod(k,1) % odd
    iip = ceil(real( (k+k^2 + sqrt(-k^3-2*k^2-k + 4*k*(k+1)*L.^2/b^2)) /2/k/(k+1) ))-1;
else % even
    iip = ceil(real( (k^2 + sqrt(-k^3 + 4*k*(k+1)*L.^2/b^2)) /2/k/(k+1) ))-1;
end
jjp = ceil( sqrt(L.^2-(iip*a).^2)/b )-1;
S = iip+jjp+3
%unique([iip(:)+2 jjp(:)+2], 'rows')
plot(L, S)


%% #26. Testing randomly with rectangles and real-valued lengths: sampling the segment and using the canonical rectangle
k = 4;
L = 50;
b = 1;
a = sqrt(k*b^2);
N = 1e4;
P = 1e3; % number of points to sample the segment
result1 = NaN(1,N);
result2 = NaN(1,N);
disp('Simulating...')
for n = 1:N
    z0 = rand + 1j*rand;
    theta = rand*pi/2; % pi/4 is no longer valid for a rectangular grid; pi/2 is needed
    z1 = z0 + L*exp(1j*theta);
    points = linspace(0,L,P*L)*exp(1j*theta)+z0; % sample segment finely. Endpoints are always included
    points = real(points)/a + 1j*imag(points)/b; % normalize by cell dimensions
    num_cells = numel(unique(floor(points)));
    result1(n) = num_cells;
    result2(n) = ceil(real(z1)/a)-floor(real(z0)/a) + ceil(imag(z1)/b)-floor(imag(z0)/b) - 1; % using canonical rectangle
    % The theoretical probability that the result is less than this (because the segment passes
    % through a grid point) is 0
end
%histogram(result1)
max(result1), max(result2) % same results: checked (by trying several times). result1 sometimes gives
% less than result2; never more. The occurrrence of this decreases with P. This seems logical


%% #27. Testing randomly with rectangles and real-valued lengths: using the canonical rectangle. This allows vectorization
clear
k = 3;
L = 27;
b = 1;
a = sqrt(k*b^2);
N = 1e8;
disp('Simulating...')
z0 = rand(1,N) + 1j*rand(1,N);
theta = rand(1,N)*pi/2; % pi/4 is no longer valid for a rectangular grid; pi/2 is needed
z1 = z0 + L*exp(1j*theta);
result = ceil(real(z1)/a)-floor(real(z0)/a) + ceil(imag(z1)/b)-floor(imag(z0)/b) - 1; % using canonical rectangle
%histogram(result, 'binmethod', 'integer')
%ind = result==max(result);
%unique(ceil(real(z1(ind))/a)-floor(real(z0(ind))/a) + 1j*ceil(imag(z1(ind))/b)-floor(imag(z0(ind))/b) )
max(result)
% Sometimes this gives less than than #25, specially for large L. But checking iip, jjp from #25,
% those iip and jjp are seen to be correct: the simulation simply didn't cover that case, but it can
% be created manually with the iip, jjp from #25, see #28. Or the simulation can be forced to
% sample around those valuesm see #29


%% #28: plotting the result from #25
clear
clf
k = 3;
L = 27;
b = 1;
a = sqrt(k*b^2);
iip = 8; jjp = 23; % from #25
z0 = 0;
theta = atan2(jjp*b, iip*a);
z1 = z0 + L*exp(1j*theta);
plot([z0 z1])
grid on
xticks(0:a:max(real(z1))+a), yticks(0:b:max(imag(z1))+b)


%% #29. Testing randomly but forcing values around the iip, jpp from #25, to check that the number given by #25 is actually reached
clear
k = 3;
L = 27;
b = 1;
a = sqrt(k*b^2);
N = 1e8;
disp('Simulating...')
z0 = (.98+.02*rand(1,N))*a + (.98j+.02j*rand(1,N))*b;
iip = 8; jjp = 23; % insert values from #25
theta = atan2(jjp*b, iip*a) -.05 + .1*rand(1,N)*pi/2;
z1 = z0 + L*exp(1j*theta);
result = ceil(real(z1)/a)-floor(real(z0)/a) + ceil(imag(z1)/b)-floor(imag(z0)/b) - 1; % using canonical rectangle
%histogram(result, 'binmethod', 'integer')
%ind = result==max(result);
%unique(ceil(real(z1(ind))/a)-floor(real(z0(ind))/a) + 1j*ceil(imag(z1(ind))/b)-floor(imag(z0(ind))/b) )
max(result)
% Sometimes this gives less than than #25, specially for large L. But checking iip, jjp from #25,
% those iip and jjp are seen to be correct: the simulation simply didn't cover that case, but it can
% be created manually with the iip, jjp from #25, as in #28:


%% #30. Formulas for iip for k=a^2/b^2 integer, without using k:
clear, clc
k = 4.3;
L = 345;
b = 1.3;
a = sqrt(k*b^2);
if mod(k,2) % odd, without ceil. From #25
    iip = real( (k+k^2 + sqrt(-k^3-2*k^2-k + 4*k*(k+1)*L.^2/b^2)) /2/k/(k+1) ) - 1
    iip_new = (a*(a^2+b^2) + b*real(sqrt(-(a^2+b^2)^2+4*L^2*(a^2+b^2)))) / 2/a/(a^2+b^2) - 1
    iip_new2 = (a + b*real(sqrt(4*L^2/(a^2+b^2)-1))) / 2/a - 1
else % even, without ceil. From #25
    iip = real( (k^2 + sqrt(-k^3 + 4*k*(k+1)*L.^2/b^2)) /2/k/(k+1) ) -1
    iip_new = (a^3 + b*real(sqrt(-a^4+4*L^2*(a^2+b^2)))) / 2/a/(a^2+b^2) - 1 % checked
end


%% #31. Figure: i,j,L,S, rectangular, 1.35
clf
a = 1.35; b = 1;
hold on
set(gca, 'layer', 'top')
grid on
plot([2*a+2j*b 2*a+8.7j*b], 'k') % "secondary" axes
plot([2*a+2j*b 7.7*a+2j*b], 'k')
chosen = [2*a+2j*b 2*a+3j*b 3*a+3j*b 3*a+4j*b 3*a+5j*b 4*a+5j*b 4*a+6j*b 4*a+7j*b 5*a+7j*b 5*a+8j*b 6*a+8j*b];
plot(chosen, 'ko', 'markerfacecolor', 'k', 'markersize', 5) % chosen (i,j) pairs
plot(setdiff((2:6)*a + 1j*(2:8).'*b, chosen), 'ko', 'markerfacecolor', 'none', 'markersize', 5) % non-chosen (i,j) pairs
axis equal %pbaspect([1 1 1])
axis([0 8*a 0 9*b])
theta = linspace(0,pi/2,200);
for r = [ b hypot(a,(1:3)*b) hypot(2*a,(3:5)*b) hypot(3*a,(5:6)*b) hypot(4*a,6*b) ] 
    set(gca, 'colororderindex', 1)
    z = 2*a+2j*b+r*exp(1j*theta);
    z(real(z)> 7.4*a | imag(z)>8.5*b) = [];
    plot(z)
end
plot([1.8 2.4]*a, [2.2 1.6]*b, 'k--')
plot([1.8 3.4]*a, [3.2 1.6]*b, 'k--')
plot([1.8 4.4]*a, [4.2 1.6]*b, 'k--')
plot([1.8 5.4]*a, [5.2 1.6]*b, 'k--')
plot([1.8 6.4]*a, [6.2 1.6]*b, 'k--')
plot([1.8 7.4]*a, [7.2 1.6]*b, 'k--')
plot([1.8 7.4]*a, [8.2 2.6]*b, 'k--') % limit to secondary x-axis
plot([2.5 7.4]*a, [8.5 3.6]*b, 'k--') % limit to secondary x-axis and y-axis
plot([3.5 7.4]*a, [8.5 4.6]*b, 'k--')
plot([4.5 7.4]*a, [8.5 5.6]*b, 'k--')
plot([5.5 7.4]*a, [8.5 6.6]*b, 'k--')
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4, 'XTick', xl(1):xl(2), 'YTick', yl(1):yl(2))
xticks((0:7)*a), yticks((0:9)*b)
xticklabels([{'0'} arrayfun(@(n)[num2str(n) char(8201) 'a'], 1:7, 'UniformOutput', false)]) % thin space
yticklabels([{'0'} arrayfun(@(n)[num2str(n) char(8201) 'b'], 1:9, 'UniformOutput', false)])
xlabel(['i' char(8201) 'a']), ylabel(['j' char(8201) 'b']) % thin space
text(2.11, 1.35, ['i+j' char(8211) '1 ' char(8201) '= ' char(8201) '3            4            5            6            7           8'])
plot(xl*a, ( a^2/b^2*xl - 5*a^2/b^2/2 + 3/2 )*b, '--') % lower-limit line, with axes i*a and i*b


%% #32. Figure: i,j,L,S, "rectangular", 1
clf
hold on
set(gca, 'layer', 'top')
axis equal
grid on
plot([2+2j 2+7.6j], 'k')
plot([2+2j 7.6+2j], 'k')
plot([2+2j 3+2j 3+3j 4+3j 4+4j 5+4j], 'ko', 'markerfacecolor', 'k', 'markersize', 5)
plot([4+2j 2+4j 5+2j 2+5j 3+5j 5+3j 2+6j 6+2j 3+6j 6+3j 7+2j 2+7j 2+3j 3+4j 4+5j], 'ko', 'markerfacecolor', 'none', 'markersize', 5)
theta = linspace(0,pi/2,200);
for r = [1 sqrt(2) sqrt(5) sqrt(8) sqrt(13)] 
    set(gca, 'colororderindex', 1)
    plot(2+2j+r*exp(1j*theta))
end
plot([1.8 2.4], [2.2 1.6], 'k--')
plot([1.8 3.4], [3.2 1.6], 'k--')
plot([1.8 4.4], [4.2 1.6], 'k--')
plot([1.8 5.4], [5.2 1.6], 'k--')
plot([1.8 6.4], [6.2 1.6], 'k--')
plot([1.8 7.4], [7.2 1.6], 'k--')
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4)
xticks((0:7)), yticks((0:9))
xticklabels([{'0'} arrayfun(@(n)[num2str(n) char(8201) 'a'], 1:7, 'UniformOutput', false)]) % thin space
yticklabels([{'0'} arrayfun(@(n)[num2str(n) char(8201) 'b'], 1:9, 'UniformOutput', false)])
xlabel(['i' char(8201) 'a']), ylabel(['j' char(8201) 'b']) % thin space
xlabel(['i' char(8201) 'a']), ylabel(['j' char(8201) 'b']) % thin space
text(1.4333, 1.3597, ['i+j' char(8211) '1 = 3         4          5         6          7         8'])
axis([0 8 0 8])


%% #32b. Figure: i,j,L,S, "rectangular", 1. More points than #32: like #31
clf
hold on
set(gca, 'layer', 'top')
axis equal
grid on
plot([2+2j 2+6.6j], 'k')
plot([2+2j 7.6+2j], 'k')
plot([2+2j 3+2j 3+3j 4+3j 4+4j 5+4j 5+5j 6+5j 6+6j], 'ko', 'markerfacecolor', 'k', 'markersize', 5)
plot([4+2j 2+4j 5+2j 2+5j 3+5j 5+3j 2+6j 6+2j 3+6j 6+3j 2+3j 3+4j 4+5j 4+6j 6+4j 5+6j], 'ko', 'markerfacecolor', 'none', 'markersize', 5)
theta = linspace(0,pi/2,200);
for r = [1 sqrt(2) sqrt(5) sqrt(8) sqrt(13) sqrt(18) 5 sqrt(32)] 
    set(gca, 'colororderindex', 1)
    z = 2+2j+r*exp(1j*theta);
    z(real(z)> 7.4*a | imag(z)>6.5*b) = [];
    plot(z)
end
plot([1.8 2.4], [2.2 1.6], 'k--')
plot([1.8 3.4], [3.2 1.6], 'k--')
plot([1.8 4.4], [4.2 1.6], 'k--')
plot([1.8 5.4], [5.2 1.6], 'k--')
plot([1.8 6.4], [6.2 1.6], 'k--')
plot([2.5 7.4]*a, [6.5 1.6]*b, 'k--') % limit to secondary x-axis and y-axis
plot([3.5 7.4]*a, [6.5 2.6]*b, 'k--')
plot([4.5 7.4]*a, [6.5 3.6]*b, 'k--')
plot([5.5 7.4]*a, [6.5 4.6]*b, 'k--')
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4)
xticks((0:7)), yticks((0:9))
xticklabels([{'0'} arrayfun(@(n)[num2str(n) char(8201) 'a'], 1:7, 'UniformOutput', false)]) % thin space
yticklabels([{'0'} arrayfun(@(n)[num2str(n) char(8201) 'b'], 1:9, 'UniformOutput', false)])
xlabel(['i' char(8201) 'a']), ylabel(['j' char(8201) 'b']) % thin space
xlabel(['i' char(8201) 'a']), ylabel(['j' char(8201) 'b']) % thin space
text(1.4333, 1.3597, ['i+j' char(8211) '1 = 3         4          5         6          7         8'])
axis([0 8 0 7])


%% #33. Figure: touched squares, no grid points / one grid point
clf
hold on
a = 1.35; b = 1; 
%sq = [0+0j 1+0j 1+1j 2+1j 3+1j 3+2j 4+2j 4+3j]; % no grid points
sq = [0+0j 1+0j 1+1j 2+1j 3+2j 4+2j 4+3j]; % grid points
for k = 1:numel(sq)
    patch(real([sq(k) sq(k)+1 sq(k)+1 sq(k) sq(k)])*a, ...
          imag([sq(k) sq(k) sq(k)+1j sq(k)+1j sq(k)])*b, ...
          .88*[1 1 1], 'edgecolor', 'none')
end
set(gca, 'layer', 'top')
patch([0 5 5 0 0]*a, [0 0 4 4 0]*b, 'k', 'linestyle', '--', 'facecolor', 'none', 'linewidth', 1);
z1 = .6*a +.2j*b ;
%z2 = 4.8*a + 3.1j*b; % no grid points
z2 = 4.8*a + 3.35j*b; % grid points
plot([z1 z2], '-', 'linewidth', 1)
axis equal
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4)
xticks((-1:7)*a), yticks((-1:9)*b)
grid on, axis([-1*a 6*a -1*b 5*b])
xlabel x
ylabel y


%% #34. Figures with examples, for the paper; rectangular
clf
a = 1.35; b = 1;
hold on
sq = [1+1j 1+2j 2+2j   4+4j 4+3j 4+2j 5+2j 5+1j   6+5j 6+4j 7+4j 7+3j 8+3j 8+2j 8+1j 9+1j];
for k = 1:numel(sq)
    patch(real([sq(k) sq(k)+1 sq(k)+1 sq(k) sq(k)])*a, ...
          imag([sq(k) sq(k) sq(k)+1j sq(k)+1j sq(k)])*b, ...
          .88*[1 1 1], 'edgecolor', 'none')
end
set(gca, 'layer', 'top')
xticks((-1:11)*a), yticks((-1:7)*b)
xlabel x
ylabel y

z1 = 1.49*a+1.8j*b;
r = 1;
%z1 = 1.55+1.9j;
%r = .8;
theta = 30/180*pi;
plot(z1+[0 r*exp(1j*theta)], '-', 'linewidth', 1)

z1 = 4.75*a+4.12j*b;
r = 2.4;
theta = -78.5/180*pi;
set(gca, 'colororderindex', 1)
plot(z1+[0 r*exp(1j*theta)], '-', 'linewidth', 1)

z1 = 6.8*a+5.1j*b;
r = 4.7;
theta = -48.5/180*pi;
set(gca, 'colororderindex', 1)
plot(z1+[0 r*exp(1j*theta)], '-', 'linewidth', 1)

axis equal
grid on, axis([0 11*a 0 7])
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4)


%% #35. Figure: L ineq sqrt, i=4, j=3; rectangular
clf
a = 1.35; b = 1;
hold on
sq = [0+0j  3+2j];
for k = 1:numel(sq)
    patch(real([sq(k) sq(k)+1 sq(k)+1 sq(k) sq(k)])*a, ...
          imag([sq(k) sq(k) sq(k)+1j sq(k)+1j sq(k)])*b, ...
          .88*[1 1 1], 'edgecolor', 'none')
end
plot([0*a+1j*b 1*a+1j*b], 'k:', 'linewidth', .75)
plot([1*a+0j*b 1*a+1j*b], 'k:', 'linewidth', .75)
plot([3*a+2j*b 4*a+2j*b], 'k:', 'linewidth', .75)
plot([3*a+2j*b 3*a+3j*b], 'k:', 'linewidth', .75)
set(gca, 'layer', 'top')
patch([0 4 4 0 0]*a, [0 0 3 3 0]*b, 'k', 'linestyle', '--', 'facecolor', 'none', 'linewidth', 1);
z1 = .8*a +.95j*b ;
z2 = 3.2*a + 2.05j*b;
plot([z1 z2], '-', 'linewidth', 1)
z1 = .05*a+.1j*b;
z2 = 3.9*a+2.95j*b;
set(gca, 'colororderindex', 1)
plot([z1 z2], '-', 'linewidth', 1)
axis equal
grid on, axis([-1*a 5*a -1*b 4*b])
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4, 'XTick', xl(1):a:xl(2), 'YTick', yl(1):b:yl(2))
xlabel x
ylabel y


%% #36. Figure: L ineq sqrt, i=4, j=2; rectangular
clf
a = 1.35; b = 1;
hold on
sq = [0+0j  3+1j];
for k = 1:numel(sq)
    patch(real([sq(k) sq(k)+1 sq(k)+1 sq(k) sq(k)])*a, ...
          imag([sq(k) sq(k) sq(k)+1j sq(k)+1j sq(k)])*b, ...
          .88*[1 1 1], 'edgecolor', 'none')
end
plot([0*a+1j*b 1*a+1j*b], 'k:', 'linewidth', .75)
plot([1*a+0j*b 1*a+1j*b], 'k:', 'linewidth', .75)
plot([3*a+1j*b 4*a+1j*b], 'k:', 'linewidth', .75)
plot([3*a+1j*b 3*a+2j*b], 'k:', 'linewidth', .75)
set(gca, 'layer', 'top')
patch([0 4 4 0 0]*a, [0 0 2 2 0]*b, 'k', 'linestyle', '--', 'facecolor', 'none', 'linewidth', 1);
z1 = .85*a +.95j*b ;
z2 = 3.1*a + 1.15j*b;
plot([z1 z2], '-', 'linewidth', 1)
z1 = .04*a+.06j*b;
z2 = 3.93*a+1.95j*b;
set(gca, 'colororderindex', 1)
plot([z1 z2], '-', 'linewidth', 1)
axis equal
grid on, axis([-1*a 5*a -1*b 4*b])
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4, 'XTick', xl(1):a:xl(2), 'YTick', yl(1):b:yl(2))
xlabel x
ylabel y


%% #37. Similar to #21 but we can choose if in case of equality we increase i or j
clear
clf
a = 1.35; b = 1;
d_max = 21;
[ii, jj] = ndgrid(0:d_max/min(a,b)); % With this, if there are coincidences only one point is taken: that with larger i 
%[jj, ii] = ndgrid(0:d_max/min(a,b)); % With this, if there are coincidences only one point is taken: that with larger j
ii = ii(:)+2; jj = jj(:)+2;
d2 = ((ii-2)*a).^2 + ((jj-2)*b).^2;
ind = d2<=d_max^2;
ii = ii(ind); jj = jj(ind); d2 = d2(ind);
ij1 = ii+jj-1;
data = [ii jj ij1 d2];
data = sortrows(data, [3 4]);
ind = find(diff([0; data(:, 3)])); % first (minimum) d2 for each ij1
data_chosen = data(ind, :);
%data_chosen
plot(data_chosen(:,1), data_chosen(:,2), 'o'), %axis equal
hold on
xl = xlim;
%plot(xl, a^2/b^2*(xl-2)+(a^2/b^2-1)/2 + 1 + 2) % checked
set(gca, 'colororderindex', 2)
plot(xl, a^2/b^2*xl - 3*a^2/b^2/2 + 5/2, '--') % limiting line (upper): checked
%plot(xl, a^2/b^2*(xl-2) - a^2/b^2/2 -1/2 + 2) % limiting lower line derived for a^2/b^2 odd
set(gca, 'colororderindex', 2)
plot(xl, a^2/b^2*xl - 5*a^2/b^2/2 + 3/2, '-') % checked with previous one
grid on
%figure % to see the "direction" of (i,j) the pairs
%plot(data_chosen(:,3), data_chosen(:,2)./data_chosen(:,1), 'o'), %axis equal
%hold on, plot(xlim, [k k]), ylim(ylim.*[0 1]) % force lower limit 0


%% #37bis. Like #37 but axes are i*a, j*b instead of i, j, and with some adjustments
clear
clf
a = sqrt(2); b = 1;
d_max = 21; % choose manually
[ii, jj] = ndgrid(0:d_max/min(a,b)); % With this, if there are coincidences only one point is taken: that with larger i 
%[jj, ii] = ndgrid(0:d_max/min(a,b)); % With this, if there are coincidences only one point is taken: that with larger j
ii = ii(:)+2; jj = jj(:)+2;
d2 = ((ii-2)*a).^2 + ((jj-2)*b).^2;
ind = d2<=d_max^2;
ii = ii(ind); jj = jj(ind); d2 = d2(ind);
ij1 = ii+jj-1;
data = [ii jj ij1 d2];
data = sortrows(data, [3 4]);
ind = find(diff([0; data(:, 3)])); % first (minimum) d2 for each ij1
data_chosen = data(ind, :);
%data_chosen
plot(data_chosen(:,1)*a, data_chosen(:,2)*b, 'ko', 'markersize', 5)
hold on
xl = xlim;
set(gca, 'colororderindex', 1)
plot(xl, (a^2/b^2*xl/a - 3*a^2/b^2/2 + 5/2)*b, '--') % limiting line, upper
set(gca, 'colororderindex', 1)
plot(xl, (a^2/b^2*xl/a - 5*a^2/b^2/2 + 3/2)*b, '-') % limiting line, lower
grid on
axis equal
xt = 2:10; yt = 0:20; % choose manually
xlim(xt([1 end])*a), ylim(yt([1 end])*b)
xticks(xt*a), yticks(yt*b)
xtl = arrayfun(@(n)[num2str(n) char(8201) 'a'], xt, 'UniformOutput', false); % thin space
xtl(2:2:end) = {''};
xticklabels(xtl)
ytl = arrayfun(@(n)[num2str(n) char(8201) 'b'], yt, 'UniformOutput', false);
ytl(2:2:end) = {''};
yticklabels(ytl)
xlabel(['i' char(8201) 'a'])
ylabel(['j' char(8201) 'b']) % thin space


%% #38. Computation of funl(G). Checked (mostly for a=1 b=1, with #19) 
clear
a = 1; b = 1;
G = 20;

ij = 2+2j;
fprintf('%i: %f\n', 3, hypot((real(ij)-2)*a, (imag(ij)-2)*b))
for g = 4:G
    if imag(ij) <= real(ij)*a^2/b^2 - 3*a^2/b^2/2 + 3/2
        ij = ij + 1j;
    else
        ij = ij + 1;
    end
    fprintf('%i: %f\n', g, hypot((real(ij)-2)*a, (imag(ij)-2)*b))
end


%% #39. Computation of funt(L) (sequential) Checked (a=1, b=1, with #16; and a^2/b^2 integer, with #25)
clc
a = sqrt(2.45); b = 1;
L_step = .001;
L_max = 19;
L = L_step:L_step:L_max;
result = NaN(size(L));
for k = 1:numel(L)
    ij = 2+2j;
    g = 3;
    while L(k) > hypot((real(ij)-2)*a, (imag(ij)-2)*b)
        g = g+1;
        if imag(ij) <= real(ij)*a^2/b^2 - 3*a^2/b^2/2 + 3/2
            ij = ij + 1j;
        else
            ij = ij + 1;
        end
    end
    result(k) = g-1;
end
plot(L, result, '-')
hold on


%% #40. Like #25 but only "odd" case: it always works, even for a^2/b^2 non-integer (comparing with #39):
%clear
%L = .001:.001:19;
%a = sqrt(2.45); b = 1;
L = 3.1;
a = 1.35; b = 1;
k = a^2/b^2;
iip = ceil(real( (k+k^2 + sqrt(-k^3-2*k^2-k + 4*k*(k+1)*L.^2/b^2)) /2/k/(k+1) ))-1;
jjp = ceil( sqrt(L.^2-(iip*a).^2)/b )-1;
S = iip+jjp+3
%plot(L, S, '--')
x = (k+k^2 + sqrt(-k^3-2*k^2-k + 4*k*(k+1)*L.^2/b^2)) /2/k/(k+1) % unquantized solution for i-2
4*a^2*x^2 - 4*a^2*x + a^2+b^2-4*L^2*b^2/(a^2+b^2) % check with its quadratic equation


%% #41. Like #31 but a smaller section is shown, with the lower bound for j.
clf
flag_in = false; % whether the pair is in the minimal sufficient set or not
flag_hexagram = false;
flag_fun = false;
flag_pairs = true;
a = 1.35; b = 1;
hold on
set(gca, 'layer', 'top')
grid on
plot([2*a+2j*b 2*a+5.7j*b], 'k') % "secondary" axes
plot([2*a+2j*b 5.7*a+2j*b], 'k')
chosen = [2*a+2j*b 2*a+3j*b 3*a+3j*b 3*a+4j*b 3*a+5j*b 4*a+5j*b];
chosen_2 = 4*a+4j*b;
if flag_hexagram ms = 4.5; else ms = 5; end
plot(chosen, 'ko', 'markerfacecolor', 'k', 'markersize', ms) % chosen (i,j) pairs
plot(chosen_2, 'ko', 'markerfacecolor', 'w', 'markersize', ms) % chosen 2
plot(setdiff((2:5)*a + 1j*(2:5).'*b, [chosen chosen_2]), 'ko', 'markerfacecolor', 'none', 'markersize', ms) % non-chosen (i,j) pairs
axis equal %pbaspect([1 1 1])
axis([1*a 6*a 1*b 6*b])
theta = linspace(0,pi/2,200);
if flag_in r = 3.1; else r = 3.7; end 
set(gca, 'colororderindex', 1)
z = 2*a+2j*b+r*exp(1j*theta);
z(real(z)> 6.4*a | imag(z)>5.5*b) = [];
plot(z, '-')
plot([1.8 2.4]*a, [2.2 1.6]*b, 'k--')
plot([1.8 3.4]*a, [3.2 1.6]*b, 'k--')
plot([1.8 4.4]*a, [4.2 1.6]*b, 'k--')
plot([1.8 5.4]*a, [5.2 1.6]*b, 'k--')
plot([2.5 5.4]*a, [5.5 2.6]*b, 'k--') % limit to secondary axes
plot([3.5 5.4]*a, [5.5 3.6]*b, 'k--')
plot([4.5 5.4]*a, [5.5 4.6]*b, 'k--')
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4, 'XTick', xl(1):xl(2), 'YTick', yl(1):yl(2))
xticks((0:7)*a), yticks((0:9)*b)
xticklabels([{'0'} arrayfun(@(n)[num2str(n) char(8201) 'a'], 1:7, 'UniformOutput', false)]) % thin space
yticklabels([{'0'} arrayfun(@(n)[num2str(n) char(8201) 'b'], 1:9, 'UniformOutput', false)])
xlabel(['i' char(8201) 'a']), ylabel(['j' char(8201) 'b']) % thin space
set(gca, 'colororderindex', 1)
if flag_fun
    if flag_in text(6.98, 1.43, [char(hex2dec('03C4')) '(' char(hex2dec('2113')) ')' char(hex2dec('2009')) '=' char(hex2dec('2009')) '6']);
    else text(6.98, 2.42, [char(hex2dec('03C4')) '(' char(hex2dec('2113')) ')' char(hex2dec('2009')) '=' char(hex2dec('2009')) '7']);
    end
else
    if flag_in text(6.644, 1.43, ['i+j' char(8211) '1 = 6']);
    else text(6.96, 2.41, ['i+j' char(8211) '1 = 7']);
    end
end
xl = [2.5549 4.6948] + [-.15 .15];
plot(xl*a, ( a^2/b^2*xl - 5*a^2/b^2/2 + 3/2 )*b, '-') % lower-limit line, with axes i*a and i*b
if flag_in text(2.99,5.26,char(hex2dec('2113'))); else text(3.68,5.56,char(hex2dec('2113'))); end
set(gca, 'colororderindex', 1)
if flag_in xl = 3.81568; else xl = 4.08875; end
if flag_hexagram ms = 8; else ms = 6; end
plot(xl*a, ( a^2/b^2*xl - 5*a^2/b^2/2 + 3/2 )*b, 's', 'markersize', ms)
if flag_hexagram
    if flag_in plot(3*a, 4*b, 'kh', 'markersize', 11.5)
    else plot(4*a, 4*b, 'kh', 'markersize', 11.5)
    end
end
if flag_pairs
    if flag_in text(4.5, 3.89, '(i^+,j^+)'), text(4.16-~flag_hexagram*.05, 4.2, '(i^\ast,j^\ast)')
    else text(5.61, 4.42 , '(i^+,j^+)'), text(5.13-~flag_hexagram*.02, 3.76+~flag_hexagram*.03, '(i^\ast,j^\ast)')
    end
end


%% # 42. Inverse problem: funl, for general grids
clear
a = 1.35; b = 1;
k = a^2/b^2;
G = 8;
%ii = (G+5*k/2-1/2)/(1+k); % checked
ii = ((2*G-1)*b^2+5*a^2)/2/(a^2+b^2) % checked
%jj = (G*k + 3/2*(1-k))/(1+k); % checked
jj = (2*G*a^2 + 3*(b^2-a^2))/2/(a^2+b^2); % checked
L = hypot(max(floor(ii)-2,0)*a, (max(ceil(jj)-2,0))*b) % checked


%% #43. Figure for inverse problem, analogous to #41
clf
flag_hexagram = false;
flag_fun = false;
flag_pairs = true;
a = 1.35; b = 1;
hold on
set(gca, 'layer', 'top')
grid on
plot([2*a+2j*b 2*a+5.7j*b], 'k') % "secondary" axes
plot([2*a+2j*b 5.7*a+2j*b], 'k')
chosen = [2*a+2j*b 2*a+3j*b 3*a+3j*b 3*a+4j*b 3*a+5j*b 4*a+5j*b];
chosen_2 = 4*a+4j*b;
if flag_hexagram ms = 4.5; else ms = 5; end
plot(chosen, 'ko', 'markerfacecolor', 'k', 'markersize', ms) % chosen (i,j) pairs
plot(chosen_2, 'ko', 'markerfacecolor', 'w', 'markersize', ms) % chosen 2
plot(setdiff((2:5)*a + 1j*(2:5).'*b, [chosen chosen_2]), 'ko', 'markerfacecolor', 'none', 'markersize', ms) % non-chosen (i,j) pairs
axis equal %pbaspect([1 1 1])
axis([1*a 6*a 1*b 6*b])
plot([1.8 2.4]*a, [2.2 1.6]*b, 'k--')
plot([1.8 3.4]*a, [3.2 1.6]*b, 'k--')
plot([1.8 4.4]*a, [4.2 1.6]*b, 'k--')
plot([1.8 5.4]*a, [5.2 1.6]*b, 'k--')
plot([2.5 5.4]*a, [5.5 2.6]*b, 'k--') % limit to secondary axes
plot([3.5 5.4]*a, [5.5 3.6]*b, 'k--')
plot([4.5 5.4]*a, [5.5 4.6]*b, 'k--')
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4, 'XTick', xl(1):xl(2), 'YTick', yl(1):yl(2))
xticks((0:7)*a), yticks((0:9)*b)
xticklabels([{'0'} arrayfun(@(n)[num2str(n) char(8201) 'a'], 1:7, 'UniformOutput', false)]) % thin space
yticklabels([{'0'} arrayfun(@(n)[num2str(n) char(8201) 'b'], 1:9, 'UniformOutput', false)])
xlabel(['i' char(8201) 'a']), ylabel(['j' char(8201) 'b']) % thin space
set(gca, 'colororderindex', 1)
xl = [2.5549 4.6948] + [-.15 .15];
plot(xl*a, ( a^2/b^2*xl - 5*a^2/b^2/2 + 3/2 )*b, '-') % lower-limit line, with axes i*a and i*b
xl = [2.75 4.15] + [-.15 .15];
set(gca, 'colororderindex', 1)
plot((xl-1)*a, ( a^2/b^2*xl - 5*a^2/b^2/2 + 3/2 )*b+1, '--') % upper-limit line, with axes i*a and i*b
if flag_fun
    ang = 75; quiver(2*a,2*b,(3 -2)*a*.97,(4-2)*b*.97,0,'k')
    text(3.44,3.8,[char(hex2dec('03BB')) '(t)'])
    text(7.1, 1.45, ['t' char(hex2dec('2009')) '=' char(hex2dec('2009')) '6'])
else
    if flag_pairs text(6.78, 1.43, ['i+j' char(8211) '1 = t']), else text(6.68, 1.43, ['i+j' char(8211) '1 = 6']), end
end
set(gca, 'colororderindex', 1)
xl = 3.562887511071745; % computed
if flag_hexagram ms = 8; else ms = 6; end
plot(xl*a, ( a^2/b^2*xl - 5*a^2/b^2/2 + 3/2 )*b, 's', 'markersize', ms)
if flag_hexagram plot(3*a, 4*b, 'kh', 'markersize', 11.5), end
if flag_pairs text(4.94, 3.53, '(i^+,j^+)'), text(4.17-~flag_hexagram*.07, 4.22, '(i_t,j_t)'), end


%% #44. General formulas for funt and funl

% funt:
clear
a = 2.35^3; b = 1;
L = .001:.001:50;
ii = ceil((a + b*real(sqrt(4*L.^2/(a^2+b^2)-1))) / 2/a) + 1;
jj = ceil(sqrt( L.^2 - (ii-2).^2*a^2 )/b) + 1;
funt = ii+jj-1;
plot(L, funt)
xlabel L, ylabel G, hold on, grid on

% funl
G = 3:max(funt);
ii = ((2*G-1)*b^2+5*a^2)/2/(a^2+b^2);
jj = (2*G*a^2 + 3*(b^2-a^2))/2/(a^2+b^2);
funl = hypot(max(floor(ii)-2,0)*a, (max(ceil(jj)-2,0))*b);
plot(funl, G, '.') % checked.
% Also checked that funl is the same if a and b are swapped

% Approx to funt
plot(L, L*hypot(1/a,1/b)+5/2, '--')

% Approximation error
figure
approx_error = funt-(L*hypot(1/a,1/b)+5/2);
plot(L, approx_error)
[min(approx_error) max(approx_error)]
ind = ~mod(L,1); % integer lengths
hold on
plot(L(ind), approx_error(ind), '.')
plot(xlim, -.5*[1 1], 'k:'), plot(xlim, .5*[1 1], 'k:')
% The error does not seem to be bounded between -.5 and .5, even for integer L only and discarding the first
% values


%% #45. Examples of funt and of funl
clear
ft = figure; hold on
fl = figure; hold on
L = .001:.001:35;
G = 1:35;
a_all = [ 1  1.35   5   5  10 ];
b_all = [ 1     1   1 1.5  3 ];
for k = 1:numel(a_all)
    figure(ft)
    a = a_all(k); b = b_all(k);
    ii = ceil((a + b*real(sqrt(4*L.^2/(a^2+b^2)-1))) / 2/a) + 1;
    jj = ceil(sqrt( L.^2 - (ii-2).^2*a^2 )/b) + 1;
    funt = ii+jj-1;
    ind = funt<=max(G);
    %set(gca, 'colororderindex', 1)
    plot(L(ind), funt(ind))
    
    figure(fl)
    ii = ((2*G-1)*b^2+5*a^2)/2/(a^2+b^2);
    jj = (2*G*a^2 + 3*(b^2-a^2))/2/(a^2+b^2);
    funl = hypot(max(floor(ii)-2,0)*a, (max(ceil(jj)-2,0))*b);
    ind = funl<=max(L);
    %set(gca, 'colororderindex', 1)
    plot(G(ind), funl(ind), 'o', 'markersize', 4)
    
end
figure(ft), xlabel L, ylabel G, hold on, grid on, axis equal, axis([0 max(L) 0 max(G)]) 
figure(fl), xlabel G, ylabel L, hold on, grid on, axis equal, axis([0 max(L) 0 max(G)])


%% #46. Check formulas for square grid with real-valued lengths
clear
clc
L = .001:.001:35;
G = 1:35;
a = 1.5754;
b = a;

ii = ceil((a + b*real(sqrt(4*L.^2/(a^2+b^2)-1))) / 2/a) + 1;
jj = ceil(sqrt( L.^2 - (ii-2).^2*a^2 )/b) + 1;
funt_rect = ii+jj-1;

ii = ceil((1 + real(sqrt(2*L.^2/a^2-1)))/2) + 1;
jj = ceil(sqrt(L.^2/a^2 - (ii-2).^2 )) + 1;
funt_sq = ii+jj-1;

ii = ceil(L/a/sqrt(2)) + 1;
jj = ceil(sqrt(L.^2/a^2 - (ii-2).^2 )) + 1;
funt_sq_2 = ii+jj-1;
isequal(funt_rect, funt_sq, funt_sq_2)

ii = ((2*G-1)*b^2+5*a^2)/2/(a^2+b^2);
jj = (2*G*a^2 + 3*(b^2-a^2))/2/(a^2+b^2);
ii = round(2*ii)/2; jj = round(2*jj)/2; % fix numerical issues
funl_rect = hypot(max(floor(ii)-2,0)*a, (max(ceil(jj)-2,0))*b);

ind = G<3;
funl_sq(ind) = 0;
ind = mod(G,2) & (G>=3);
funl_sq(ind) = (G(ind)-3)/sqrt(2) * a;
ind = ~mod(G,2) & (G>=3);
funl_sq(ind) = sqrt(((G(ind)-3).^2+1)/2) * a;
max(abs(funl_rect-funl_sq))


%% #47. Check formulas for unit square grid with integer-valued lengths

clear, clc

L = 1:1e7;
ii = ceil(L/sqrt(2)) + 1;
jj = ceil(sqrt(L.^2 - (ii-2).^2 )) + 1;
funt_real = ii+jj-1; % formula valid for real-valued lengths

funt_int = floor(sqrt(2*L.^2-2))+3;
isequal(funt_real, funt_int)

G = 1:1e7;
ind = G<3;
funl_real(ind) = 0;
ind = mod(G,2) & (G>=3);
funl_real(ind) = (G(ind)-3)/sqrt(2);
ind = ~mod(G,2) & (G>=3);
funl_real(ind) = sqrt(((G(ind)-3).^2+1)/2);
funl_real = floor(funl_real)+1; % formula obtained from the case with real-valued lengths

ind = G<3;
funl_int(ind) = 1;
funl_int(~ind) = ceil(sqrt((G(~ind)-3).^2/2+1));
isequal(funl_real, funl_int)


%%  #48. Figures of the direct-problem and inverse problem sequences (funti and funli)
clear
close all
figure
L = 1:25;
funti = floor(sqrt(2*L.^2-2))+3;
stem(L, funti, 'o')
grid on, % axis equal
axis([0 25 0 40])
xlabel(char(hex2dec('2113'))), ylabel(['T(' char(hex2dec('2113')) ')'])

clear
figure
G = 1:25;
ind = G<3;
funli(ind) = 1;
funli(~ind) =  ceil(sqrt((G(~ind)-3).^2/2+1));
stem(G, funli, 'o')
grid on, % axis equal
axis([0 25 0 16])
xlabel t, ylabel([char(hex2dec('039B')) '(t)'])


%% #49. Average number of touched squares, experimentally

clear, clc
L = .1:.1:9;
a = 1.35; b = 1;
N = 1e5;
z0 = rand(N, numel(L))*a + 1j*rand(N, numel(L))*b;
z1 = z0 + L.*exp(1j*2*pi*rand(N, numel(L)));
ii = ceil(abs(real(z1/a))) + (real(z1)<0); % add 1 for negative
jj = ceil(abs(imag(z1/b))) + (imag(z1)<0); % add 1 for negative
result = mean(ii + jj, 1) - 1;
plot(L, result, '.-'), grid on, axis equal


%% #50. Selectable range for theta
% I see that theta on (0,pi/2) doesn't give the same result. Plottin with something like
% plot([z0(1:600,60), z1(1:600,60)].', '.-'), samples seem to be missing near the extremes

L = .1:.1:9;
a = 1.35; b = 1;
N = 1e5;
theta_max = 2*pi;   
unfix = @(x) sign(x).*ceil(abs(x)); % round away from 0; real values only
z0 = rand(N, numel(L))*a + 1j*rand(N, numel(L))*b;
z1 = z0 + L.*exp(1j*theta_max*rand(N, numel(L)));
ii = unfix(real(z1/a)) - (real(z1)<0); % subtract 1 for negative; it's signed
jj = unfix(imag(z1/b)) - (imag(z1)<0); % subtract 1 for negative; it's signed
result = mean(abs(ii) + abs(jj), 1) - 1;
plot(L, result, '.-'), grid on, axis equal


%% #51. I define ki, kj signed (page 17, (*)), from which ii, jj are obtained 

L = .1:.1:5;
a = 1.35; b = 1;
N = 3e5;%1e5;
theta_max = 2*pi;
%z0 = repmat(.42+.73j, N, numel(L));
z0 = rand(N, numel(L))*a + 1j*rand(N, numel(L))*b;
z1 = z0 + L.*exp(1j*theta_max*rand(N, numel(L)));
ki = floor(real(z1/a));
kj = floor(imag(z1/b));
ii = 1+abs(ki);
jj = 1+abs(kj);
result = mean(ii + jj, 1) - 1;
plot(L, result, '.-'), grid on, axis equal

ind_L_test = 25;
ki_test = 2;
mean(ki(:,ind_L_test)>=ki_test)
mean(real(acos((a*ki_test-real(z0(:,ind_L_test)))/L(ind_L_test))/pi)) % checked
ki_test = -2;
mean(ki(:,ind_L_test)<ki_test+1)
mean(real(acos((a*(abs(ki_test)-1)+real(z0(:,ind_L_test)))/L(ind_L_test))/pi)) % checked

ind_L_test = 25;
kj_test = 2;
mean(kj(:,ind_L_test)>=kj_test)
mean(real(acos((b*kj_test-imag(z0(:,ind_L_test)))/L(ind_L_test))/pi)) % checked
kj_test = -2;
mean(kj(:,ind_L_test)<kj_test+1)
mean(real(acos((b*(abs(kj_test)-1)+imag(z0(:,ind_L_test)))/L(ind_L_test))/pi)) % checked

ind_L_test = 65;
ii_test = 2;
mean(ii(:,ind_L_test)>=ii_test)
mean(real(acos((a*(ii_test-1)-real(z0(:,ind_L_test)))/L(ind_L_test))/pi)) + ...
    mean(real(acos((a*(ii_test-2)+real(z0(:,ind_L_test)))/L(ind_L_test))/pi)) % checked
f = @(t,len) -t.*acos(t./len) + sqrt(len.^2-t.^2);
2/pi/a*(f(a*ii_test-2*a, L(ind_L_test)) - f(a*ii_test-a, L(ind_L_test))) % checked
g = @(t) t.*acos(t) - sqrt(1-t.^2);
2/pi/a*L(ind_L_test)*(g(a*(ii_test-1)/L(ind_L_test)) - g(a*(ii_test-2)/L(ind_L_test))) % checked

ind_L_test = 65;
ii_test = 4; % valid for ii_test >= 2
mean(ii(:,ind_L_test)==ii_test)
mean(real(acos((a*(ii_test-1)-real(z0(:,ind_L_test)))/L(ind_L_test))/pi)) - ...
    mean(real(acos((a*ii_test-real(z0(:,ind_L_test)))/L(ind_L_test))/pi)) + ...
    mean(real(acos((a*(ii_test-2)+real(z0(:,ind_L_test)))/L(ind_L_test))/pi)) - ...
    mean(real(acos((a*(ii_test-1)+real(z0(:,ind_L_test)))/L(ind_L_test))/pi)) % checked
mean(ii(:,ind_L_test)==1)
1 - mean(real(acos((a-real(z0(:,ind_L_test)))/L(ind_L_test))/pi)) - ...
    mean(real(acos((real(z0(:,ind_L_test)))/L(ind_L_test))/pi)) % checked

mean(ii(:,ind_L_test))
2/pi/a*L(ind_L_test)+1 % checked!

jj_test = 4;
mean(jj(:,ind_L_test))
2/pi/b*L(ind_L_test)+1

ind_L_test = 65;
mean(ii(:,ind_L_test)+jj(:,ind_L_test)-1)
2*L(ind_L_test)/pi*(1/a+1/b)+1 % checked!!


%% #52. Plot of formula 2*L/pi*(1/a+1/b)+1

L = .2:.2:10;
%a = 1.35; b = 2;
%a = 4/pi; b = a;
a = 1; b = 1000;
N = 1e5;
theta_max = 2*pi;
%z0 = repmat(.42+.73j, N, numel(L));
z0 = rand(N, numel(L))*a + 1j*rand(N, numel(L))*b;
z1 = z0 + L.*exp(1j*theta_max*rand(N, numel(L)));
ki = floor(real(z1/a));
kj = floor(imag(z1/b));
ii = 1+abs(ki);
jj = 1+abs(kj);
result = mean(ii + jj, 1) - 1;
plot(L, result, 'o'), grid on, axis equal
hold on
plot(L, 2*L/pi*(1/a+1/b)+1)


%% #53. Figures for proof of Pr[i>=n]

clear, clf
a = 1.35; b = 1;
x1 = .6*a; y1 = .69*b;
n = 3;
min_x = -2; max_x = 3; min_y = -3; max_y = 4;
flag_full = false;
flag_segment = true;
if flag_full
    len = 2.1*a;
else
    len = 1.78*a;
end
m = max(0, a*(n-1)-len);
patch([m a a m m], [0 0 b b 0], .88*[1 1 1], 'edgecolor', 'none')
hold on
if ~flag_full
    plot([m m+1j*b], 'k:')
end
patch([0 a a 0 0], [0 0 b b 0], 1*[1 1 1], 'edgecolor', 'k', 'FaceAlpha', 0)
plot(x1 + 1j*y1, 'k.', 'markersize', 7)
grid on, axis equal
set(gca, 'GridAlpha', .4)
xticks((min_x:max_x)*a), yticks((min_y:max_y)*b)
xticklabels([arrayfun(@(n)[num2str(n) char(8201) repmat('a', 1, n~=0)], min_x:max_x, 'UniformOutput', false)]) % thin space
yticklabels([arrayfun(@(n)[num2str(n) char(8201) repmat('b', 1, n~=0)], min_y:max_y, 'UniformOutput', false)])
xlabel x, ylabel y
axis([min_x*a max_x*a min_y*b max_y*b])
set(gca, 'layer', 'top')
plot([min_y*b*1j max_y*b*1j], 'k')
plot(complex([min_x*a max_x*a]), 'k')
plot(a*(n-1)+[(min_y*b+.7)*1j max_y*b*1j], 'k--', 'linewidth', .75)
plot(-a*(n-2)+[(min_y*b+.7)*1j max_y*b*1j], 'k--', 'linewidth', .75)
theta_range = [-63 65]/180*pi;
theta = linspace(theta_range(1), theta_range(2), 1000);
z = x1+1j*y1 + len*exp(1j*theta);
set(gca, 'colororderindex', 1)
plot(z, '--', 'linewidth', .75)
ind = real(z)>=a*(n-1);
set(gca, 'colororderindex', 1)
plot(z(ind), '-', 'linewidth', .75)
if flag_segment
    if flag_full
        ang_segment = 42.5/180*pi;
    else
        ang_segment = 19/180*pi;
    end
    plot([x1+1j*y1 x1+1j*y1+len*exp(1j*ang_segment)],'k')
    plot(x1+1j*y1+len*exp(1j*ang_segment), 'k.', 'markersize', 7)
else
    ang_quiver = 59.5/180*pi;
    quiver(x1,y1,len*cos(ang_quiver),len*sin(ang_quiver),.99,'k')
end
if flag_full
    text(0.44, 0.43, 'x_1,y_1')
    if flag_segment
        text(2.0, 2.28, char(hex2dec('2113')))
        text(3.13, 2.75, 'x_2,y_2')
    else
        text(1.52, 2.63, char(hex2dec('2113')))
    end
else
    text(0.52, 0.45, 'x_1,y_1')
    if flag_segment
        text(1.94, 1.41, char(hex2dec('2113')))
        text(3.2, 1.73, 'x_2,y_2')
    else
        text(1.43, 2.44, char(hex2dec('2113')))
    end
end
text(2.35, -2.56, ['x=a(n' char(8211) '1)'])
text(-1.87, -2.56, ['x=' char(8211) 'a(n' char(8211) '2)'])


%% #54. Plot of average number of tiles

clear
clf, hold on
L = .001:.001:35;
G = 1:35;
a_all = [ 1.35 ];
b_all = [ 1  ];
for k = 1:numel(a_all)
    a = a_all(k); b = b_all(k);
    ii = ceil((a + b*real(sqrt(4*L.^2/(a^2+b^2)-1))) / 2/a) + 1;
    jj = ceil(sqrt( L.^2 - (ii-2).^2*a^2 )/b) + 1;
    funt = ii+jj-1;
    ind = funt<=max(G);
    set(gca, 'colororderindex', k)
    plot(L(ind), funt(ind), '-')
    set(gca, 'colororderindex', k)
    plot(L, 2*L/pi*(1/a+1/b)+1, '-')
    set(gca, 'colororderindex', k)
    plot(L, ceil(sqrt(max(L.^2-b^2, 1)))/a, '-')
end
xlabel(char(hex2dec('2113'))), ylabel([char(hex2dec('03C4')) '(' char(hex2dec('2113')) ')']), grid on, axis equal, axis([0 max(L) 0 max(G)]) 


%% #54. Plot of ratio of average to maximum number of tiles

clear
clf, hold on
L = .001:.001:50;
G = 1:35;
a_all = [ 1 1.35 ];
b_all = [ 1 1  ];
for k = 1:numel(a_all)
    a = a_all(k); b = b_all(k);
    ii = ceil((a + b*real(sqrt(4*L.^2/(a^2+b^2)-1))) / 2/a) + 1;
    jj = ceil(sqrt( L.^2 - (ii-2).^2*a^2 )/b) + 1;
    funt = ii+jj-1;
    funta = 2*L/pi*(1/a+1/b)+1;
    set(gca, 'colororderindex', k)
    plot(L, funt./funta, '-')
end
xlabel(char(hex2dec('2113'))), ylabel([char(hex2dec('03C4')) '(' char(hex2dec('2113')) ')']), grid on


%% #55. Like #49 but for arbitrary number of dimensions: average number of visited cells, experimentally

clear, clc
L = .1:.1:9;
a = [1.35 1 .8];
N = 1e4;
ndim = numel(a);
result = NaN(size(L));
for iter_L = 1:numel(L)
    dir = 2*rand(N, ndim)-1;
    ind = sum(dir.^2, 2) <= 1;
    dir = dir(ind,:);
    dir = dir./sqrt(sum(dir.^2, 2));
    %plot3(dir(:,1), dir(:,2), dir(:,3), '.', 'markersize', 2), axis equal
    z0 = rand(size(dir)).*a(:).';
    z1 = z0 + L(iter_L)*dir;
    ii = ceil(abs(z1./a)) + (z1<0); % add 1 for negative
    result(iter_L) = mean(sum(ii, 2)) - ndim + 1;
end
plot(L, result, 'o'), grid on, axis equal


%% #56. Like #55 but cells are counted using deterministic uniform sampling along the segment
% Checked: gives simular results as #55

clear, clc
L = .1:.3:5;
a = [1.35 1 .8];
N = 3e3;
n_points_along = 1000;
ndim = numel(a);
t = linspace(0, 1, n_points_along);
result = NaN(size(L));
for iter_L = 1:numel(L)
    dir = 2*rand(N, ndim)-1;
    ind = sum(dir.^2, 2) <= 1;
    dir = dir(ind,:);
    dir = dir./sqrt(sum(dir.^2, 2));
    %plot3(dir(:,1), dir(:,2), dir(:,3), '.', 'markersize', 2), axis equal
    z0 = rand(size(dir)).*a;
    z1 = z0 + L(iter_L)*dir;
    size_unique  = 0;
    for n = 1:size(z0,1)
        size_unique = size_unique + ...
            size(unique(floor((z0(n,:).*(1-t(:)) + z1(n,:).*t(:))./a), 'rows'), 1);
    end
    result(iter_L) = size_unique/size(z0,1);
end
plot(L, result, '.'), grid on, axis equal, hold on
switch ndim
case 2
    plot(L, 1+2*L/pi*sum(1./a), '-')
case 3
    plot(L, 1+L/2*sum(1./a), '-')
end
% There's a general formula for the coefficient (Jacuwes, 18; found by Alex)


%% #57. Probability of visiting the maximum number of tiles, experimentally. Based on #47 and #49

clear, clc
L = 3.8:.1:4.2; %3.8:.1:4.2;%3:.1:3.5;
N = 1e7;
ii = ceil(L/sqrt(2)) + 1;
jj = ceil(sqrt(L.^2 - (ii-2).^2 )) + 1;
funt_real = ii+jj-1 % formula valid for real-valued lengths
z0 = rand(N, numel(L)) + 1j*rand(N, numel(L));
z1 = z0 + L.*exp(1j*2*pi*rand(N, numel(L)));
ii = ceil(abs(real(z1))) + (real(z1)<0); % add 1 for negative
jj = ceil(abs(imag(z1))) + (imag(z1)<0); % add 1 for negative
result = mean(ii+jj-1 == funt_real , 1)


%% #58. Computation for #57, for odd values of the maximum number of tiles. Page 21. 
% This should be run after #57. funt_real must contain odd values
% Checked with #57

clc
result_computed = NaN(size(result));
for n = 1:numel(L)
    k = funt_real(n)-1;
    assert(mod(k,2)==0)
    %alpha = acos((k/2-1)/L(n))-pi/4;
    %result_computed(n) = (alpha*(k/2-1)^2 - L(n)^2/4*(cos(2*(pi/4+alpha))-cos(2*pi/4)) - L(n)*(k/2-1)*(-cos(pi/4+alpha)+sin(pi/4+alpha)+cos(pi/4)-sin(pi/4))) * 4/pi;
    %result_computed(n) = integral(@(theta) (1-k/2+L(n)*sin(theta)).*(1-k/2+L(n)*cos(theta)), pi/4, pi/4+alpha) * 4/pi;
    alphap = acos((k/2-1)/L(n));
    %result_computed(n) = ((alphap-pi/4)*(k/2-1)^2 - L(n)^2/4*cos(2*alphap) + L(n)*(k/2-1)*(cos(alphap)-sin(alphap))) * 4/pi;
    %result_computed(n) = ((alphap-pi/4)*(k/2-1)^2 - L(n)^2/4*(2*((k-2)/2/L(n))^2-1) + L(n)*(k-2)/2*((k-2)/2/L(n)-sqrt(1-((k-2)/2/L(n))^2))) * 4/pi;
    %result_computed(n) = ((alphap-pi/4)*(k-2)^2/4 - ((k-2)^2-2*L(n)^2)/8 + (k-2)/4*(k-2-sqrt(4*L(n)^2-(k-2)^2))) * 4/pi;
    result_computed(n) = ((alphap-pi/4)*(k-2)^2/4 + (k-2)^2/8 + L(n)^2/4 - (k-2)/4*sqrt(4*L(n)^2-(k-2)^2)) * 4/pi;
    if L(n)^2 > (k/2-2)^2+(k/2)^2 % add neighbours
        %result_computed(n) = result_computed(n) + integral(@(theta) (2-k/2+L(n)*sin(theta)).*(-k/2+L(n)*cos(theta)),  asin((k/2-2)/L(n)), acos(k/2/L(n))) * 4/pi;
        alphap = acos(k/2/L(n)); betap =  asin((k/2-2)/L(n));
        %result_computed(n) = result_computed(n) + ( (alphap-betap)*k*(k-4)/4 - L(n)^2/4*( 2*(k/2/L(n))^2 - 2 + 2*((k-4)/2/L(n))^2) + L(n)*k/2*(k/2/L(n)-sqrt(1-((k-4)/2/L(n))^2)) - L(n)*(k-4)/2*(sqrt(1-(k/2/L(n))^2)-(k-4)/2/L(n)) ) * 4/pi;
        %result_computed(n) = result_computed(n) + ( (alphap-betap)*k*(k-4)/4 - (k^2+(k-4)^2-4*L(n)^2)/8 + k/4*(k-sqrt(4*L(n)^2-(k-4)^2)) +(k-4)/4*(k-4-sqrt(4*L(n)^2-k^2)) ) * 4/pi;
        result_computed(n) = result_computed(n) + ( (alphap-betap)*k*(k-4)/4 + k^2/8 + (k-4)^2/8 + L(n)^2/2 - k/4*sqrt(4*L(n)^2-(k-4)^2) - (k-4)/4*sqrt(4*L(n)^2-k^2) ) * 4/pi;
    end
end
[result; result_computed].'


%% #59. Computation for #57, for even values of the maximum number of tiles. Page 23. 
% This should be run after #57. funt_real must contain even values

result_computed = NaN(size(result));
for n = 1:numel(L)
    k = funt_real(n)-1;
    assert(mod(k,2)==1)
    alphap = acos((k/2-1/2)/L(n));
    betap =  asin((k/2-3/2)/L(n));
    %result_computed(n) = ((alphap-betap)*(k-1)*(k-3)/4 - L(n)^2/4*(cos(2*alphap)-cos(2*betap)) + L(n)*(k-1)/2*(cos(alphap)-cos(betap)) - L(n)*(k-3)/2*(sin(alphap)-sin(betap)) ) * 4/pi;
    result_computed(n) = ( (alphap-betap)*(k-1)*(k-3)/4 + (k-1)^2/8 + (k-3)^2/8 + L(n)^2/2 - (k-1)/4*sqrt(4*L(n)^2-(k-3)^2) - (k-3)/4*sqrt(4*L(n)^2-(k-1)^2) ) *4/pi;
    if L(n)^2 > (k/2-5/2)^2+(k/2+1/2)^2 % add neighbours
        alphap = acos((k+1)/2/L(n)); betap =  asin((k-5)/2/L(n));
        result_computed(n) = result_computed(n) + ( (alphap-betap)*(k+1)*(k-5)/4 + (k+1)^2/8 + (k-5)^2/8 + L(n)^2/2 - (k+1)/4*sqrt(4*L(n)^2-(k-5)^2) - (k-5)/4*sqrt(4*L(n)^2-(k+1)^2) ) * 4/pi;
    end
end
[result; result_computed].'


%% #60. All together now

clear, clc
L = .1:.1:10; %3.8:.1:4.2;%3:.1:3.5;
N = 1e6;
ii = ceil(L/sqrt(2)) + 1;
jj = ceil(sqrt(L.^2 - (ii-2).^2 )) + 1;
funt = ii+jj-1 % formula valid for real-valued lengths
z0 = rand(N, numel(L)) + 1j*rand(N, numel(L));
z1 = z0 + L.*exp(1j*2*pi*rand(N, numel(L)));
ii = ceil(abs(real(z1))) + (real(z1)<0); % add 1 for negative
jj = ceil(abs(imag(z1))) + (imag(z1)<0); % add 1 for negative
result = mean(ii+jj-1 == funt , 1);
plot(L, result, 'o', 'markersize', 4), hold on, grid on

L_fine = linspace(.01,10,1e4);
ii = ceil(L_fine/sqrt(2)) + 1;
jj = ceil(sqrt(L_fine.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
ind_odd = mod(funt_fine,2)==1; % odd maximum number of squares
ind_odd2 = ind_odd & (2*L_fine.^2>funt_fine.^2-6*funt_fine+13);
alpha_odd = acos((funt_fine(ind_odd)-3)/2./L_fine(ind_odd));
alphap_odd2 = acos((funt_fine(ind_odd2)-1)/2./L_fine(ind_odd2));
betap_odd2 = asin((funt_fine(ind_odd2)-5)/2./L_fine(ind_odd2));
result_computed(ind_odd) = ( (alpha_odd-pi/4).*(funt_fine(ind_odd)-3).^2 + (funt_fine(ind_odd)-3).^2/2 + ...
    L_fine(ind_odd).^2 - (funt_fine(ind_odd)-3).*sqrt(4*L_fine(ind_odd).^2-(funt_fine(ind_odd)-3).^2) ) /pi;
result_computed(ind_odd2) = result_computed(ind_odd2) + ...
    ( (alphap_odd2-betap_odd2).*(funt_fine(ind_odd2)-1).*(funt_fine(ind_odd2)-5) + (funt_fine(ind_odd2)-1).^2/2 + (funt_fine(ind_odd2)-5).^2/2 + ...
    2*L_fine(ind_odd2).^2 - (funt_fine(ind_odd2)-1).*sqrt(4*L_fine(ind_odd2).^2 - (funt_fine(ind_odd2)-5).^2) - (funt_fine(ind_odd2)-5).*sqrt(4*L_fine(ind_odd2).^2 - (funt_fine(ind_odd2)-1).^2) ) / pi;
ind_even = mod(funt_fine,2)==0; % even maximum number of squares
ind_even2 = ind_even & (2*L_fine.^2>funt_fine.^2-6*funt_fine+18);
alpha_even = acos((funt_fine(ind_even)-2)/2./L_fine(ind_even));
beta_even = asin((funt_fine(ind_even)-4)/2./L_fine(ind_even));
alphap_even2 = acos(funt_fine(ind_even2)/2./L_fine(ind_even2));
betap_even2 = asin((funt_fine(ind_even2)-6)/2./L_fine(ind_even2));
result_computed(ind_even) = ( (alpha_even-beta_even).*(funt_fine(ind_even)-2).*(funt_fine(ind_even)-4) + (funt_fine(ind_even)-2).^2/2 + (funt_fine(ind_even)-4).^2/2 + ...
    2*L_fine(ind_even).^2 - (funt_fine(ind_even)-2).*sqrt(4*L_fine(ind_even).^2 - (funt_fine(ind_even)-4).^2) - (funt_fine(ind_even)-4).*sqrt(4*L_fine(ind_even).^2 - (funt_fine(ind_even)-2).^2) ) / pi;
result_computed(ind_even2) = result_computed(ind_even2) + ...
    ( (alphap_even2-betap_even2).*funt_fine(ind_even2).*(funt_fine(ind_even2)-6) + funt_fine(ind_even2).^2/2 + (funt_fine(ind_even2)-6).^2/2 + ...
    2*L_fine(ind_even2).^2 - funt_fine(ind_even2).*sqrt(4*L_fine(ind_even2).^2 - (funt_fine(ind_even2)-6).^2) - (funt_fine(ind_even2)-6).*sqrt(4*L_fine(ind_even2).^2 - funt_fine(ind_even2).^2) ) / pi;
%plot(L_fine(ind_odd), result_computed(ind_odd), '.', 'markersize', 3), hold on, grid on
%plot(L_fine(ind_even), result_computed(ind_even), '.', 'markersize', 3)
plot(L_fine.*ind_odd./ind_odd, result_computed.*ind_odd./ind_odd, '-', 'linewidth', .7), hold on, grid on
plot(L_fine.*ind_even./ind_even, result_computed.*ind_even./ind_even, '-', 'linewidth', .7)
legend({'Simulated' 'Computed, odd max number of squares' 'Computed, even max number of squares'})


%% #61. Probabilities at the end of the intervals with odd maximum number of segments

t = 3:2:7001; L_fine = sqrt((t.^2-4*t+5)/2)-1e-10;
ii = ceil(L_fine/sqrt(2)) + 1;
jj = ceil(sqrt(L_fine.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
ind_odd = mod(funt_fine,2)==1; % odd maximum number of squares
ind_odd2 = ind_odd & (2*L_fine.^2>funt_fine.^2-6*funt_fine+13);
alpha_odd = acos((funt_fine(ind_odd)-3)/2./L_fine(ind_odd));
alphap_odd2 = acos((funt_fine(ind_odd2)-1)/2./L_fine(ind_odd2));
betap_odd2 = asin((funt_fine(ind_odd2)-5)/2./L_fine(ind_odd2));
result_computed(ind_odd) = ( (alpha_odd-pi/4).*(funt_fine(ind_odd)-3).^2 + (funt_fine(ind_odd)-3).^2/2 + ...
    L_fine(ind_odd).^2 - (funt_fine(ind_odd)-3).*sqrt(4*L_fine(ind_odd).^2-(funt_fine(ind_odd)-3).^2) ) /pi;
result_computed(ind_odd2) = result_computed(ind_odd2) + ...
    ( (alphap_odd2-betap_odd2).*(funt_fine(ind_odd2)-1).*(funt_fine(ind_odd2)-5) + (funt_fine(ind_odd2)-1).^2/2 + (funt_fine(ind_odd2)-5).^2/2 + ...
    2*L_fine(ind_odd2).^2 - (funt_fine(ind_odd2)-1).*sqrt(4*L_fine(ind_odd2).^2 - (funt_fine(ind_odd2)-5).^2) - (funt_fine(ind_odd2)-5).*sqrt(4*L_fine(ind_odd2).^2 - (funt_fine(ind_odd2)-1).^2) ) / pi;
plot(L_fine, result_computed.*funt_fine), grid on


%% #61. Probabilities at the end of the intervals with even maximum number of segments

t = 4:2:7002; L_fine = (t-2)/sqrt(2)-1e-10;
ii = ceil(L_fine/sqrt(2)) + 1;
jj = ceil(sqrt(L_fine.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
ind_even = mod(funt_fine,2)==0; % even maximum number of squares
ind_even2 = ind_even & (2*L_fine.^2>funt_fine.^2-6*funt_fine+18);
alpha_even = acos((funt_fine(ind_even)-2)/2./L_fine(ind_even));
beta_even = asin((funt_fine(ind_even)-4)/2./L_fine(ind_even));
alphap_even2 = acos(funt_fine(ind_even2)/2./L_fine(ind_even2));
betap_even2 = asin((funt_fine(ind_even2)-6)/2./L_fine(ind_even2));
result_computed(ind_even) = ( (alpha_even-beta_even).*(funt_fine(ind_even)-2).*(funt_fine(ind_even)-4) + (funt_fine(ind_even)-2).^2/2 + (funt_fine(ind_even)-4).^2/2 + ...
    2*L_fine(ind_even).^2 - (funt_fine(ind_even)-2).*sqrt(4*L_fine(ind_even).^2 - (funt_fine(ind_even)-4).^2) - (funt_fine(ind_even)-4).*sqrt(4*L_fine(ind_even).^2 - (funt_fine(ind_even)-2).^2) ) / pi;
result_computed(ind_even2) = result_computed(ind_even2) + ...
    ( (alphap_even2-betap_even2).*funt_fine(ind_even2).*(funt_fine(ind_even2)-6) + funt_fine(ind_even2).^2/2 + (funt_fine(ind_even2)-6).^2/2 + ...
    2*L_fine(ind_even2).^2 - funt_fine(ind_even2).*sqrt(4*L_fine(ind_even2).^2 - (funt_fine(ind_even2)-6).^2) - (funt_fine(ind_even2)-6).*sqrt(4*L_fine(ind_even2).^2 - funt_fine(ind_even2).^2) ) / pi;
plot(L_fine, result_computed.*funt_fine), grid on


%% #62. Figure for proof of probmax: possible tiles
clf
hold on
set(gca, 'layer', 'top')
axis equal
grid on
plot([-3.7 4.7], [0 0], 'k')
plot([-3.7j 4.7j], 'k')
set(gca, 'GridAlpha', .4)
xticks(-3:4), yticks(-3:4)
%xticklabels([arrayfun(@(n)[num2str(n) char(8201) 'a'], -3:-1, 'UniformOutput', false) {'0'}  arrayfun(@(n)[num2str(n) char(8201) 'a'], 1:4, 'UniformOutput', false)]) % thin space
%yticklabels([arrayfun(@(n)[num2str(n) char(8201) 'a'], -3:-1, 'UniformOutput', false) {'0'}  arrayfun(@(n)[num2str(n) char(8201) 'a'], 1:4, 'UniformOutput', false)]) % thin space
xlabel x 
xlabel y
axis([-4 5 -4 5])
for k = [0 1+3j 2+2j 3+1j 1-3j 2-2j 3-1j -3+1j -2+2j -1+3j -3-1j -2-2j -1-3j]
    if k==0
        plot(real(k)+[0 1 1 0 0], imag(k)+[0 0 1 1 0], 'k-', 'linewidth', 1)
    else
        plot(real(k)+[0 1 1 0 0], imag(k)+[0 0 1 1 0], 'k--', 'linewidth', 1.4)
    end
end


%% #63. Figure for proof of probmax: funt odd, diagonal
clf
hold on
set(gca, 'layer', 'top')
axis equal
grid on
plot([-.4 5], [0 0], 'k')
plot([-.4j 5j], 'k')
set(gca, 'GridAlpha', .4)
xticks(-1:5), yticks(-1:5)
xticklabels({'' '0' '1' '' ['(t' char(hex2dec('2212')) '1)/2'] '' ''})
yticklabels({'' '0' '1' '' ['(t' char(hex2dec('2212')) '1)/2'] '' ''})
xlabel x 
ylabel y
axis([-1 5 -1 5])
L = 3.37;
theta_0 = asin(2/L);
theta = 3*pi/2 - linspace(30*pi/180,theta_0,50);
plot(3+3j + L*exp(1j*theta), '--', 'linewidth', .75)
theta_1 = acos(2/L);
theta = 3*pi/2 - linspace(theta_1,60*pi/180,50);
set(gca, 'colororderindex', 1)
plot(3+3j + L*exp(1j*theta), '--', 'linewidth', .75)
theta = 3*pi/2 - linspace(theta_0,theta_1,50);
set(gca, 'colororderindex', 1)
plot(3+3j + L*exp(1j*theta), '-', 'linewidth', .75)
plot([3+3j], 'k.', 'markersize', 7)
theta_example = theta_0 + (theta_1-theta_0)*.67;
endpoint = 3+3j + L*exp(1j*(3*pi/2 - theta_example));
plot(endpoint, 'k.', 'markersize', 7)
plot([endpoint real(endpoint)+1j], 'k:')
plot([endpoint 1+1j*imag(endpoint)], 'k:')
patch([real(endpoint) 1 1 real(endpoint) real(endpoint)], [imag(endpoint) imag(endpoint) 1 1 imag(endpoint)], .88*[1 1 1], 'edgecolor', 'none')
plot([3+3j endpoint], 'k')
for k = [0 3+3j]
    if k==0
        plot(real(k)+[0 1 1 0 0], imag(k)+[0 0 1 1 0], 'k-', 'linewidth', 1)
    else
        plot(real(k)+[0 1 1 0 0], imag(k)+[0 0 1 1 0], 'k--', 'linewidth', 1)
    end
end
text(1.65, 2.27, char(hex2dec('2113')))
plot([3+3j 3+2.45j], 'k:')
theta = 3*pi/2 - linspace(0,theta_example);
plot(3+3j + .35*exp(1j*theta), 'k-')
text(2.75, 2.51, char(hex2dec('03B8')))


%% #64. Figure for proof of probmax: funt odd, off-diagonal
clf
hold on
set(gca, 'layer', 'top')
axis equal
grid on
plot([-.4 5], [0 0], 'k')
plot([-.4j 5j], 'k')
set(gca, 'GridAlpha', .4)
xticks(-1:5), yticks(-1:5)
xticklabels({'' '0' '1' ['(t' char(hex2dec('2212')) '3)/2'] '' '' ''})
yticklabels({'' '0' '1' '' '' ['(t' '+' '1)/2'] ''})
xlabel x 
ylabel y
axis([-1 5 -1 5])
L = 3.5;
theta_0 = asin(1/L);
theta = 3*pi/2 - linspace(10*pi/180,theta_0,50);
plot(2+4j + L*exp(1j*theta), '--', 'linewidth', .75)
theta_1 = acos(3/L);
theta = 3*pi/2 - linspace(theta_1,39*pi/180,50);
set(gca, 'colororderindex', 1)
plot(2+4j + L*exp(1j*theta), '--', 'linewidth', .75)
theta = 3*pi/2 - linspace(theta_0,theta_1,50);
set(gca, 'colororderindex', 1)
plot(2+4j + L*exp(1j*theta), '-', 'linewidth', .75)
plot([2+4j], 'k.', 'markersize', 7)
theta_example = theta_0 + (theta_1-theta_0)*.55;
endpoint = 2+4j + L*exp(1j*(3*pi/2 - theta_example));
plot(endpoint, 'k.', 'markersize', 7)
plot([endpoint real(endpoint)+1j], 'k:')
plot([endpoint 1+1j*imag(endpoint)], 'k:')
patch([real(endpoint) 1 1 real(endpoint) real(endpoint)], [imag(endpoint) imag(endpoint) 1 1 imag(endpoint)], .88*[1 1 1], 'edgecolor', 'none')
plot([2+4j endpoint], 'k')
for k = [0 2+4j]
    if k==0
        plot(real(k)+[0 1 1 0 0], imag(k)+[0 0 1 1 0], 'k-', 'linewidth', 1)
    else
        plot(real(k)+[0 1 1 0 0], imag(k)+[0 0 1 1 0], 'k--', 'linewidth', 1)
    end
end
text(1.09, 2.77, char(hex2dec('2113')))
plot([2+4j 2+3.45j], 'k:')
theta = 3*pi/2 - linspace(0,theta_example);
plot(2+4j + .35*exp(1j*theta), 'k-')
text(1.81, 3.45, char(hex2dec('03B8')))


%% #65. Like #60 but with the g function
% Checked

clear, clc
a = 1;
L = .1:.1:10; %3.8:.1:4.2;%3:.1:3.5;
N = 3e5;
ii = ceil(L/sqrt(2)) + 1;
jj = ceil(sqrt(L.^2 - (ii-2).^2 )) + 1;
funt = ii+jj-1 % formula valid for real-valued lengths
z0 = rand(N, numel(L)) + 1j*rand(N, numel(L));
z1 = z0 + L.*exp(1j*2*pi*rand(N, numel(L)));
ii = ceil(abs(real(z1))) + (real(z1)<0); % add 1 for negative
jj = ceil(abs(imag(z1))) + (imag(z1)<0); % add 1 for negative
result = mean(ii+jj-1 == funt , 1);
plot(L, result, 'o', 'markersize', 4), hold on, grid on

g = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;
L_fine = linspace(.01,10,1e4);
ii = ceil(L_fine/sqrt(2)) + 1;
jj = ceil(sqrt(L_fine.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
ind_odd = (mod(funt_fine,2)==1) & (L_fine > (funt_fine-3)/sqrt(2)) & (L_fine <= sqrt((funt_fine-5).^2+(funt_fine-1).^2)/2);
result_computed(ind_odd) = g(L_fine(ind_odd)/a, funt_fine(ind_odd)-3, funt_fine(ind_odd)-3);
ind_odd2 = (mod(funt_fine,2)==1) & (L_fine > sqrt((funt_fine-5).^2+(funt_fine-1).^2)/2) & (L_fine <= sqrt((funt_fine-3).^2+(funt_fine-1).^2)/2);
result_computed(ind_odd2) = g(L_fine(ind_odd2)/a, funt_fine(ind_odd2)-3, funt_fine(ind_odd2)-3) + 2*g(L_fine(ind_odd2)/a, funt_fine(ind_odd2)-5, funt_fine(ind_odd2)-1);
ind_even = (mod(funt_fine,2)==0) & (L_fine > sqrt((funt_fine-4).^2+(funt_fine-2).^2)/2) & (L_fine <= sqrt((funt_fine-6).^2+funt_fine.^2)/2);
result_computed(ind_even) = 2*g(L_fine(ind_even)/a, funt_fine(ind_even)-4, funt_fine(ind_even)-2);
ind_even2 = (mod(funt_fine,2)==0) & (L_fine > sqrt((funt_fine-6).^2+funt_fine.^2)/2) & (L_fine <= (funt_fine-2)/sqrt(2));
result_computed(ind_even2) = 2*g(L_fine(ind_even2)/a, funt_fine(ind_even2)-4, funt_fine(ind_even2)-2) + 2*g(L_fine(ind_even2)/a, funt_fine(ind_even2)-6, funt_fine(ind_even2));
plot(L_fine.*(ind_odd|ind_odd2)./(ind_odd|ind_odd2), result_computed.*(ind_odd|ind_odd2)./(ind_odd|ind_odd2), '-', 'linewidth', .7), hold on, grid on
plot(L_fine.*(ind_even|ind_even2)./(ind_even|ind_even2), result_computed.*(ind_even|ind_even2)./(ind_even|ind_even2), '-', 'linewidth', .7), hold on, grid on
legend({'Simulated' 'Computed, odd max number of squares' 'Computed, even max number of squares'})


%% # 66. Tests, symbolic. Not useful: Matlab cannot compute the limit. I move to Maple, which can

syms t
L = sqrt((t-3).^2+(t-1).^2)/2;
G = (acos((t-3)/2./sqrt((t-3).^2+(t-1).^2)/2)-asin((t-3)/2./sqrt((t-3).^2+(t-1).^2)/2)).*(t-3).*(t-3) + ((t-3).^2+(t-1).^2)/2 + (t-3).^2/2 + (t-3).^2/2 - (t-3).*sqrt((t-3).^2+(t-1).^2-(t-3).^2) - (t-3).*sqrt((t-3).^2+(t-1).^2-(t-3).^2);
G = (acos((t-3)/2/sqrt((t-3)^2+(t-1)^2)/2)-asin((t-3)/2/sqrt((t-3)^2+(t-1)^2)/2))*(t-3)*(t-3) + ((t-3)^2+(t-1)^2)/2 + (t-3)^2/2 + (t-3)^2/2 - (t-3)*sqrt((t-3)^2+(t-1)^2-(t-3)^2) - (t-3)*sqrt((t-3)^2+(t-1)^2-(t-3)^2);
limit(G*t, t, inf)
syms t, assume(t, 'real'), pretty(simplify(diff((acos((t-3)/hypot(t-3, t-1)) - asin((t-3)/hypot(t-1, t-3))), t)))


%% #67. Figure

clear, clc
a = 1;
g = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;
L_fine = .01:.01:16;
ii = ceil(L_fine/sqrt(2)) + 1;
jj = ceil(sqrt(L_fine.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
ind_odd = (mod(funt_fine,2)==1) & (L_fine > (funt_fine-3)/sqrt(2)) & (L_fine <= sqrt((funt_fine-5).^2+(funt_fine-1).^2)/2);
result_computed(ind_odd) = g(L_fine(ind_odd)/a, funt_fine(ind_odd)-3, funt_fine(ind_odd)-3);
ind_odd2 = (mod(funt_fine,2)==1) & (L_fine > sqrt((funt_fine-5).^2+(funt_fine-1).^2)/2) & (L_fine <= sqrt((funt_fine-3).^2+(funt_fine-1).^2)/2);
result_computed(ind_odd2) = g(L_fine(ind_odd2)/a, funt_fine(ind_odd2)-3, funt_fine(ind_odd2)-3) + 2*g(L_fine(ind_odd2)/a, funt_fine(ind_odd2)-5, funt_fine(ind_odd2)-1);
ind_even = (mod(funt_fine,2)==0) & (L_fine > sqrt((funt_fine-4).^2+(funt_fine-2).^2)/2) & (L_fine <= sqrt((funt_fine-6).^2+funt_fine.^2)/2);
result_computed(ind_even) = 2*g(L_fine(ind_even)/a, funt_fine(ind_even)-4, funt_fine(ind_even)-2);
ind_even2 = (mod(funt_fine,2)==0) & (L_fine > sqrt((funt_fine-6).^2+funt_fine.^2)/2) & (L_fine <= (funt_fine-2)/sqrt(2));
result_computed(ind_even2) = 2*g(L_fine(ind_even2)/a, funt_fine(ind_even2)-4, funt_fine(ind_even2)-2) + 2*g(L_fine(ind_even2)/a, funt_fine(ind_even2)-6, funt_fine(ind_even2));
hold on, grid on, box on, set(gca, 'colororderindex', 2)
plot(L_fine.*(ind_odd|ind_odd2)./(ind_odd|ind_odd2), result_computed.*(ind_odd|ind_odd2)./(ind_odd|ind_odd2), '-')
set(gca, 'colororderindex', 1)
plot(L_fine.*(ind_even|ind_even2)./(ind_even|ind_even2), result_computed.*(ind_even|ind_even2)./(ind_even|ind_even2), '-')
%plot(L_fine, result_computed, '-')
legend({'Odd maximum number of squares' 'Even maximum number of squares'})
xlabel(char(hex2dec('2113'))), ylabel([char(hex2dec('03C1')) '(' char(hex2dec('2113')) ')'])


%% #68. Like #61 but including "a" (grid spacing)

a = 2;
t = 4:2:7002; L_fine = a*(t-2)/sqrt(2)-1e-10;
ii = ceil(L_fine/sqrt(2)/a) + 1;
jj = ceil(sqrt(L_fine.^2/a^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
ind_even = mod(funt_fine,2)==0; % even maximum number of squares
ind_even2 = ind_even & (2*L_fine.^2/a^2>funt_fine.^2-6*funt_fine+18);
alpha_even = acos(a*(funt_fine(ind_even)-2)/2./L_fine(ind_even));
beta_even = asin(a*(funt_fine(ind_even)-4)/2./L_fine(ind_even));
alphap_even2 = acos(a*funt_fine(ind_even2)/2./L_fine(ind_even2));
betap_even2 = asin(a*(funt_fine(ind_even2)-6)/2./L_fine(ind_even2));
result_computed(ind_even) = ( (alpha_even-beta_even).*(funt_fine(ind_even)-2).*(funt_fine(ind_even)-4) + (funt_fine(ind_even)-2).^2/2 + (funt_fine(ind_even)-4).^2/2 + ...
    2*L_fine(ind_even).^2/a^2 - (funt_fine(ind_even)-2).*sqrt(4*L_fine(ind_even).^2/a^2 - (funt_fine(ind_even)-4).^2) - (funt_fine(ind_even)-4).*sqrt(4*L_fine(ind_even).^2/a^2 - (funt_fine(ind_even)-2).^2) ) / pi;
result_computed(ind_even2) = result_computed(ind_even2) + ...
    ( (alphap_even2-betap_even2).*funt_fine(ind_even2).*(funt_fine(ind_even2)-6) + funt_fine(ind_even2).^2/2 + (funt_fine(ind_even2)-6).^2/2 + ...
    2*L_fine(ind_even2).^2/a^2 - funt_fine(ind_even2).*sqrt(4*L_fine(ind_even2).^2/a^2 - (funt_fine(ind_even2)-6).^2) - (funt_fine(ind_even2)-6).*sqrt(4*L_fine(ind_even2).^2/a^2 - funt_fine(ind_even2).^2) ) / pi;
plot(L_fine, result_computed.*L_fine), grid on


%% #69. Like the first part of #46 (that is, funt) but including formula found by Alex for funt in the square case with real-valued lengths (e-mail 26viii21)

clear
clc
L = .001:.001:300;
G = 1:35;
a = 1.5754;
b = a;

ii = ceil((a + b*real(sqrt(4*L.^2/(a^2+b^2)-1))) / 2/a) + 1;
jj = ceil(sqrt( L.^2 - (ii-2).^2*a^2 )/b) + 1;
funt_rect = ii+jj-1;

ii = ceil((1 + real(sqrt(2*L.^2/a^2-1)))/2) + 1;
jj = ceil(sqrt(L.^2/a^2 - (ii-2).^2 )) + 1;
funt_sq = ii+jj-1;

ii = ceil(L/a/sqrt(2)) + 1;
jj = ceil(sqrt(L.^2/a^2 - (ii-2).^2 )) + 1;
funt_sq_2 = ii+jj-1;

funt_sq_3 = floor(sqrt(2*ceil(L.^2/a^2) - 2)) + 3; % formula found by Alex
isequal(funt_rect, funt_sq, funt_sq_2, funt_sq_3)


%% #70. Figure for proof of {the proposition that gives the formula of the minimum sufficient set} based on a continuous solution that is then rounded (Alex)
% Based on #43 (with flag_hexagram = false; flag_fun = false; flag_pairs = true)

clear, clf, clc
flag_pairs = true;
a = 1.35; b = 1;
hold on
set(gca, 'layer', 'top')
grid on
plot([2*a+2j*b 2*a+5.7j*b], 'k') % "secondary" axes
plot([2*a+2j*b 5.7*a+2j*b], 'k')
chosen = [2*a+2j*b 2*a+3j*b 3*a+3j*b 3*a+4j*b 3*a+5j*b 4*a+5j*b];
chosen_2 = 4*a+4j*b;
ms = 5;
plot(chosen, 'ko', 'markerfacecolor', 'k', 'markersize', ms) % chosen (i,j) pairs
plot(chosen_2, 'ko', 'markerfacecolor', 'w', 'markersize', ms) % chosen 2
plot(setdiff((2:5)*a + 1j*(2:5).'*b, [chosen chosen_2]), 'ko', 'markerfacecolor', 'none', 'markersize', ms) % non-chosen (i,j) pairs
axis equal %pbaspect([1 1 1])
axis([1*a 6*a 1*b 6*b])
plot([1.8 2.4]*a, [2.2 1.6]*b, 'k--')
plot([1.8 3.4]*a, [3.2 1.6]*b, 'k--')
plot([1.8 4.4]*a, [4.2 1.6]*b, 'k--')
plot([1.8 5.4]*a, [5.2 1.6]*b, 'k--')
plot([2.5 5.4]*a, [5.5 2.6]*b, 'k--') % limit to secondary axes
plot([3.5 5.4]*a, [5.5 3.6]*b, 'k--')
plot([4.5 5.4]*a, [5.5 4.6]*b, 'k--')
xl = xlim; yl = ylim;
set(gca, 'GridAlpha', .4, 'XTick', xl(1):xl(2), 'YTick', yl(1):yl(2))
xticks((0:7)*a), yticks((0:9)*b)
xticklabels([{'0'} arrayfun(@(n)[num2str(n) char(8201) 'a'], 1:7, 'UniformOutput', false)]) % thin space
yticklabels([{'0'} arrayfun(@(n)[num2str(n) char(8201) 'b'], 1:9, 'UniformOutput', false)])
xlabel(['i' char(8201) 'a']), ylabel(['j' char(8201) 'b']) % thin space
set(gca, 'colororderindex', 1)
xl = [1.85 3.95];
plot(xl*a, (xl-2)*a^2/b^2+2, '-') % line
text(6.96, 2.43, ['i+j' char(8211) '1 = t'])
set(gca, 'colororderindex', 1)
t = 7; xl = (t-3)*b^2/(a^2+b^2) + 2; yl = (t-3)*a^2/(a^2+b^2) + 2; % computed
ms = 6;
plot(xl*a, yl*b, 's', 'markersize', ms)
text(4.78, 4.69, '(i^+,j^+)')
text(3.52, 4.82, '(i_t,j_t)')

