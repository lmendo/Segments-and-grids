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
flag_pairs = 2; % 0: no; 1: yes, normalized labels, 2: yes, not normalized
a = 1.35; b = 1;
hold on
set(gca, 'layer', 'top')
grid on
plot([2*a+2j*b 2*a+5.7j*b], 'k') % "secondary" axes
plot([2*a+2j*b 5.7*a+2j*b], 'k')
chosen = [2*a+2j*b 2*a+3j*b 3*a+3j*b 3*a+4j*b 3*a+5j*b 4*a+5j*b];
chosen_2 = 4*a+4j*b;
if flag_hexagram ms = 4.5; else ms = 5; end
plot([1.8 2.4]*a, [2.2 1.6]*b, 'k--')
plot([1.8 3.4]*a, [3.2 1.6]*b, 'k--')
plot([1.8 4.4]*a, [4.2 1.6]*b, 'k--')
plot([1.8 5.4]*a, [5.2 1.6]*b, 'k--')
plot([2.5 5.4]*a, [5.5 2.6]*b, 'k--') % limit to secondary axes
plot([3.5 5.4]*a, [5.5 3.6]*b, 'k--')
plot([4.5 5.4]*a, [5.5 4.6]*b, 'k--')
if flag_pairs==1
    if flag_in text(4.5, 3.89, '(i^+,j^+)'), text(4.16-~flag_hexagram*.05, 4.2, '(i^\ast,j^\ast)')
    else text(5.61, 4.42 , '(i^+,j^+)'), text(5.13-~flag_hexagram*.02, 3.76+~flag_hexagram*.03, '(i^\ast,j^\ast)')
    end
elseif flag_pairs==2
    fs = 10; % fontsize
    if flag_in text(4.24, 3.9, '(i^+a, j^+b)', 'fontsize', fs, 'Backgroundcolor', 'w'), text(3.28-~flag_hexagram*.05, 4.18, '(i^\ast{}a, j^\ast{}b)', 'fontsize', fs, 'Backgroundcolor', 'w')
    else
        text(5.65, 4.45 , '(i^+a, j^+b)', 'fontsize', fs, 'Backgroundcolor', 'w'), text(5.15-~flag_hexagram*.02, 3.76+~flag_hexagram*.03, '(i^\ast{}a, j^\ast{}b)', 'fontsize', fs,'Backgroundcolor', 'w')
        %text(5.61, 4.42 , ['(i' char(8314) 'a, j' char(8314) 'a)'], 'fontsize', fs, 'Backgroundcolor', 'w'), text(5.13-~flag_hexagram*.02, 3.76+~flag_hexagram*.03, '(i^\ast{}a, j^\ast{}a)', 'fontsize', fs,'Backgroundcolor', 'w')
    end
end
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
set(gcf, 'Position', [680 558 560 400])


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


%% #67. Figure: probmax

clear, clc
a = 1;
%g= @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi; % #65
g = @(r,u,v) ((acos(v/2./r)-asin(u/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi; % #74
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
legend({'Odd maximum number of tiles' 'Even maximum number of tiles'})
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
text(6.96, 2.43, ['i+j' char(8211) '1 = t'])
t = 7; xl = (t-3)*b^2/(a^2+b^2) + 2; yl = (t-3)*a^2/(a^2+b^2) + 2; % computed
%text(4.78, 4.69, '(i^+,j^+)')
text(4.71, 4.5, '(i^+a, j^+b)', 'BackgroundColor', 'w')
%text(3.52, 4.82, '(i_t,j_t)')
text(3.27, 4.83, '(i_ta, j_tb)', 'BackgroundColor', 'w')
ms = 6;
plot(xl*a, yl*b, 's', 'markersize', ms)
set(gca, 'colororderindex', 1)
xl = [1.85 3.95];
plot(xl*a, (xl-2)*a^2/b^2+2, '-') % line
chosen = [2*a+2j*b 2*a+3j*b 3*a+3j*b 3*a+4j*b 3*a+5j*b 4*a+5j*b];
chosen_2 = 4*a+4j*b;
ms = 5;
plot(chosen, 'ko', 'markerfacecolor', 'k', 'markersize', ms) % chosen (i,j) pairs
plot(chosen_2, 'ko', 'markerfacecolor', 'w', 'markersize', ms) % chosen 2
plot(setdiff((2:5)*a + 1j*(2:5).'*b, [chosen chosen_2]), 'ko', 'markerfacecolor', 'none', 'markersize', ms) % non-chosen (i,j) pairs


%% #71. Check for new formula for funt:

clear, clc
t = 18;
a = 1.56;
b = 1.15;
iplus = max((t-3)*b^2/(a^2+b^2), 0)+2;
jplus = max((t-3)*a^2/(a^2+b^2), 0)+2;
istar = floor(iplus+1/2);
jstar = ceil(jplus-1/2);
lambda = sqrt((istar-2)^2*a^2 + (jstar-2)^2*b^2) % original formula
%lambda_2 = sqrt( (t-3)^2*a^2*b^2/(a^2+b^2) + (istar-iplus)^2*a^2 + (jstar-jplus)^2*b^2 ) % new formula, corrected
lambda_2 = sqrt( (t-3)^2*a^2*b^2/(a^2+b^2) + (istar-iplus)^2*(a^2+b^2) ) % new formula, corrected


%% #72. Ratio of asymptotic slopes (with respect to length) for maximum and average numbers of visited tiles, as a function of a/b

clear, clc, close all
ratio_ab = logspace(-2,2,300);
as_ratio_slopes = 2/pi*(1+ratio_ab)./hypot(1, ratio_ab);
plot(ratio_ab, as_ratio_slopes)
set(gca, 'xscale', 'log'), grid on, xticklabels(xticks), xlabel('a/b'), ylabel([char(hex2dec('03C3')) '(a/b)']), set(gcf, 'Position', [360 300 560 310])


%% #73. Maximum and average numbers of visited tiles as a function of length, with (a,b) as a parameter.
% Part taken from #42

clear, clc, close all
L = .01:.01:200;
a = 2.35; b = 1;
ii = ceil((a + b*real(sqrt(4*L.^2/(a^2+b^2)-1))) / 2/a) + 1;
jj = ceil(sqrt( L.^2 - (ii-2).^2*a^2 )/b) + 1;
funt = ii+jj-1;
funta = 2*L/pi*(1/a+1/b)+1;
as_ratio_slopes = 2/pi*(1+a/b)./hypot(1, a/b)
%plot(L, funt), hold on, plot(L, funta);
plot(funta./funt), hold on, plot(xlim, repmat(as_ratio_slopes, 1, 2))


%% #74. Like #65 but defining the g function differently (gives the same result)
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

%g= @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi; % #65
g = @(r,u,v) ((acos(v/2./r)-asin(u/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi; % #74
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


%% #75. 27iii23: a reviewer finds a mistake in the probability of visiting the maximum number of tiles

% From #65

% Simulation
clear, clc
a = 1;
L = 6.35; %3.8:.1:4.2;%3:.1:3.5;
N = 2e7;
ii = ceil(L/sqrt(2)) + 1;
jj = ceil(sqrt(L.^2 - (ii-2).^2 )) + 1;
funt = ii+jj-1 % formula valid for real-valued lengths
z0 = rand(N, numel(L)) + 1j*rand(N, numel(L));
z1 = z0 + L.*exp(1j*2*pi*rand(N, numel(L)));
ii = ceil(abs(real(z1))) + (real(z1)<0); % add 1 for negative
jj = ceil(abs(imag(z1))) + (imag(z1)<0); % add 1 for negative
result = mean(ii+jj-1 == funt , 1)

% Theoretical
g = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;
L_fine = L;
ii_theo = ceil(L_fine/sqrt(2)) + 1;
jj_theo = ceil(sqrt(L_fine.^2 - (ii_theo-2).^2 )) + 1;
funt_fine = ii_theo+jj_theo-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
ind_odd = (mod(funt_fine,2)==1) & (L_fine > (funt_fine-3)/sqrt(2)) & (L_fine <= sqrt((funt_fine-5).^2+(funt_fine-1).^2)/2);
result_computed(ind_odd) = g(L_fine(ind_odd)/a, funt_fine(ind_odd)-3, funt_fine(ind_odd)-3);
ind_odd2 = (mod(funt_fine,2)==1) & (L_fine > sqrt((funt_fine-5).^2+(funt_fine-1).^2)/2) & (L_fine <= sqrt((funt_fine-3).^2+(funt_fine-1).^2)/2);
result_computed(ind_odd2) = g(L_fine(ind_odd2)/a, funt_fine(ind_odd2)-3, funt_fine(ind_odd2)-3) + 2*g(L_fine(ind_odd2)/a, funt_fine(ind_odd2)-5, funt_fine(ind_odd2)-1);
ind_even = (mod(funt_fine,2)==0) & (L_fine > sqrt((funt_fine-4).^2+(funt_fine-2).^2)/2) & (L_fine <= sqrt((funt_fine-6).^2+funt_fine.^2)/2);
result_computed(ind_even) = 2*g(L_fine(ind_even)/a, funt_fine(ind_even)-4, funt_fine(ind_even)-2);
ind_even2 = (mod(funt_fine,2)==0) & (L_fine > sqrt((funt_fine-6).^2+funt_fine.^2)/2) & (L_fine <= (funt_fine-2)/sqrt(2));
result_computed(ind_even2) = 2*g(L_fine(ind_even2)/a, funt_fine(ind_even2)-4, funt_fine(ind_even2)-2) + 2*g(L_fine(ind_even2)/a, funt_fine(ind_even2)-6, funt_fine(ind_even2));
result_computed

% Histogram of ii, jj from simulation
figure
ind = ii+jj-1==11 & real(z1)>0 & imag(z1)>0; % ii, jj that produce the maximum number of visited tiles (11).
% And only positive coordinates for simplicity
mean(ind)
histogram2(ii(ind), jj(ind), 'binmethod', 'integer', 'DisplayStyle' ,'tile', 'ShowEmptyBins', 'off')

% : there are values with abs(ii-jj) greater than 2, for example ii=4, jj=8, and ii=8, ii=4. Their probability is so
% small (but non-zero) that the simulated values seem to match with the theorerical computations.
% That's why I didn't catch this mistake with the simulations.


%% # 76. Like #65 but correcting the mistake by including more terms

% The computations here seem to be correct. For L around 10 the difference between the old (wrong) and new
% computation is already noticeable, and the simulation results match the new, not the old.
%  Computing for L_fine large enough (around 1000, with fine enough sampling) and using log-log
% scale the 1/sqrt(L) asymptotic variation of the probability maxima is clearly visible.

clear, clc
L = .1:.1:10; %3.8:.1:4.2;%3:.1:3.5;
N = 3e5;
ii = ceil(L/sqrt(2)) + 1;
jj = ceil(sqrt(L.^2 - (ii-2).^2 )) + 1;
funt = ii+jj-1; % formula valid for real-valued lengths
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
num_k =  NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_fine(ind_fine);
    if mod(t,2) % odd
        res = g(L, t-3, t-3); % this term is always present
        potential_k = 1:(t-3)/2;
        %%potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-3)/2-potential_k).^2 + ((t-3)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-3-2*k, t-3+2*k);
        end
    else % even
        res = 2*g(L, t-4, t-2); % this term is always present
        potential_k = 1:(t-4)/2;
        %%potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-4)/2-potential_k).^2 + ((t-2)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-4-2*k, t-2+2*k);
        end
    end
    result_computed(ind_fine) = res;
    num_k(ind_fine) = numel(actual_k);
end
ind_odd = mod(funt_fine,2)==1; % for plotting
ind_even = ~ind_odd;
plot(L_fine.*ind_odd./ind_odd, result_computed.*ind_odd./ind_odd, '-', 'linewidth', .7), hold on, grid on
plot(L_fine.*ind_even./ind_even, result_computed.*ind_even./ind_even, '-', 'linewidth', .7), hold on, grid on
xlabel length, ylabel probability
legend({'Simulated' 'Computed, odd max number of squares' 'Computed, even max number of squares'})


%% # 77. Like #76 but looking into the asymptotic behaviour

g = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;
L_fine = linspace(10,1000,1e5);
ii = ceil(L_fine/sqrt(2)) + 1;
jj = ceil(sqrt(L_fine.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
num_k =  NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_fine(ind_fine);
    if mod(t,2) % odd
        res = g(L, t-3, t-3); % this term is always present
        potential_k = 1:(t-3)/2;
        %%potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-3)/2-potential_k).^2 + ((t-3)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-3-2*k, t-3+2*k);
        end
    else % even
        res = 2*g(L, t-4, t-2); % this term is always present
        potential_k = 1:(t-4)/2;
        %%potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-4)/2-potential_k).^2 + ((t-2)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-4-2*k, t-2+2*k);
        end
    end
    result_computed(ind_fine) = res;
    num_k(ind_fine) = numel(actual_k);
end
ind_odd = mod(funt_fine,2)==1; % for plotting
ind_even = ~ind_odd;
plot(L_fine.*ind_odd./ind_odd, result_computed.*ind_odd./ind_odd, '-', 'linewidth', .7), hold on, grid on
plot(L_fine.*ind_even./ind_even, result_computed.*ind_even./ind_even, '-', 'linewidth', .7), hold on, grid on
xlabel length, ylabel probability

ind = 984891; % manually choose: near a maximum, with even t for example
t = funt_fine(ind);
L = L_fine(ind);
potential_k = 1:(t-4)/2;
ind_actual_k = ((t-4)/2-potential_k).^2 + ((t-2)/2+potential_k).^2 < L^2;
actual_k = potential_k(ind_actual_k);
numel(actual_k), sqrt(t/2), sqrt(L)/2^.25 % checked: similar
k = 0; g(L, t-4-k, t-2+k)*t, 2/3/pi % checked: similar
k = 4; g(L, t-4-k, t-2+k)*t, 2/3/pi % checked: similar

result_computed(ind), sqrt(2)/3/pi/sqrt(t), 2^.25/3/pi/sqrt(L) % checked: similar, and similar to maxima seen in curves; a little less, around 91.5%

plot(L_fine, 2^.25/3/pi./sqrt(L_fine)) % this line is a little above the maxima. Perhaps it gets a little closer as L increses
plot(L_fine, 2^.5/3/pi./sqrt(funt_fine)) % similar to the above line
% I think most of the difference comes from g() being different from 2/3/pi; the difference there is noticeable too


%% # 78. Similar to the above, but focusing on L near the maxima
% The ratio tends to about .915

clear, close all
g = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;
t = 11:9999;
epsilon = 1e-6;
ind = mod(t,2)==1;
L_fine(ind) = hypot(t(ind)-3,t(ind)-1)/2 - epsilon;
ind = mod(t,2)==0;
L_fine(ind) = (t(ind)-2)/sqrt(2)- epsilon;
ii = ceil(L_fine/sqrt(2)) + 1;
jj = ceil(sqrt(L_fine.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
num_k =  NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_fine(ind_fine);
    if mod(t,2) % odd
        res = g(L, t-3, t-3); % this term is always present
        potential_k = 1:(t-3)/2;
        %%potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-3)/2-potential_k).^2 + ((t-3)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-3-2*k, t-3+2*k);
        end
    else % even
        res = 2*g(L, t-4, t-2); % this term is always present
        potential_k = 1:(t-4)/2;
        %%potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-4)/2-potential_k).^2 + ((t-2)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-4-2*k, t-2+2*k);
        end
    end
    result_computed(ind_fine) = res;
    num_k(ind_fine) = numel(actual_k);
end

ind_odd = mod(funt_fine,2)==1; % for plotting
ind_even = ~ind_odd;
plot(L_fine.*ind_odd./ind_odd, result_computed.*ind_odd./ind_odd, '.', 'linewidth', .7), hold on, grid on
plot(L_fine.*ind_even./ind_even, result_computed.*ind_even./ind_even, '.', 'linewidth', .7), hold on, grid on
xlabel length, ylabel probability
plot(L_fine, 2^.25/3/pi./sqrt(L_fine))
set(gca, 'XScale', 'log', 'YScale', 'log') % a gap is visible

%figure
%plot(L_fine, result_computed./( 2^.25/3/pi./sqrt(L_fine))) % tends to about .915

%figure
%plot(L_fine, num_k./(sqrt(L_fine)/2^.25)) % checked: tends to 1


%% #79. Let's see the limit of function g
% I see that, even if t*g seems to tend to 1 for any k (up to numerical precision), it tends more
% slowly to that value as k is increased.
%   So maybe the 0.915 is caused by the fact that new values of k are included as t grows, and for
% those new values t*g is further from convergence. So for any t there are always summands t*g that
% are far from convergence
%   Yes. The sum of terms divided by the estimated number of terms sqrt(t/2) gives around 0.4571,
% which multiplied by 2 is 0.914. There was a factor of 2 or almost 2 missing, because the terms for
% k>0 should appear twice in the sum, as should that for k=0 too for t even. So the estimated number of
% terms should be sqrt(2*t), not sqrt(t/2).

clear, close all
g = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;
t_all = 1000:100000;
k = 0;

%plot(t_all, t_all.*g(hypot(t_all-3,t_all-1)/2, t_all-3-2*k, t_all-3+2*k) / (2/3/pi)); % seems to be 1, up to numerical precision
plot(t_all, t_all.*g((t_all-2)/sqrt(2), t_all-4-2*k, t_all-2+2*k) / (2/3/pi)); % seems to be 1, up to numerical precision

plot(t_all, g((t_all-2)/sqrt(2), t_all-4-2*k, t_all-2+2*k) ./ t_all.*g((t_all-2)/sqrt(2), t_all-4, t_all-2));

t = 10000;
if mod(t,2) % odd
    L = hypot(t-3,t-1)/2;
    potential_k = 1:(t-3)/2;
    ind_actual_k = ((t-3)/2-potential_k).^2 + ((t-3)/2+potential_k).^2 < L^2;
    actual_k = potential_k(ind_actual_k);
else % even
    L = (t-2)/sqrt(2);
    potential_k = 1:(t-4)/2;
    ind_actual_k = ((t-4)/2-potential_k).^2 + ((t-2)/2+potential_k).^2 < L^2;
    actual_k = potential_k(ind_actual_k);
end
actual_k = [0 actual_k];
g_actual_k = NaN(1,numel(actual_k));
c = 0;
for k = actual_k
    c = c+1;
    if mod(t,2)
        g_actual_k(c) =  t.*g(hypot(t-3,t-1)/2, t-3-2*k, t-3+2*k) / (2/3/pi);
    else
        g_actual_k(c) =  t.*g((t-2)/sqrt(2), t-4-2*k, t-2+2*k) / (2/3/pi);
    end
end
plot(actual_k, g_actual_k, 'o-')
sum(g_actual_k/sqrt(t/2))
% The result is 0.4571, which multiplied by 2 is 0.914


%% #80. Aproximations for g*L

clear, close all
g_exact = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;

t = 10000;
k = 3;

m = mod(t,1);
L = hypot(t-2-m,t-2+m)/2;
L*sqrt(2) * g_exact(L, t-4+m-2*k, t-2-m+2*k) / (2/3/pi)
g_aprox = @(r,u,v) (sqrt(r^2-v^2/4)-u/2)*(sqrt(r^2-u^2/4)-v/2) * hypot(sqrt(r^2-v^2/4)-u/2,sqrt(r^2-u^2/4)-v/2)/r /3/pi;
% It would be: g_aprox = @(r,u,v) x*y * hypot(x,y)/r /3/pi, with x = sqrt(r^2-v^2/4)-u/2, y = sqrt(r^2-u^2/4)-v/2
L*sqrt(2) * g_aprox(L, t-4+m-2*k, t-2-m+2*k) / (2/3/pi) % similar to the result using g_exact

% Drawing similar to that of the reviewer
n = 40; L = hypot(n,n+1) - 1e-6; %L = sqrt(41) - 1e-6;
theta = linspace(0,2*pi,1000);
fL = floor(L);
plot(1+1j+L*exp(2*j*pi*theta))
xl = xlim; yl = ylim; xticks(xl(1):xl(end)), yticks(yl(1):yl(end)), axis square, grid on
hold on
t = floor(sqrt(2*ceil(L^2)-2))+3;
plot([1 t-1-yl(2):xl(2)], [1 yl(2):-1:t-1-xl(2)], '.')


%% #81

clear, close all
g_exact = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;

t = 1001; % odd
k=3;
L = hypot(t-3,t-1)/2;
g_exact(L, t-3-2*k, t-3+2*k) / (2/3/pi)
v = t-3+2*k;
u = t-3-2*k;
x = sqrt(L^2-v^2/4)-u/2
x = sqrt((t-1)^2/4-k^2-k*(t-3)) - (t-3)/2 + k % checked
y = sqrt(L^2-u^2/4)-v/2
y = sqrt((t-1)^2/4-k^2+k*(t-3)) - (t-3)/2 - k % checked
x^2 + y^2
2*L^2 - u*sqrt(L^2-v^2/4) - v*sqrt(L^2-u^2/4) % checked
x = (t-1)/2 * sqrt(1-(4*k^2/(t-1)^2+4*k/(t-1)^2*(t-3))) - (t-3)/2 + k % checked
x = (t-1)/2 * sqrt(1-(4*k^2+4*k*(t-3))/(t-1)^2) - (t-3)/2 + k % checked
x = (t-1)/2 * (1-2*(k^2+k*(t-3))/(t-1)^2) - (t-3)/2 + k % sqrt(1-z) ~ 1-z/2: not good aprox here


%% #82

clear
t = 1e4; % the result for t = 1e5 should give around 1.37e-3, and for t = 1e6 around 0.44e-3, computed from #78
xmax = sqrt(t);
a = sqrt(2)/2/t;
c = sqrt(2)/2;

f = @(m) sqrt(1+m.^2)./(1+m).^2./(1-m).^2.*(-m.^2/4/a+c).^3;
I = 2*sqrt(2)/pi/t/6 * 4*sqrt(2)/2/a*integral(f, 0, 2*a*xmax) % checked! Result matches the
% probability rho(ell_t) in the paper
I = 4/3/pi/t/a * integral(f, 0, 2*a*xmax) % checked

f = @(m) (1+m.^2/2)./(1+m).^2./(1-m).^2.*(-m.^2/4/a+c).^3; % approximating sqrt(1+m.^2) as 1+m.^2/2
I = 4/3/pi/t/a * integral(f, 0, 2*a*xmax) % checked


%% #83 The same integral but symbolically

clear
syms t m pi % pi needes so it treats it a symbolic
xmax = sqrt(t);
a2 = 1/2/t^2;
c2 = 1/2;
f = sqrt(1+m.^2)./(1+m).^2./(1-m).^2.*(-m.^2/4/sqrt(a2)+sqrt(c2)).^3;
%f = sqrt(1+m.^2)./(1+m).^2./(1-m).^2;
I = 4/3/pi/t/sqrt(a2) * int(f, m); % doesn't give an expression

f = (1+m.^2)./(1+m).^2./(1-m).^2.*(-m.^2/4/sqrt(a2)+sqrt(c2)).^3; % approximating sqrt(1+m.^2) as 1+m.^2/2
I = 4/3/pi/t/sqrt(a2) * int(f, m); % does give an expression
double(eval(subs(subs(I,m,sqrt(2/t_val)),t,t_val))) % checked!
pretty(simplify(I))

I = int(f, m);
double(4*sqrt(2)/3/pi * eval(subs(subs(I,m,sqrt(2/t_val)),t,t_val))) % checked

f = (1+5/2*m.^2).*(-m.^2/4/sqrt(a2)+sqrt(c2)).^3; % approximating more terms
I = int(f, m);
double(4*sqrt(2)/3/pi * eval(subs(subs(I,m,sqrt(2/t_val)),t,t_val))) % checked

f = (-m.^2/4/sqrt(a2)+sqrt(c2)).^3; % approximating more terms by 1
I = int(f, m);
assert(double(subs(I,m,0))==0)
double(4*sqrt(2)/3/pi * eval(subs(subs(I,m,sqrt(2/t_val)),t,t_val))) % checked

f = (1-m.^2/2*t).^3; % approximating more terms
I = int(f, m);
assert(double(subs(I,m,0))==0)
double(2/3/pi * eval(subs(subs(I,m,sqrt(2/t_val)),t,t_val))) % checked

syms n
f = (1-n).^3/sqrt(n)/sqrt(t); % approximating more terms
I = int(f, n);
assert(double(subs(I,n,0))==0)
double(sqrt(2)/3/pi * eval(subs(subs(I,n,1),t,t_val))) % checked

syms p
f = (1-p^2).^3/sqrt(t); % approximating more terms
I = int(f, p);
assert(double(subs(I,p,0))==0)
double(2*sqrt(2)/3/pi * eval(subs(subs(I,p,1),t,t_val))) % checked

pi_val = double(pi)
32*sqrt(2)/105/pi_val/sqrt(t_val)


%% #84. Computing probability non-asymptotically: number of terms in the sum

clear, close all

% From #77:
g = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;
L_fine = linspace(1,1000,1e6);
ii = ceil(L_fine/sqrt(2)) + 1;
jj = ceil(sqrt(L_fine.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
num_k =  NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_fine(ind_fine);
    if mod(t,2) % odd
        res = g(L, t-3, t-3); % this term is always present
        potential_k = 1:(t-3)/2;
        %%potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-3)/2-potential_k).^2 + ((t-3)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-3-2*k, t-3+2*k);
        end
    else % even
        res = 2*g(L, t-4, t-2); % this term is always present
        potential_k = 1:(t-4)/2;
        %%potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-4)/2-potential_k).^2 + ((t-2)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-4-2*k, t-2+2*k);
        end
    end
    result_computed(ind_fine) = res;
    num_k(ind_fine) = numel(actual_k);
end
plot(num_k)

% New: formula for number of extra terms:
num_k2 = ceil(sqrt(L_fine.^2/2-(funt_fine-3).^2/4)-(1-mod(funt_fine,2))/2)-1;
hold on
plot(num_k2, '--')
all(num_k==num_k2) % checked


%% #85. Computing probability non-asymptotically using the number of terms from #84

clear, close all
L_fine = linspace(10,1000,1e5);

% From #77
g = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;
ii = ceil(L_fine/sqrt(2)) + 1;
jj = ceil(sqrt(L_fine.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
num_k =  NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_fine(ind_fine);
    if mod(t,2) % odd
        res = g(L, t-3, t-3); % this term is always present
        potential_k = 1:(t-3)/2;
        %%potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-3)/2-potential_k).^2 + ((t-3)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-3-2*k, t-3+2*k);
        end
    else % even
        res = 2*g(L, t-4, t-2); % this term is always present
        potential_k = 1:(t-4)/2;
        %%potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-4)/2-potential_k).^2 + ((t-2)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-4-2*k, t-2+2*k);
        end
    end
    result_computed(ind_fine) = res;
    num_k(ind_fine) = numel(actual_k);
end
ind_odd = mod(funt_fine,2)==1; % for plotting
ind_even = ~ind_odd;
plot(L_fine.*ind_odd./ind_odd, result_computed.*ind_odd./ind_odd, '-', 'linewidth', .7), hold on, grid on
plot(L_fine.*ind_even./ind_even, result_computed.*ind_even./ind_even, '-', 'linewidth', .7), hold on, grid on
xlabel length, ylabel probability

% New:
g = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;
ind_odd = mod(funt_fine,2)==1;
ind_even = ~ind_odd;
num_k2 = ceil( sqrt(L_fine.^2/2-(funt_fine-3).^2/4) - (1+ind_even/2) );
result_computed2 = NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_fine(ind_fine);
    if ind_odd(ind_fine)
        result_computed2(ind_fine) = g(L, t-3, t-3); % this term is always present
        for k = 1:num_k2(ind_fine) % extra terms
            result_computed2(ind_fine) = result_computed2(ind_fine) + 2*g(L, t-3-2*k, t-3+2*k);
        end
    else
        result_computed2(ind_fine) = 2*g(L, t-4, t-2); % this term is always present
        for k = 1:num_k2(ind_fine) % extra terms
            result_computed2(ind_fine) = result_computed2(ind_fine) + 2*g(L, t-4-2*k, t-2+2*k);
        end
    end
end
plot(L_fine.*ind_odd./ind_odd, result_computed2.*ind_odd./ind_odd, 'k--', 'linewidth', .7)
plot(L_fine.*ind_even./ind_even, result_computed2.*ind_even./ind_even, 'k--', 'linewidth', .7)
% checked

% New, with number of terms not explicit
g = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;
ind_odd = mod(funt_fine,2)==1;
ind_even = ~ind_odd;
result_computed3 = NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_fine(ind_fine);
    potential_i = 0:t-3;
    potential_j = t-1-potential_i;
    ind_actual_ij = (potential_i-1).^2 + (potential_j-1).^2 < L^2;
    actual_i = potential_i(ind_actual_ij);
    actual_j = potential_j(ind_actual_ij);
    res = 0;
    for ind_ij = 1:numel(actual_i)
        ii = actual_i(ind_ij);
        jj = actual_j(ind_ij);
        res = res + g(L, 2*(ii-1), 2*(jj-1));
    end
    result_computed3(ind_fine) = res;
end
% checked

% Simulation, from #76
% I see the above computations match the simulation results
clear, clc
L = .1:.05:15; %3.8:.1:4.2;%3:.1:3.5;
N = 3e5;
ii = ceil(L/sqrt(2)) + 1;
jj = ceil(sqrt(L.^2 - (ii-2).^2 )) + 1;
funt = ii+jj-1; % formula valid for real-valued lengths
z0 = rand(N, numel(L)) + 1j*rand(N, numel(L));
z1 = z0 + L.*exp(1j*2*pi*rand(N, numel(L)));
ii = ceil(abs(real(z1))) + (real(z1)<0); % add 1 for negative
jj = ceil(abs(imag(z1))) + (imag(z1)<0); % add 1 for negative
result = mean(ii+jj-1 == funt , 1);
plot(L, result, 'o', 'markersize', 4), hold on, grid on

% New, like _computed2 but redefining g: I redefined u, v so that input (t-3) now becomes (t-3)/2, etc
g = @(r,u,v) (2*(acos(u./r)-asin(v./r)).*u.*v + r.^2 + u.^2 + v.^2 - 2*u.*sqrt(r.^2-v.^2) - 2*v.*sqrt(r.^2-u.^2))/pi;
ind_odd = mod(funt_fine,2)==1;
ind_even = ~ind_odd;
num_k2 = ceil( sqrt(L_fine.^2/2-(funt_fine-3).^2/4) - (1+ind_even/2) );
result_computed4 = NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_fine(ind_fine);
    if ind_odd(ind_fine)
        result_computed4(ind_fine) = g(L, (t-3)/2, (t-3)/2); % this term is always present
        for k = 1:num_k2(ind_fine) % extra terms
            result_computed4(ind_fine) = result_computed4(ind_fine) + 2*g(L, (t-3)/2-k, (t-3)/2+k);
        end
    else
        result_computed4(ind_fine) = 2*g(L, (t-4)/2, (t-2)/2); % this term is always present
        for k = 1:num_k2(ind_fine) % extra terms
            result_computed4(ind_fine) = result_computed4(ind_fine) + 2*g(L, (t-4)/2-k, (t-2)/2+k);
        end
    end
end
plot(L_fine.*ind_odd./ind_odd, result_computed4.*ind_odd./ind_odd, 'k--', 'linewidth', .7), hold on
plot(L_fine.*ind_even./ind_even, result_computed4.*ind_even./ind_even, 'k--', 'linewidth', .7)
% checked

% New, like _computed3 but redefining g: I redefined u, v so that input (t-3) now becomes (t-3)/2, etc
g = @(r,u,v) (2*(acos(u./r)-asin(v./r)).*u.*v + r.^2 + u.^2 + v.^2 - 2*u.*sqrt(r.^2-v.^2) - 2*v.*sqrt(r.^2-u.^2))/pi;
ind_odd = mod(funt_fine,2)==1;
ind_even = ~ind_odd;
result_computed5 = NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_fine(ind_fine);
    potential_i = 0:t-3;
    potential_j = t-1-potential_i;
    ind_actual_ij = (potential_i-1).^2 + (potential_j-1).^2 < L^2;
    actual_i = potential_i(ind_actual_ij);
    actual_j = potential_j(ind_actual_ij);
    res = 0;
    for ind_ij = 1:numel(actual_i)
        ii = actual_i(ind_ij);
        jj = actual_j(ind_ij);
        res = res + g(L, ii-1, jj-1);
    end
    result_computed5(ind_fine) = res;
end
plot(L_fine.*ind_odd./ind_odd, result_computed5.*ind_odd./ind_odd, 'k--', 'linewidth', .7), hold on
plot(L_fine.*ind_even./ind_even, result_computed5.*ind_even./ind_even, 'k--', 'linewidth', .7)
% checked


%% #86: like #64 but changing the segment shift, and with minor improvements. Figure for proof of probmax: funt odd
clf
hold on
set(gca, 'layer', 'top')
axis equal
grid on
plot([-.4 5], [0 0], 'k')
plot([-.4j 5j], 'k')
set(gca, 'GridAlpha', .4)
xticks(-1:5), yticks(-1:5)
xticklabels({'' '0' '1' ['i' char(hex2dec('2212')) '1'] '' '' ''})
yticklabels({'' '0' '1' '' '' ['j' char(hex2dec('2212')) '1'] ''})
xlabel x 
ylabel y
axis([-1 5 -1 5])
L = 3.52;
theta_0 = asin(3/L);
theta_1 = acos(1/L);
theta = theta_0 + linspace(-10*pi/180,0,50);
plot(1+1j + L*exp(1j*theta), ':', 'linewidth', .75)
theta = theta_1 + linspace(0,10*pi/180,50);
set(gca, 'colororderindex', 1)
plot(1+1j + L*exp(1j*theta), ':', 'linewidth', .75)
theta = linspace(theta_0,theta_1,50);
set(gca, 'colororderindex', 1)
plot(1+1j + L*exp(1j*theta), '-', 'linewidth', .75)
theta_example = theta_0 + (theta_1-theta_0)*.55;
endpoint = 1+1j + L*exp(1j*(theta_example));
patch([real(endpoint) 2 2 real(endpoint) real(endpoint)], [imag(endpoint) imag(endpoint) 4 4 imag(endpoint)], .88*[1 1 1], 'edgecolor', 'none')
plot([endpoint real(endpoint)+4j], 'k:')
plot([endpoint 2+1j*imag(endpoint)], 'k:')
plot([1+1j endpoint], 'k')
for k = [0 2+4j]
    if k==0
        plot(real(k)+[0 1 1 0 0], imag(k)+[0 0 1 1 0], 'k-', 'linewidth', 1)
    else
        plot(real(k)+[0 1 1 0 0], imag(k)+[0 0 1 1 0], 'k-', 'linewidth', 1)
    end
end
plot([1+1j], 'k.', 'markersize', 7)
plot(endpoint, 'k.', 'markersize', 7)
text(1.35, 2.65, char(hex2dec('2113')))
plot([1+1j 1.55+1j], 'k:')
theta = linspace(0,theta_example);
plot(1+1j + .3*exp(1j*theta), 'k-')
text(1.37, 1.26, char(hex2dec('03B8')))


%% #87: Figure for number of terms in the sum of probmax
clf
ms1 = 14;
ms2 = 10;
hold on
set(gca, 'layer', 'top')
axis equal
grid on
set(gca, 'GridAlpha', .4)
xticks(0:8), yticks(0:8)
%xticks(1:10), yticks(1:10)
xlabel x 
ylabel y
axis([min(xticks) max(xticks) min(yticks) max(yticks)])
plot([0 1 1 0 0], [0 0 1 1 0], 'k-', 'linewidth', 1) % reference tile
%plot([1+1j], 'k.', 'markersize', ms2) % dot at the corner of reference tile
t = floor(sqrt(2*ceil(L^2)-2))+3;
for L = [4.48 4.98 5.54 7.94 8.32 8.6 8.89 9.183]-2*sqrt(2) %[7.9 8.24 8.59 8.85 9.12]-sqrt(2)
    set(gca, 'colororderindex', 1)
    t = floor(sqrt(2*ceil(L^2)-2))+3;
    ii = 1:t; % corners are at (i-1,j-1)
    jj = t+1-ii;
    ind = (ii-2).^2+(jj-2).^2<L^2;
    ii = ii(ind);
    jj = jj(ind);
    ind_extra = abs(ii-jj)>1;
    theta_aux = acos((t-3)/2*sqrt(2)/L);
    theta_0 = pi/4 - theta_aux;
    theta_1 = pi/4 + theta_aux;
    theta = linspace(theta_0, theta_1, 100);
    plot(1+1j + L*exp(1j*theta), '-', 'linewidth', .5)
    plot([(t-1)*1j t-1], '--', 'Color', 0*[1 1 1])
    plot(ii(~ind_extra)-1+1j*(jj(~ind_extra)-1), 'k.', 'markersize', 14) % dot at the corner of other tile
    plot(ii(ind_extra)-1+1j*(jj(ind_extra)-1), 'ko', 'markersize', 4) % dot at the corner of other, extra tile
end
text(8.162, 1.04, 't = 10')
text(8.162, 2.04, 't = 11')
text(4.8, .5, 't = 6')
text(2.66, .5, 't = 5')


%% #88: Figure for computing the number of terms in the sum of probmax; t even
clf
ms1 = 14;
ms2 = 10;
hold on
set(gca, 'layer', 'top')
axis equal
grid on
set(gca, 'GridAlpha', .4)
xticks(0:8), yticks(0:8)
%xticks(1:10), yticks(1:10)
xlabel x 
ylabel y
axis([min(xticks) max(xticks) min(yticks) max(yticks)])
plot([0 1 1 0 0], [0 0 1 1 0], 'k-', 'linewidth', 1) % reference tile
%plot([1+1j], 'k.', 'markersize', ms2) % dot at the corner of reference tile
for L = 8.32-2*sqrt(2) %[7.9 8.24 8.59 8.85 9.12]-sqrt(2)
    set(gca, 'colororderindex', 1)
    t = floor(sqrt(2*ceil(L^2)-2))+3;
    ii = 1:t; % corners are at (i-1,j-1)
    jj = t+1-ii;
    ind = (ii-2).^2+(jj-2).^2<L^2;
    ii = ii(ind);
    jj = jj(ind);
    ind_extra = abs(ii-jj)>1;
    for n = 1:numel(ii)
        patch(ii(n)-1+[0 1 1 0 0], jj(n)-1+[0 0 1 1 0], .88*[1 1 1], 'edgecolor', 'none')
        lw = .5;
        plot(ii(n)-1 + [0 1], jj(n)-1 + [0 0], 'k-', 'linewidth', lw) % horizontal sides of tiles
        plot(ii(n)-1 + [0 0], jj(n)-1 + [0 1], 'k-', 'linewidth', lw) % vertical sides of tiles
        plot(ii(n)-1 + [0 1], jj(n)   + [0 0], 'k-', 'linewidth', lw) % farther horizontal sides
        plot(ii(n)   + [0 0], jj(n)-1 + [0 1], 'k-', 'linewidth', lw) % farther vertical sides
    end
    theta_aux = acos((t-3)/2*sqrt(2)/L);
    theta_0 = pi/4 - theta_aux;
    theta_1 = pi/4 + theta_aux;
    theta = linspace(theta_0, theta_1, 100);
    plot(1+1j + L*exp(1j*theta), '-', 'linewidth', .75) % arc
    plot([(t-1)*1j t-1], '--', 'Color', 0*[1 1 1]) % diagonal line, -45º slope
    plot(ii(~ind_extra)-1+1j*(jj(~ind_extra)-1), 'k.', 'markersize', 14) % dot at the corner of other tile
    plot(ii(ind_extra)-1+1j*(jj(ind_extra)-1), 'ko', 'markersize', 4) % dot at the corner of other, extra tile
    plot([1+1j 1+1j+L*exp(1j*theta_0)], 'k-') % radius labelled "l"
    plot([1+1j (t-1)/2*(1+1j)], 'k-') % line with 45º slope
    plot([(t-1)/2*(1+1j) 1+1j+L*exp(1j*theta_0)], 'k-', 'linewidth', 1.25) % segment labelled "w"
end
text(8.162, 1.04, 't = 10')
text(3.47, 1.7, char(hex2dec('2113')))
text(5.06, 3.54, 'w')
text(2, 4.54, ['((t' char(hex2dec('2212')) '1)/2,(t' char(hex2dec('2212')) '1)/2)'])
set(gcf, 'Position', [680 558 560 400])


%% #89: Figure for computing the number of terms in the sum of probmax; t odd
clf
ms1 = 14;
ms2 = 10;
hold on
set(gca, 'layer', 'top')
axis equal
grid on
set(gca, 'GridAlpha', .4)
xticks(0:8), yticks(0:8)
%xticks(1:10), yticks(1:10)
xlabel x 
ylabel y
axis([min(xticks) max(xticks) min(yticks) max(yticks)])
plot([0 1 1 0 0], [0 0 1 1 0], 'k-', 'linewidth', 1) % reference tile
%plot([1+1j], 'k.', 'markersize', ms2) % dot at the corner of reference tile
for L = 8.88-2*sqrt(2) %[7.9 8.24 8.59 8.85 9.12]-sqrt(2)
    set(gca, 'colororderindex', 1)
    t = floor(sqrt(2*ceil(L^2)-2))+3;
    ii = 1:t; % corners are at (i-1,j-1)
    jj = t+1-ii;
    ind = (ii-2).^2+(jj-2).^2<L^2+4; % añado 4 para que haya otro par de puntos
    ii = ii(ind);
    jj = jj(ind);
    ind_extra = abs(ii-jj)>1;
    for n = 1:numel(ii)
        if (ii(n)-2).^2+(jj(n)-2).^2<L^2 % only plot tiles that not covered
            patch(ii(n)-1+[0 1 1 0 0], jj(n)-1+[0 0 1 1 0], .88*[1 1 1], 'edgecolor', 'none')
            lw = .5;
            plot(ii(n)-1 + [0 1], jj(n)-1 + [0 0], 'k-', 'linewidth', lw) % horizontal sides of tiles
            plot(ii(n)-1 + [0 0], jj(n)-1 + [0 1], 'k-', 'linewidth', lw) % vertical sides of tiles
            plot(ii(n)-1 + [0 1], jj(n)   + [0 0], 'k-', 'linewidth', lw) % farther horizontal sides
            plot(ii(n)   + [0 0], jj(n)-1 + [0 1], 'k-', 'linewidth', lw) % farther vertical sides
        end
    end
    theta_aux = acos((t-3)/2*sqrt(2)/L);
    theta_0 = pi/4 - theta_aux;
    theta_1 = pi/4 + theta_aux;
    theta = linspace(theta_0, theta_1, 100);
    plot(1+1j + L*exp(1j*theta), '-', 'linewidth', .75) % arc
    plot([(t-1)*1j t-1], '--', 'Color', 0*[1 1 1]) % diagonal line, -45º slope
    plot(ii(~ind_extra)-1+1j*(jj(~ind_extra)-1), 'k.', 'markersize', 14) % dot at the corner of other tile
    plot(ii(ind_extra)-1+1j*(jj(ind_extra)-1), 'ko', 'markersize', 4) % dot at the corner of other, extra tile
    plot([1+1j 1+1j+L*exp(1j*theta_0)], 'k-') % radius labelled "l"
    plot([1+1j (t-1)/2*(1+1j)], 'k-') % line with 45º slope
    plot([(t-1)/2*(1+1j) 1+1j+L*exp(1j*theta_0)], 'k-', 'linewidth', 1.25) % segment labelled "w"
end
text(8.162, 2.04, 't = 11')
text(3.6, 1.96, char(hex2dec('2113')))
text(5.25, 4.31, 'w')
text(2.47, 5.03, ['((t' char(hex2dec('2212')) '1)/2,(t' char(hex2dec('2212')) '1)/2)'])
set(gcf, 'Position', [680 558 560 400])


%% #90. Figure: probmax, incorporating the correction from #85
% I use #85 computed4, additionally changing the order of second and third inputs of g

clear, clc
L_fine = .001:.001:16; % spacing should be small because curves are almost vertical at some points
ii = ceil(L_fine/sqrt(2)) + 1;
jj = ceil(sqrt(L_fine.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths

g = @(r,u,v) (2*(acos(v./r)-asin(u./r)).*u.*v + r.^2 + u.^2 + v.^2 - 2*u.*sqrt(r.^2-v.^2) - 2*v.*sqrt(r.^2-u.^2))/pi;
ind_odd = mod(funt_fine,2)==1;
ind_even = ~ind_odd;
num_k2 = ceil( sqrt(L_fine.^2/2-(funt_fine-3).^2/4) - (1+ind_even/2) );
result_computed = NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_fine(ind_fine);
    if ind_odd(ind_fine)
        result_computed(ind_fine) = g(L, (t-3)/2, (t-3)/2); % this term is always present
        for k = 1:num_k2(ind_fine) % extra terms
            result_computed(ind_fine) = result_computed(ind_fine) + 2*g(L, (t-3)/2-k, (t-3)/2+k);
        end
    else
        result_computed(ind_fine) = 2*g(L, (t-4)/2, (t-2)/2); % this term is always present
        for k = 1:num_k2(ind_fine) % extra terms
            result_computed(ind_fine) = result_computed(ind_fine) + 2*g(L, (t-4)/2-k, (t-2)/2+k);
        end
    end
end

hold on, grid on, box on, set(gca, 'colororderindex', 2)
plot(L_fine.*ind_odd./ind_odd, result_computed.*ind_odd./ind_odd, '-')
set(gca, 'colororderindex', 1)
plot(L_fine.*ind_even./ind_even, result_computed.*ind_even./ind_even, '-')
%plot(L_fine, result_computed, '-')
legend({'Odd maximum number of tiles' 'Even maximum number of tiles'})
xlabel(char(hex2dec('2113'))), ylabel([char(hex2dec('03C1')) '(' char(hex2dec('2113')) ')'])
set(gcf, 'Position', [680 558 560 380])


%% #91. sqrt(t)*probmax(l_t)

clear, clf
t_all = 10:1000;
ind = mod(t_all,2)==1;
L_all(ind) = hypot(t_all(ind)-3,t_all(ind)-1)/2;
ind = mod(t_all,2)==0;
L_all(ind) = (t_all(ind)-2)/sqrt(2);
L_all = L_all - 1e-6;
ii = ceil(L_all/sqrt(2)) + 1;
jj = ceil(sqrt(L_all.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths

g = @(r,u,v) (2*(acos(v./r)-asin(u./r)).*u.*v + r.^2 + u.^2 + v.^2 - 2*u.*sqrt(r.^2-v.^2) - 2*v.*sqrt(r.^2-u.^2))/pi;
ind_odd = mod(funt_fine,2)==1;
ind_even = ~ind_odd;
num_k2 = ceil( sqrt(L_all.^2/2-(funt_fine-3).^2/4) - (1+ind_even/2) );
result_computed = NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_all(ind_fine);
    if ind_odd(ind_fine)
        result_computed(ind_fine) = g(L, (t-3)/2, (t-3)/2); % this term is always present
        for k = 1:num_k2(ind_fine) % extra terms
            result_computed(ind_fine) = result_computed(ind_fine) + 2*g(L, (t-3)/2-k, (t-3)/2+k);
        end
    else
        result_computed(ind_fine) = 2*g(L, (t-4)/2, (t-2)/2); % this term is always present
        for k = 1:num_k2(ind_fine) % extra terms
            result_computed(ind_fine) = result_computed(ind_fine) + 2*g(L, (t-4)/2-k, (t-2)/2+k);
        end
    end
end

marker = '.'; ms = 5;
hold on, grid on, box on, set(gca, 'colororderindex', 2)
plot(t_all.*ind_odd./ind_odd, result_computed.*sqrt(t_all).*ind_odd./ind_odd, marker, 'markersize', ms)
set(gca, 'colororderindex', 1)
plot(t_all.*ind_even./ind_even, result_computed.*sqrt(t_all).*ind_even./ind_even, marker, 'markersize', ms)
legend({'Odd t' 'Even t'})
xlabel('t')
ylabel(['t^{1/2}' char(hex2dec('03C1')) '(' char(hex2dec('2113')) '_t)'])
set(gcf, 'Position', [680 558 560 380])
%set(gca, 'XScale', 'log', 'YScale', 'log')
set(gca, 'colororderindex', 3)
%plot(t_all, 32*sqrt(2)/105/pi)
set(gcf, 'Position', [680 558 560 380])


%% #92. Figure for proof of asymptotic probmax

clear, clf, clc
t = 9; % 10
if mod(t,2)==1
    L = hypot(t-3,t-1)/2;
else
    L = (t-2)/sqrt(2);
end
y_aux = (t-3)/sqrt(2);
y0 = L-y_aux;
x0 = sqrt(L^2-y_aux^2);
theta_aux = acos(y_aux/L);
theta = linspace(pi/2-theta_aux, pi/2+theta_aux, 200);
plot(-y_aux*1j+L*exp(1j*theta), 'linewidth', 0.5)
hold on, axis equal
x_dots = (-t+(1-mod(t,2))/2:t)*sqrt(2);
x_dots = x_dots(abs(x_dots)<=x0);
xlim([-3.2 3.2]), ylim([-4.6 1.09])
plot(xlim, [0 0], 'k--')
plot([0 0], ylim, 'k--')
plot([0-1j*y_aux x0+0j], 'k-')
text(1.03, -1.96, [char(hex2dec('2113')) '_t'])
%xt = x_dots;
%xt = [xt x0]; [xt, ind] = sort(xt);
%xtl = cellfun(@(s) [s char(hex2dec('221A')) '2'], {[char(hex2dec('2212')) '1.5'], [char(hex2dec('2212')) '0.5'], '0.5', '1.5'}, 'UniformOutput', false);
%xtl = [xtl 'x_0']; %xrl = xtl(ind);
xt = [-x0 -sqrt(2) 0 sqrt(2) x0]; xticks(xt)
%xtl = {[ hex2dec('2212') 'w_t'] 'x_{t,n}' '0' '' 'w_t'};
xtl = {[ hex2dec('2212') 'w_t'] 'x_{t,n}' '0' [hex2dec('221A') '2/2'] 'w_t'};
xticks(xt)
xticklabels(xtl)
%yt = [-y_aux y0]; ytl = {['h' hex2dec('2212') hex2dec('2113')], 'h'};
yt = [-y_aux 0 y0];
ytl = {[hex2dec('2212') '(t' hex2dec('2212') '3)/' char(hex2dec('221A')) '2'] '0' 'h_t'};
yticks(yt)
yticklabels(ytl)
ind_dot = 1; % select a tile
theta = 1.815; % pi/2*1.06; % choose manually, within the selected tile
z0 = -1j*y_aux+L*exp(1j*theta);
m = (x_dots(ind_dot)+z0)/2;
d = (x_dots(ind_dot)-z0)/2;
z1 = exp(-2j*(angle(d)-5/4*pi))*d+m;
z2 = -exp(-2j*(angle(d)-5/4*pi))*d+m;
patch(real([z0 z1 x_dots(ind_dot) z2 z0]), imag([z0 z1 x_dots(ind_dot) z2 z0]), .88*[1 1 1], 'edgecolor', 'none')
plot([z0 z1], 'k:')
plot([z0 z2], 'k:')
for x_dots_each = x_dots
    if abs(x_dots_each)<.8
        plot(x_dots_each, 0, 'k.', 'markersize', 14)
    else
        plot(x_dots_each, 0, 'ko', 'markersize', 4)
    end
    
end
%text(-1,-5, ['h' hex2dec('2212') hex2dec('2113')])
set(gca, 'colororderindex', 1)
%theta = linspace(pi/2, pi/2+.2264, 20); % manually
theta = linspace(1.7127, 1.9447, 50); % manually
plot(-y_aux*1j+L*exp(1j*theta), 'linewidth', 1.25)
for ind = 1:3
    plot(x_dots(ind)+[-sqrt(.5) 0], [sqrt(.5) 0], 'k-', 'linewidth', .5)
    plot(x_dots(ind)+[sqrt(.5) 0], [sqrt(.5) 0], 'k-', 'linewidth', .5)
end
set(gca, 'colororderindex', 1)
plot([-1j*y_aux z0], '-', 'linewidth', .5)
theta = linspace(0,angle(z0+1j*y_aux), 80);
plot(-1j*y_aux + .35*exp(1j*theta), '-')
text(-.35, -3.88, char(hex2dec('03B8')))
text(-1.86, .77, 'c_{t,n}')
plot(-1j*y_aux+[0 .55], 'k:')
s = .5825; plot(x_dots(ind_dot)+[-sqrt(.5) 0]*s, [sqrt(.5) 0]*s, 'k-', 'linewidth', 1) % s manually
s = 1; plot(x_dots(ind_dot)+[sqrt(.5) 0], [sqrt(.5) 0], 'k-', 'linewidth', 1) % s manually
grid on
xlabel x, ylabel y
box off
set(gcf, 'Position', [680 558 560 390]) % smaller figure, so font is larger in document
%figure(gcf)


%% #93. Detail of #92, single term

clear, clf, clc
t = 9; % 10
if mod(t,2)==1
    L = hypot(t-3,t-1)/2;
else
    L = (t-2)/sqrt(2);
end
y_aux = (t-3)/sqrt(2);
y0 = L-y_aux;
x0 = sqrt(L^2-y_aux^2);
hold on, axis equal
x_dots = (-t+(1-mod(t,2))/2:t)*sqrt(2);
x_dots = x_dots(abs(x_dots)<=x0);
ind_dot = 1; % select a tile
xlim([-2.29 -.285])
ylim([0 .77])
%plot(xlim, min(ylim)*[1 1], 'color', 'w', 'linewidth', 1.5)
%plot(min(xlim)*[1 1], ylim, 'color', 'w', 'linewidth', 1.5)
%plot(xlim, [0 0], 'k--')
theta = 1.815; % pi/2*1.06; % choose manually, within the selected tile
z0 = -1j*y_aux+L*exp(1j*theta);
z0 = real(z0) + 1j*.5747; % manually
m = (x_dots(ind_dot)+z0)/2;
d = (x_dots(ind_dot)-z0)/2;
z1 = exp(-2j*(angle(d)-5/4*pi))*d+m;
z2 = -exp(-2j*(angle(d)-5/4*pi))*d+m;
patch(real([z0 z1 x_dots(ind_dot) z2 z0]), imag([z0 z1 x_dots(ind_dot) z2 z0]), .88*[1 1 1], 'edgecolor', 'none')
if abs(x_dots(ind_dot))<.8
    plot(x_dots(ind_dot), 0, 'k.', 'markersize', 14)
else
    plot(x_dots(ind_dot), 0, 'ko', 'markersize', 4)
end
set(gca, 'colororderindex', 1)
%theta = linspace(pi/2, pi/2+.2264, 20); % manually
theta = linspace(1.7127, 1.9447, 40); % manually
plot(-y_aux*1j+L*exp(1j*theta), 'linewidth', 1)
theta_1 = linspace(theta(1)-.043, theta(1), 20); % manually
set(gca, 'colororderindex', 1)
plot(-y_aux*1j+L*exp(1j*theta_1), 'linewidth', .5)
theta_2 = linspace(theta(end), theta(end)+.050, 20); % manually
set(gca, 'colororderindex', 1)
plot(-y_aux*1j+L*exp(1j*theta_2), 'linewidth', .5)
plot(x_dots(ind_dot)+[-sqrt(.5) 0], [sqrt(.5) 0], 'k-', 'linewidth', .5)
plot(x_dots(ind_dot)+[sqrt(.5) 0], [sqrt(.5) 0], 'k-', 'linewidth', .5)
s1 = .5825; plot(x_dots(ind_dot)+[-sqrt(.5) 0]*s1, [sqrt(.5) 0]*s1, 'k-', 'linewidth', 1) % s manually
s2 = 1; plot(x_dots(ind_dot)+[sqrt(.5) 0]*s2, [sqrt(.5) 0]*s2, 'k-', 'linewidth', 1) % s manually
set(gca, 'colororderindex', 1)
plot([x_dots(ind_dot)-sqrt(.5)*s1 x_dots(ind_dot)+sqrt(.5)*s2], [sqrt(.5)*s1 sqrt(.5)*s2], 'linewidth', 1)
set(gca, 'colororderindex', 1)
plot([x_dots(ind_dot) x_dots(ind_dot)+1j*.520557], '-', 'linewidth', .5) % manually
plot([z0 z1], 'k:')
plot([z0 z2], 'k:')
%grid on
xticks([-1.5 -1 -.5]*sqrt(2))
xticklabels({['x_{t,n}' hex2dec('2212') hex2dec('221A') '2/2'] 'x_{t,n}' ['x_{t,n}' '+' hex2dec('221A') '2/2']})
yticks([0 .5*sqrt(2)])
yticklabels({'0' [hex2dec('221A') '2/2']})
%yticks([])
%xlabel x, ylabel y
%set(gca, 'XColor', 'none', 'Ycolor', 'none')
%figure(gcf)
text(-1.79,.24,'a_{t,n}')
text(-.945,.4,'b_{t,n}')
text(-1.39,.665,'c_{t,n}')
text(-1.533,.437,'y_{t,n}')
set(gcf, 'Position', [680 558 500 300]) % smaller figure, so font is larger in document
%grid on


%% #93bis. Like #93 but including tangent and some other changes

clear, clf, clc
t = 9; % 10
if mod(t,2)==1
    L = hypot(t-3,t-1)/2;
else
    L = (t-2)/sqrt(2);
end
y_aux = (t-3)/sqrt(2);
y0 = L-y_aux;
x0 = sqrt(L^2-y_aux^2);
hold on, axis equal
x_dots = (-t+(1-mod(t,2))/2:t)*sqrt(2);
x_dots = x_dots(abs(x_dots)<=x0);
ind_dot = 1; % select a tile
xlim([-2.29 -.285])
ylim([0 .77])
%plot(xlim, min(ylim)*[1 1], 'color', 'w', 'linewidth', 1.5)
%plot(min(xlim)*[1 1], ylim, 'color', 'w', 'linewidth', 1.5)
%plot(xlim, [0 0], 'k--')
theta = 1.815; % pi/2*1.06; % choose manually, within the selected tile
z0 = -1j*y_aux+L*exp(1j*theta);
z0 = real(z0) + 1j*.5747; % manually
m = (x_dots(ind_dot)+z0)/2;
d = (x_dots(ind_dot)-z0)/2;
z1 = exp(-2j*(angle(d)-5/4*pi))*d+m;
z2 = -exp(-2j*(angle(d)-5/4*pi))*d+m;
patch(real([z0 z1 x_dots(ind_dot) z2 z0]), imag([z0 z1 x_dots(ind_dot) z2 z0]), .88*[1 1 1], 'edgecolor', 'none')
if abs(x_dots(ind_dot))<.8
    plot(x_dots(ind_dot), 0, 'k.', 'markersize', 14)
else
    plot(x_dots(ind_dot), 0, 'ko', 'markersize', 4)
end
set(gca, 'colororderindex', 1)
%theta = linspace(pi/2, pi/2+.2264, 20); % manually
theta = linspace(1.7127, 1.9447, 100); % manually
plot(-y_aux*1j+L*exp(1j*theta), 'linewidth', 1)
theta_1 = linspace(theta(1)-.043, theta(1), 20); % manually
set(gca, 'colororderindex', 1)
plot(-y_aux*1j+L*exp(1j*theta_1), 'linewidth', .5)
theta_2 = linspace(theta(end), theta(end)+.050, 20); % manually
set(gca, 'colororderindex', 1)
plot(-y_aux*1j+L*exp(1j*theta_2), 'linewidth', .5)
plot(-y_aux*1j+L*( exp(1j*theta(end)) + 35*[0 exp(1j*theta_2(1))-exp(1j*theta_2(2))] ), 'k--', 'linewidth', .5) % tangent. Length manually
plot(x_dots(ind_dot)+[-sqrt(.5) 0], [sqrt(.5) 0], 'k-', 'linewidth', .5)
plot(x_dots(ind_dot)+[sqrt(.5) 0], [sqrt(.5) 0], 'k-', 'linewidth', .5)
s1 = .5825; plot(x_dots(ind_dot)+[-sqrt(.5) 0]*s1, [sqrt(.5) 0]*s1, 'k-', 'linewidth', 1) % s manually
s2 = 1; plot(x_dots(ind_dot)+[sqrt(.5) 0]*s2, [sqrt(.5) 0]*s2, 'k-', 'linewidth', 1) % s manually
set(gca, 'colororderindex', 1)
plot([x_dots(ind_dot)-sqrt(.5)*s1 x_dots(ind_dot)+sqrt(.5)*s2], [sqrt(.5)*s1 sqrt(.5)*s2], 'linewidth', 1)
set(gca, 'colororderindex', 1)
plot([x_dots(ind_dot) x_dots(ind_dot)+1j*.520557], '-', 'linewidth', .5) % manually
plot([z0 z1], 'k:')
plot([z0 z2], 'k:')
%grid on
xticks([-1.5 -1 -.5]*sqrt(2))
xticklabels({['x_{t,n}' hex2dec('2212') hex2dec('221A') '2/2'] 'x_{t,n}' ['x_{t,n}' '+' hex2dec('221A') '2/2']})
yticks([0 .5*sqrt(2)])
yticklabels({'0' [hex2dec('221A') '2/2']})
%yticks([])
xlabel x, ylabel y
%set(gca, 'XColor', 'none', 'Ycolor', 'none')
%figure(gcf)
text(-1.79,.24,'a_{t,n}')
text(-.945,.4,'b_{t,n}')
text(-1.17,.715,'c_{t,n}')
%text(-1.533,.437,'y_{t,n}')
%text(-1.533,.436,[hex2dec('1EF9') '_{t,n}']) % "y" with tilde
text(-1.532,.436,['y_{t,n}']), text(-1.532,.4682,'~','fontsize',8)
set(gcf, 'Position', [680 558 500 300]) % smaller figure, so font is larger in document
text(-1.053,.5658,'c_{t,n}'), text(-1.053,.598,'~','fontsize',8) % there is no "c" with tilde
%grid on


%% #93ter. Like #93bis but including some changes

clear, clf, clc
t = 9; % 10
if mod(t,2)==1
    L = hypot(t-3,t-1)/2;
else
    L = (t-2)/sqrt(2);
end
y_aux = (t-3)/sqrt(2);
y0 = L-y_aux;
x0 = sqrt(L^2-y_aux^2);
hold on, axis equal
x_dots = (-t+(1-mod(t,2))/2:t)*sqrt(2);
x_dots = x_dots(abs(x_dots)<=x0);
ind_dot = 1; % select a tile
xlim([-2.29 -.295])
ylim([0 .77])
%plot(xlim, min(ylim)*[1 1], 'color', 'w', 'linewidth', 1.5)
%plot(min(xlim)*[1 1], ylim, 'color', 'w', 'linewidth', 1.5)
%plot(xlim, [0 0], 'k--')
theta = 1.815; % pi/2*1.06; % choose manually, within the selected tile
z0 = -1j*y_aux+L*exp(1j*theta);
z0 = real(z0) + 1j*.5747; % manually
m = (x_dots(ind_dot)+z0)/2;
d = (x_dots(ind_dot)-z0)/2;
z1 = exp(-2j*(angle(d)-5/4*pi))*d+m;
z2 = -exp(-2j*(angle(d)-5/4*pi))*d+m;
patch(real([z0 z1 x_dots(ind_dot) z2 z0]), imag([z0 z1 x_dots(ind_dot) z2 z0]), .88*[1 1 1], 'edgecolor', 'none')
if abs(x_dots(ind_dot))<.8
    plot(x_dots(ind_dot), 0, 'k.', 'markersize', 14)
else
    plot(x_dots(ind_dot), 0, 'ko', 'markersize', 4)
end
set(gca, 'colororderindex', 1)
%theta = linspace(pi/2, pi/2+.2264, 20); % manually
theta = linspace(1.7127, 1.9447, 100); % manually
plot(-y_aux*1j+L*exp(1j*theta), 'linewidth', 1)
theta_1 = linspace(theta(1)-.043, theta(1), 20); % manually
set(gca, 'colororderindex', 1)
plot(-y_aux*1j+L*exp(1j*theta_1), 'linewidth', .5)
theta_2 = linspace(theta(end), theta(end)+.050, 20); % manually
set(gca, 'colororderindex', 1)
plot(-y_aux*1j+L*exp(1j*theta_2), 'linewidth', .5)
plot(-y_aux*1j+L*( exp(1j*theta(end)) + 38*[0 exp(1j*theta_2(1))-exp(1j*theta_2(2))] ), 'k--', 'linewidth', .5) % tangent. Length manually
plot(x_dots(ind_dot)+[-sqrt(.5) 0], [sqrt(.5) 0], 'k-', 'linewidth', .5)
plot(x_dots(ind_dot)+[sqrt(.5) 0], [sqrt(.5) 0], 'k-', 'linewidth', .5)
s1 = .5825; plot(x_dots(ind_dot)+[-sqrt(.5) 0]*s1, [sqrt(.5) 0]*s1, 'k-', 'linewidth', 1) % s manually
s2 = 1; plot(x_dots(ind_dot)+[sqrt(.5) 0]*s2, [sqrt(.5) 0]*s2, 'k-', 'linewidth', 1) % s manually
set(gca, 'colororderindex', 1)
plot([x_dots(ind_dot)-sqrt(.5)*s1 x_dots(ind_dot)+sqrt(.5)*s2], [sqrt(.5)*s1 sqrt(.5)*s2], 'linewidth', 1)
plot([x_dots(ind_dot)-sqrt(.5)*s1 x_dots(ind_dot)+sqrt(.5)*s2], [0 0], 'k--', 'linewidth', 1.2) % projection
text(-1.205, .068 , 'I_{t,n}') % its text
set(gca, 'colororderindex', 1)
plot([x_dots(ind_dot) x_dots(ind_dot)+1j*.520557], 'k-', 'linewidth', .5) % manually
plot([z0 z1], 'k:')
plot([z0 z2], 'k:')
%grid on
xticks([-1.5 -1 -.5]*sqrt(2))
xticklabels({['x_{t,n}' hex2dec('2212') hex2dec('221A') '2/2'] 'x_{t,n}' ['x_{t,n}' '+' hex2dec('221A') '2/2']})
yticks([0 .520557 .553 .5*sqrt(2)])
yticklabels({'0' '' '' [hex2dec('221A') '2/2']})
%yticks([])
xlabel x, ylabel y
%set(gca, 'XColor', 'none', 'Ycolor', 'none')
%figure(gcf)
%text(-1.79,.254,'a_{t,n}')
text(-1.787,.254,[hex2dec('03B1') '_{t,n}'])
%text(-.945,.41,'b_{t,n}')
text(-.945,.41,[hex2dec('03B2') '_{t,n}'])
text(-1.16,.705,'c_{t,n}')
%text(-1.533,.437,'y_{t,n}')
%text(-1.533,.436,[hex2dec('1EF9') '_{t,n}']) % "y" with tilde
x = -2.41; y = .495; text(x,y,['y_{t,n}']), text(x+.001,y+.0322,'~','fontsize',8)
x = -2.41; y = .595; text(x,y,['y_{t,n}'])
set(gcf, 'Position', [280 258 550 330]) % smaller figure, so font is larger in document
x = -1.053; y = .5658; text(x,y,'c_{t,n}'), text(x+.001,y+.0322,'~','fontsize',8) % there is no "c" with tilde
grid on



%% #94. Scaled versions of the circle

clf, hold on
for t = 4:2:16 %5:2:17
    if mod(t,2)==1
        L = hypot(t-3,t-1)/2;
    else
        L = (t-2)/sqrt(2);
    end
    w = sqrt(L^2-(t-3)^2/2);
    y_aux = (t-3)/sqrt(2);
    y0 = L-y_aux;
    %z = linspace(-w,w,300)/sqrt(t); plot(z, y0-L+sqrt(L^2-t*z.^2))
    z = linspace(-w,w,300)/w; plot(z, y0-L+sqrt(L^2-w^2*z.^2))
end


%% #95. Check of derivative of Omega normalized with respect to t, for t odd

clear, clc, close all
t = 11;
lt = hypot(t-3, t-1)/2;
wt = sqrt(t-2);
z = linspace(-1,1,200);
Omega = -(t-3)/sqrt(2) + sqrt( ((t-2).^2+1)/2 - (t-2)*z.^2 );
%plot(z, Omega), hold on

y = -(t-3)/sqrt(2) + sqrt(lt^2-wt^2*z.^2);
%plot(z, y, '--') % checked

tp = t-2;
Omega2 = -(tp-1)/sqrt(2) + sqrt( (tp.^2+1)/2 - tp.*z.^2 );
%plot(z, Omega2, '--'), % checked

clear
t = linspace(11,21,1000);
z = .008;
tp = t-2;
clear t
Omega2 = -(tp-1)/sqrt(2) + sqrt( (tp.^2+1)/2 - tp.*z.^2 );
%plot(tp(1:end-1), diff(Omega2)./diff(tp)), hold on
%plot(tp, -1/sqrt(2) + (tp-z.^2)/2./sqrt( (tp.^2+1)/2 - tp.*z.^2 ), '--') % checked

%plot(tp-z.^2), hold on
%plot(sqrt(tp.^2+1-sqrt(2)+z.^2))

clear
tp = 7;
z = linspace(0,1,200);
%plot(z, (tp-z.^2)/2./sqrt( (tp.^2+1)/2 - tp.*z.^2))
y = (tp-z.^2)/2./sqrt( (tp.^2+1)/2 - tp.*z.^2);
%y = (tp-z.^2)/2./sqrt( (tp.^2+z)/2 - tp.*z.^2);
plot(z(1:end-1), diff(y)./diff(z))



%% #96: like #62 but only one tile and its symmetric tiles
clf
hold on
set(gca, 'layer', 'bottom')
axis equal
grid on
plot([-3.7 4.7], [0 0], 'k--', 'linewidth', .75) % horizontal axis
plot([-3.7j 4.7j], 'k--', 'linewidth', .75) % vertical axis
set(gca, 'GridAlpha', .4)
xticks(-3:4), yticks(-3:4)
%xticklabels([arrayfun(@(n)[num2str(n) char(8201) 'a'], -3:-1, 'UniformOutput', false) {'0'}  arrayfun(@(n)[num2str(n) char(8201) 'a'], 1:4, 'UniformOutput', false)]) % thin space
%yticklabels([arrayfun(@(n)[num2str(n) char(8201) 'a'], -3:-1, 'UniformOutput', false) {'0'}  arrayfun(@(n)[num2str(n) char(8201) 'a'], 1:4, 'UniformOutput', false)]) % thin space
xlabel x 
ylabel y
axis([-4 5 -4 5])
for k = [0 2+3j 2-3j -2+3j -2-3j] % [0 3+2j 3-2j -3+2j -3-2j]
    if k==0
        patch(real([k k+1 k+1 k k]), imag([k k k+1j k+1j k]), .88*[1 1 1])
        plot(real(k)+[0 1 1 0 0], imag(k)+[0 0 1 1 0], 'k-', 'linewidth', 1)
    else
        patch(real([k k+1 k+1 k k]), imag([k k k+1j k+1j k]), .88*[1 1 1])
        plot(real(k)+[0 1 1 0 0], imag(k)+[0 0 1 1 0], 'k--', 'linewidth', 1)
        plot(k, 'k.', 'markersize', 12)
    end
end
set(gcf, 'Position', [680 558 560 400])
h = 1.23; v = 2.7;
text(h,   v,   ['(i' hex2dec('2212') '1,j' hex2dec('2212') '1)'])
text(h-4, v,   ['(1' hex2dec('2212') 'i,j' hex2dec('2212') '1)'])
text(h,   v-6, ['(i' hex2dec('2212') '1,1' hex2dec('2212') 'j)'])
text(h-4, v-6, ['(1' hex2dec('2212') 'i,1' hex2dec('2212') 'j)'])


%% #97: Check for N_t

clear, close all
t_all = 3:1000;
num_terms = NaN(size(t_all));
num_terms2 = NaN(size(t_all));
for c = 1:numel(t_all)
    t = t_all(c);
    potential_k = -t:t; % more than enough
    if mod(t,2) % odd
        L = hypot(t-3,t-1)/2;
        ind_actual_k = ((t-3)/2-potential_k).^2 + ((t-3)/2+potential_k).^2 < L^2;
        num_terms2(c) = 2*ceil(sqrt((t-2)/2))-1;
    else % even
        L = (t-2)/sqrt(2);
        ind_actual_k = ((t-4)/2-potential_k).^2 + ((t-2)/2+potential_k).^2 < L^2;
        num_terms2(c) = 2*ceil(sqrt((t-2.5)/2)-0.5);
    end
    num_terms(c) = sum(ind_actual_k);
end
isequal(num_terms, num_terms2) % check


%% #98. Based on #85: figure for comparing computed and simulated probability

clear, close all
L_fine = 15:.001:20;

% From #77
g = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;
ii = ceil(L_fine/sqrt(2)) + 1;
jj = ceil(sqrt(L_fine.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
num_k =  NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_fine(ind_fine);
    if mod(t,2) % odd
        res = g(L, t-3, t-3); % this term is always present
        potential_k = 1:(t-3)/2;
        %%potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-3)/2-potential_k).^2 + ((t-3)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-3-2*k, t-3+2*k);
        end
    else % even
        res = 2*g(L, t-4, t-2); % this term is always present
        potential_k = 1:(t-4)/2;
        %%potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-4)/2-potential_k).^2 + ((t-2)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-4-2*k, t-2+2*k);
        end
    end
    result_computed(ind_fine) = res;
    num_k(ind_fine) = numel(actual_k);
end
ind_odd = mod(funt_fine,2)==1; % for plotting
ind_even = ~ind_odd;
hold on
set(gca,'ColorOrderIndex',2)
plot(L_fine.*ind_odd./ind_odd, result_computed.*ind_odd./ind_odd, '-', 'linewidth', .7), hold on, grid on
set(gca,'ColorOrderIndex',1)
plot(L_fine.*ind_even./ind_even, result_computed.*ind_even./ind_even, '-', 'linewidth', .7), hold on, grid on

% From #77
g = @(r,u,v) ((acos(u/2./r)-asin(v/2./r)).*u.*v + 2*r.^2 + u.^2/2 + v.^2/2 - u.*sqrt(4*r.^2-v.^2) - v.*sqrt(4*r.^2-u.^2))/2/pi;
ii = ceil(L_fine/sqrt(2)) + 1;
jj = ceil(sqrt(L_fine.^2 - (ii-2).^2 )) + 1;
funt_fine = ii+jj-1; % formula valid for real-valued lengths
result_computed = NaN(size(funt_fine));
num_k =  NaN(size(funt_fine));
for ind_fine = 1:numel(funt_fine)
    t = funt_fine(ind_fine);
    L = L_fine(ind_fine);
    if mod(t,2) % odd
        res = g(L, t-3, t-3); % this term is always present
        %potential_k = 1:(t-3)/2;
        potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-3)/2-potential_k).^2 + ((t-3)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-3-2*k, t-3+2*k);
        end
    else % even
        res = 2*g(L, t-4, t-2); % this term is always present
        %potential_k = 1:(t-4)/2;
        potential_k = 1; % old (wrong), that is, #65
        ind_actual_k = ((t-4)/2-potential_k).^2 + ((t-2)/2+potential_k).^2 < L^2;
        actual_k = potential_k(ind_actual_k);
        for k = actual_k
            res = res + 2*g(L, t-4-2*k, t-2+2*k);
        end
    end
    result_computed(ind_fine) = res;
    num_k(ind_fine) = numel(actual_k);
end
ind_odd = mod(funt_fine,2)==1; % for plotting
ind_even = ~ind_odd;
hold on
set(gca,'ColorOrderIndex',2)
plot(L_fine.*ind_odd./ind_odd, result_computed.*ind_odd./ind_odd, '--', 'linewidth', .7), hold on, grid on
set(gca,'ColorOrderIndex',1)
plot(L_fine.*ind_even./ind_even, result_computed.*ind_even./ind_even, '--', 'linewidth', .7), hold on, grid on

% Simulation, from #76
% I see the above computations match the simulation results
clear, clc
L = 15:.01:20; %3.8:.1:4.2;%3:.1:3.5;
N = 1e5; % elegir con cuidado: memoria. Mejor NN grande que N grande
NN = 100;
ii = ceil(L/sqrt(2)) + 1;
jj = ceil(sqrt(L.^2 - (ii-2).^2 )) + 1;
funt = ii+jj-1; % formula valid for real-valued lengths
result = 0;
for nn = 1:NN
    disp(nn)
    z0 = rand(N, numel(L)) + 1j*rand(N, numel(L));
    z1 = z0 + L.*exp(1j*2*pi*rand(N, numel(L)));
    ii = ceil(abs(real(z1))) + (real(z1)<0); % add 1 for negative
    jj = ceil(abs(imag(z1))) + (imag(z1)<0); % add 1 for negative
    result = result + mean(ii+jj-1 == funt , 1);
end
result = result/NN;
plot(L, result, '.', 'color', .5*[1 1 1], 'markersize', 7)

legend({'Odd maximum number of tiles', 'Even maximum number of tiles', 'Odd maximum number of tiles, incorrect', 'Even maximum number of tiles, incorrect', 'Simulated'})
xlabel(char(hex2dec('2113'))), ylabel([char(hex2dec('03C1')) '(' char(hex2dec('2113')) ')'])
set(gcf, 'Position', [360 175 560 445])
ylim([0 .04])