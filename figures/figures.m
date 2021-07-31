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
    [ii, jj] = find((t + t.' <= L.^2) & (logical(t)==logical(t.'))); % con igualdad. No contamos las l�neas verticales u horizontales, pero s� el origen
    M = max(ii+jj)-1; % -2 por 0-based indexing; -1 por los cuadrados hasta (ii,jj); +2 por desplazamientos para coger trocitos
    [ii, jj] = find(t + t.' < L.^2);
    m = max(ii+jj)+1; % -2 por 0-based indexing; -1 por los cuadrados hasta (ii,jj); +4 por desplazamientos para coger trocitos
    result_all(end+1) = max(M,m); % M: -2 por 0-based indexing; -1 por los cuadrados hasta (ii,jj)
                                            % +4 � +2 por desplazamientos para coger trocitos. Si M y m son iguales: +4; si no (m==M-1 � M es vac�o): +2
end
% La parte de M nunca aporta nada (he comprobado hasta L=1000). Es porque realmente s�lo suma 1: el punto que se a�ade con "<=",
% el cual no estaba con "<". Pero los de al lado s� estaban: s�lo estamos sumando 1. Y luego se
% suman 2 en vez de 4: siempre es peor.


%% #4. Calculando. Pruebo f�rmula. �Funciona!
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


%% #6. Calculando. Pruebo (p�g. 5). Funciona
result_all = [];
for L = 1:1e5
    M = 1:2*L+1;
    result_all(end+1) = find(M.^2-6*M+10-mod(M,2) < 2*L^2, 1, 'last');
end


%% #7. Calculando. F�rmula de @xnor. �Funciona! He probado hasta 1e8
result_all = floor(sqrt(2*L.^2-2))+3; %


%% #8. Problema inverso
M = 3:1e3;
L_1 = floor(sqrt((M.^2-6*M+10-mod(M,2))/2))+1; % mi referencia (p�g. 3)
L_2 = ceil(sqrt((M-3).^2/2+1)); % "invirtiendo" f�rmula para M de @xnor. �Funciona!


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


%%  #47. Figures of the direct-problem and inverse problem sequences (funti and funli)
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