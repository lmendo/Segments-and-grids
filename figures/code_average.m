%% #1 Average number of touched squares, experimentally

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


%% #2 Selectable range for theta
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


%% #3 I define ki, kj signed (page 17, (*)), from which ii, jj are obtained 

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


%% #4 Plot of formula 2*L/pi*(1/a+1/b)+1

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




