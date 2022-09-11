%% Seting Impulse and Other parameters
clear;
n = 2^12; %number of grid point
twidth = 30; %width of time window [ps]
c = 299792458 * 1e9/1e12; %speed of light[nm/ps]
wavelength = 1900; %reference wavelength[nm]
T = linspace(-twidth/2, twidth/2, n); %time grid
power = 100;T1 = 0.1;
t0 = T1 / 1.665;   %impulse duration [ps]

A0 = sqrt(power) * exp(-T.^2/2/t0^2);  %impulse function
A0 = A0';
figure; plot(T,A0)
title('Impulse')

%% Find A0(w)

AW0 = fft(A0);  %colculate FFT of A0
dT = T(2) - T(1);   %grid parameter
V = 2 * pi * (-n/2:n/2-1)'/(n*dT);%frequency grid [2Pi*THZ]
w0 = (2.0 * pi * c) / wavelength; %reference frequency[2*pi*THz]
W = V + w0; %[2*pi*THz]
WL = (2 * pi * c)./ W ; %[nm]
WL = WL ./ 1000;    %[um] for N
figure; plot(V, fftshift(abs(AW0)), 'Color', 'red') %spectrum in frequency grid
title('Impulse FFT');
figure; plot(WL, fftshift(abs(AW0)), 'Color', 'green') %spectrum in wavelength grid
%% Passing through the medium
AW0 = fftshift(AW0);
N=sqrt(1+0.6961663./(1-(0.0684043./WL).^2)+0.4079426./(1-(0.1162414./WL).^2)+0.8974794./(1-(9.896161./WL).^2));
K = (W./ c).* N;    %[2Pi * THz * ps / nm] = [2Pi/nm]
L = linspace(0, 1000, 10);   %[mm]
a = exp(-1i .* K .* L .* 10^6);
AW = zeros(n,length(L));
for m = 1:n
    AW(m,:) = AW0(m) .* a(m,:);
end

%% IFFT A(w)->A(t)
A =  ifft(fftshift(AW,1));
%% Graphics

figure(); 
subplot(1,2,1);
pcolor(WL * 10^3, L, abs(AW')); % plot as pseudocolor map
xlim([1000,3000]); ylim([0,1000]);
xlabel('Wavelength / nm'); ylabel('Distance / mm');shading interp;

subplot(1,2,2);
pcolor(T', L, (abs(A'))); % plot as pseudocolor map
xlabel('Time / ps'); ylabel('Distance / mm'); shading interp;

%% Checking
figure; plot(T, abs(A(:,1))); hold on; plot(T, abs(A(:,end)))
GVD = -80.559; %[fs^2/mm]
Ld = (t0 * 1e3)^2 / abs(GVD);
T1 = t0*sqrt(1+(1000/Ld)^2); %anlitic impulse duration
figure;plot(T, abs(A(:,end).^2))
