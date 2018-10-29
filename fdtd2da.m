clear all; close all; clc;

% lower left corner = (0, 0)

I = im2double(imread("problem/ant_140x200.png"));
lambda1 = 200;
deembedd = (0.1 + 0.653)*lambda1; % de-embedding distance
xc = (0.1 + 0.1 + 0.426)*lambda1; % phase center x coordinate
yc = rows(I)/2; % phase center y coordinate
nc = round(1920/2/2); % number of columns in the computation domain
nr = round(1080/2); % number of rows in the computation domain

%I = im2double(imread("problem/ant_350x500.png"));
%lambda1 = 500;
%deembedd = (0.1 + 0.653)*lambda1;
%xc = (0.1 + 0.1 + 0.426)*lambda1;
%yc = rows(I)/2;
%nc = round(1920/2);
%nr = round(1080);

%I = im2double(imread("problem/ant0_100x500.png"));
%lambda1 = 500;
%deembedd = (0.1 + 0.018)*lambda1;
%xc = (0.1 + 0.1 + 0.027)*lambda1;
%yc = rows(I)/2;
%nc = round(1920/2);
%nr = round(1080);

xl = round(xc - nc/2);
yb = round(yc - nr/2);
column = @(x) -xl + 1 + (x - 0.5);
row = @(y) nr + yb - (y - 0.5);
portx = round(lambda1/10);
porty1 = round(lambda1*0.45);
porty2 = round(lambda1*0.55);
portr = [row(porty2 - 0.5):row(porty1 + 0.5)];
portc = column(portx + 0.5);
threshold = 0.8;
[pecr, pecc] = find(mean(I, 3) < threshold);
pec = sub2ind([nr, nc], pecr + nr - rows(I) + yb, pecc - xl);

dt = 1/sqrt(2)*0.99;
tmax = 10*lambda1;
t = [0 : dt : tmax];
f1 = 1/lambda1;
df = f1/40;
fmax = f1*6;
f = [0 : df : fmax];
if0 = length(f);
f = [-fliplr(f(2 : end)), f];
nphi = 1080/3;
dphi = 1/nphi; % phi is measured in turns
phi = -0.5 + [1 : nphi]'*dphi;

pulseWidth = lambda1/8;
stimulus = @(t) exp(-((t - pulseWidth*4)/pulseWidth).^2);
Eyinc = @(t) sqrt(dt)/sqrt(abs(porty2 - porty1))*stimulus(t);

Hz = zeros(nr, nc);
Hzx = zeros(nr, nc);
Hzy = zeros(nr, nc);
Ex = zeros(nr + 1, nc);
Ey = zeros(nr, nc + 1);
a1 = zeros(size(t));
b1 = zeros(size(t));
Ic = zeros(size(phi)(1), size(t)(2));
epsx = ones(nr, nc); % permittivity
epsy = ones(nr, nc);
mux = ones(nr, nc); % permeability
muy = ones(nr, nc);
sigmax = zeros(nr, nc); % electric conductivity
sigmay = zeros(nr, nc);
taux = zeros(nr, nc); % magnetic conductivity
tauy = zeros(nr, nc);

pmln = 20;
pmldB = -50;
pmlfactor = 2;
pmlalpha = (pmlfactor + 1)*([1 : pmln]/pmln).^pmlfactor ...
                                                      *(-pmldB/40)*log(10)/pmln;
sigmax(:, end - pmln + 1 : end) = repmat(pmlalpha, nr, 1);
sigmax(:, 1 : pmln) = repmat(fliplr(pmlalpha), nr, 1);
sigmay(end - pmln + 1 : end, :) = repmat(pmlalpha', 1, nc);
sigmay(1 : pmln, :) = repmat(fliplr(pmlalpha)', 1, nc);
taux = sigmax;
tauy = sigmay;
sigmax(pec) = Inf;
sigmay(pec) = Inf;

circleR = floor((min(nc, nr) - 2*pmln)/2) - 1;
circleX = xc + circleR*cos(2*pi*phi);
circleY = yc + circleR*sin(2*pi*phi);
circle = sub2ind([nr, nc], round(row(circleY)), round(column(circleX)));

meanx = @(f) (f(:, 2 : end) + f(:, 1 : end - 1))/2;
meany = @(f) (f(2 : end, :) + f(1 : end - 1, :))/2;

c1Hzx = 1 - 1./(mux/dt./taux + 1/2);
c2Hzx = 1./(mux/dt + taux/2);
c1Hzy = 1 - 1./(muy/dt./tauy + 1/2);
c2Hzy = 1./(muy/dt + tauy/2);
c1Ex = 1 - 1./(meany(epsy)/dt./meany(sigmay) + 1/2);
c2Ex = 1./(meany(epsy)/dt + meany(sigmay)/2);
c1Ey = 1 - 1./(meanx(epsx)/dt./meanx(sigmax) + 1/2);
c2Ey = 1./(meanx(epsx)/dt + meanx(sigmax)/2);

outDir = "out";
mkdir(outDir);
downsampling = 5;
colorN = 256;
constraint = 0.2*max(abs(Eyinc(t)));
id = tic();
for i = 1 : length(t)
  printf("i = %d (%.1f%%)\n", i, i/length(t)*100)
  Hzx = Hzx.*c1Hzx - diff(Ey, 1, 2).*c2Hzx;
  Hzy = Hzy.*c1Hzy - diff(Ex, 1, 1).*c2Hzy;
  Hz = Hzx + Hzy;
  Ex(2 : end - 1, :) = Ex(2 : end - 1, :).*c1Ex - diff(Hz, 1, 1).*c2Ex;
  Ey(:, 2 : end - 1) = Ey(:, 2 : end - 1).*c1Ey - diff(Hz, 1, 2).*c2Ey;
  Ey(portr, portc) = Eyinc(t(i)) ...
    + (-1)*Ey(portr, portc - 1) - Eyinc(t(i) - dt - 1) ...
    + (dt - 1)/(dt + 1)*( ...
      Ey(portr, portc + 1) - Eyinc(t(i) - 1) ...
      -(-1)*Ey(portr, portc) + Eyinc(t(i) - dt) ...
    );
  Ey(portr, portc - 1) = Ey(portr, portc + 1);
  a1(i) = mean(Eyinc(t(i)))*sqrt(abs(porty2 - porty1));
  b1(i) = mean(Ey(portr, portc))*sqrt(abs(porty2 - porty1)) - a1(i);
  Ic(:, i) = interp2(Hz, column(circleX), row(circleY))*sqrt(2*pi*dphi*circleR);
  if (mod(i, downsampling) == 0)
    vH = (min(1, max(-1, Hz/constraint)) + 1)/2;
    vE = min(1, sqrt(meany(Ex).^2 + meanx(Ey).^2)/constraint);
    [rH, gH, bH] = ind2rgb(gray2ind(vH, colorN), jet(colorN)); 
    [rE, gE, bE] = ind2rgb(gray2ind(vE, colorN), hot(colorN));
    rE([pec; circle]) = 0.25;
    gE([pec; circle]) = 0.25;
    bE([pec; circle]) = 0.25;
    rH([pec; circle]) = 0.25;
    gH([pec; circle]) = 0.25;
    bH([pec; circle]) = 0.25; 
    imwrite(cat(3, [rE, rH], [gE, gH], [bE, bH]), ...
      sprintf("%s/t%04d.png", outDir, i));
  endif
endfor
toc(id)

Ft = exp(1i*2*pi*f'*t)'*sqrt(df*dt); % Fourier transform --- Time -> Frequency
modeN = 20;
Ff = [ % Fourier transform --- mode <- phi
  ones(size(phi)), ...
  sqrt(2)*cos(2*pi*phi*[1 : modeN]), ...
  sqrt(2)*sin(2*pi*phi*[1 : modeN]) ...
]'*sqrt(dphi);

near2far = (
  exp(1i*pi/2*([0 : modeN] + 1/2).')*sqrt(2./(pi*2*pi*f(if0 + 1 : end)*circleR))
)./besselh([0 : modeN], 2, 2*pi*f(if0 + 1 : end)'*circleR).';
near2far = [conj(fliplr(near2far)), zeros(modeN + 1, 1), near2far];
near2far = [near2far; near2far(2 : end, :)];

% Desired directivity characteristic
goal = sqrt(6)*(cos(2*pi*phi) > 0.5).*cos(3/2*2*pi*phi)*sqrt(dphi);

Ifar = real(Ff'*(Ff*Ic*Ft.*near2far)*Ft');
S11 = (b1*Ft)./(a1*Ft);
S21 = (goal'*Ff'*(Ff*Ic*Ft.*near2far))./(a1*Ft);
T = Ff'*(Ff*Ic*Ft.*near2far)./(a1*Ft);

g = mean(abs(interp1(f.', S21.', f1*[1 : 0.2 : 2].').').^2);
printf("\nPerformance: %.2f%%\n", g*100)

IfarImg = interp1(t.', Ifar.', t(downsampling) + [2*nc - 1 : -1 : 0].', 0).';
maxIfarImg = max(max(abs(IfarImg)));
IfarImg = ind2rgb(gray2ind((IfarImg/maxIfarImg + 1)/2, colorN), jet(colorN));
imwrite(IfarImg, sprintf("%s/far.png", outDir));

dB = @(x) 20*log10(abs(x));

figw = 1920/2; figh = 1080/2;
figure(1, 'position', [200, 200, figw, figh])

subplot(2, 3, 1)
plot(t*f1, a1/sqrt(dt))
xlabel("t*f1")
title("Stimulus --- Time Domain")
grid on
subplot(2, 3, 4)
plot(f/f1, dB(a1*Ft/max(abs(a1*Ft))))
xlim([0, 4])
ylim([-30, 0])
xlabel("f/f1")
ylabel("[dB]")
title("Stimulus --- Spectrum (normalized)")
grid on

subplot(2, 3, 2)
plot(t*f1, b1/sqrt(dt))
xlabel("t*f1")
title("Backward Response --- Time Domain")
grid on
subplot(2, 3, 5)
plot(f/f1, dB(S11), ";S11;", f/f1, dB(S21), ";S21;")
xlim([0, 4])
ylim([-30, 0])
xlabel("f/f1")
ylabel("[dB]")
legend("location", "southwest")
title("Scattering Parameters")
grid on

subplot(2, 3, 3)
phiv = [0; 60]/360;
plot(t*f1, interp1(phi, Ifar, phiv)/sqrt(dphi*dt))
xlabel("t*f1")
legend(num2str(phiv*360, "phi = %.0f deg"))
title("Far-Field Response --- Time Domain")
grid on
subplot(2, 3, 6)
plot(f/f1, dB(interp1(phi, T, phiv)/sqrt(dphi))) 
xlim([0, 4])
ylim([-20, 10])
xlabel("f/f1")
ylabel("[dBi]")
legend(num2str(phiv*360, "phi = %.0f deg"), "location", "southwest")
title("Realized Gain")
grid on

print(sprintf("%s/figure1.pdf", outDir), sprintf("-S%d,%d", figw, figh))

figure(2, 'position', [200, 200, figw, figh])
for q = [0:5]
  fq = (1 + q*0.2)*f1;
  Tq = exp(1i*2*pi*fq*deembedd)*interp1(f.', T.', fq.').';
  subplot(2, 3, q + 1)
  plot(phi*360, dB(Tq/sqrt(dphi)), ...
    phi*360, dB(goal/sqrt(dphi)), 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--')
  xlim([-180, 180])
  ylim([-20, 10])
  xticks([-180 : 60 : 180])
  yticks([-20 : 5 : 10])
  xlabel("phi [deg]")
  ylabel("|T| [dBi]")
  title(sprintf("Directivity Characteristic (f/f1 = %.1f)", fq/f1))
  grid on
endfor
print(sprintf("%s/figure2.pdf", outDir), sprintf("-S%d,%d", figw, figh))

figure(3, 'position', [200, 200, figw, figh])
for q = [0 : 5]
  fq = (1 + q*0.2)*f1;
  Tq = exp(1i*2*pi*fq*deembedd)*interp1(f.', T.', fq.').';
  subplot(2, 3, q + 1)
  plot(phi*360, arg(Tq)*180/pi)
  xlim([-180, 180])
  ylim([-180, 180])
  xticks([-180 : 60 : 180])
  yticks([-180 : 45 : 180])
  xlabel("phi [deg]")
  ylabel("arg(T) [deg]")
  title(sprintf("Directivity Characteristic (f/f1 = %.1f)", fq/f1))
  grid on
endfor
print(sprintf("%s/figure3.pdf", outDir), sprintf("-S%d,%d", figw, figh))