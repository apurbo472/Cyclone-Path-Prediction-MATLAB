%% Cyclone/Tornado Path Forecast with Environmental Inputs
% Methods used: interpolation (spline), extrapolation, root-finding, and
% a lightweight env-feature regression to correct motion.

clear; clc; close all;

%% 1) Synthetic "past" data (replace later with real data)
Tpast = (0:1:48)';                       % hours
N = numel(Tpast);

% Base track (slow NW curve). Units: degrees
lat_true = 18 + 0.05*Tpast - 0.0008*Tpast.^2;
lon_true = 90 - 0.08*Tpast + 0.0006*Tpast.^2;

% Intensity proxy (e.g., max wind speed in m/s)
I_true = 20 + 0.6*Tpast - 0.007*Tpast.^2;

% Environmental vars (toy but realistic-ish trends)
Temp = 28 - 0.02*Tpast + 0.4*sin(2*pi*Tpast/24);   % °C
Hum  = 75 + 5*sin(2*pi*(Tpast+6)/24);              % %
Pres = 1002 - 0.25*Tpast + 2*sin(2*pi*(Tpast+3)/24); % hPa

% Observational noise (simulate reality)
rng(7);
lat_obs = lat_true + 0.05*randn(N,1);
lon_obs = lon_true + 0.05*randn(N,1);
I_obs   = I_true   + 0.8*randn(N,1);

%% 2) Interpolation (splines) for a smooth path over past window
lat_spline = spline(Tpast, lat_obs);
lon_spline = spline(Tpast, lon_obs);
I_spline   = spline(Tpast, I_obs);

% Smooth reconstructions
tFinePast = linspace(Tpast(1), Tpast(end), 400)';
lat_fit = ppval(lat_spline, tFinePast);
lon_fit = ppval(lon_spline, tFinePast);
I_fit   = ppval(I_spline,   tFinePast);

%% 3) Motion model: derive velocity from spline derivatives
% d(lat)/dt and d(lon)/dt from spline piecewise polynomial
latd_pp = fnder(lat_spline,1);
lond_pp = fnder(lon_spline,1);

% Sample at integer hours (aligned with env data)
v_lat = ppval(latd_pp, Tpast);
v_lon = ppval(lond_pp, Tpast);

% Design matrix from environment features (centered)
X = [Temp Hum Pres];
muX = mean(X,1); sigX = std(X,[],1);  % standardize for stability
Z = (X - muX)./sigX;
Z = [ones(N,1) Z]; % add bias

% Fit simple linear models: v_lat ~ Z, v_lon ~ Z
beta_lat = Z \ v_lat;
beta_lon = Z \ v_lon;

%% 4) Extrapolate future environment and path (next 12 hours)
Hf = 12;
Tfuture = (Tpast(end)+1 : Tpast(end)+Hf)';

% Extrapolate env via last-24h linear trend + diurnal term (toy approach)
k = max(1, N-24+1):N;
tRef = Tpast(k);
trend = @(y) polyfit(tRef, y(k), 1);

pT = trend(Temp); pH = trend(Hum); pP = trend(Pres);

TempF = polyval(pT, Tfuture) + 0.4*sin(2*pi*Tfuture/24);
HumF  = polyval(pH, Tfuture) + 5*sin(2*pi*(Tfuture+6)/24);
PresF = polyval(pP, Tfuture) + 2*sin(2*pi*(Tfuture+3)/24);

XF = [TempF HumF PresF];
ZF = ([TempF HumF PresF] - muX)./sigX;
ZF = [ones(numel(Tfuture),1) ZF];

% Predict future velocities from env
v_lat_F = ZF * beta_lat;
v_lon_F = ZF * beta_lon;

% Initialize future lat/lon from last observed (using one-step from spline)
lat0 = ppval(lat_spline, Tpast(end));
lon0 = ppval(lon_spline, Tpast(end));

% Integrate forward with Euler (1 h step)
lat_pred = zeros(Hf,1); lon_pred = zeros(Hf,1);
lat_pred(1) = lat0 + v_lat_F(1)*1;
lon_pred(1) = lon0 + v_lon_F(1)*1;
for i = 2:Hf
    lat_pred(i) = lat_pred(i-1) + v_lat_F(i)*1;
    lon_pred(i) = lon_pred(i-1) + v_lon_F(i)*1;
end

% Extrapolate intensity via spline continuation (safer: linear tail)
pI = polyfit(Tpast(end-12+1:end), I_obs(end-12+1:end), 1);
I_pred = polyval(pI, Tfuture);

%% 5) Root-finding examples
% (a) When does it cross a target latitude (e.g., 23.5°N)?
targetLat = 23.5;

% Build a continuous past+future lat(t) via a piecewise function
% We'll stitch a small interpolant over [Tpast(end), Tfuture(end)]
tAll = [tFinePast; Tfuture];
latAll = [lat_fit; lat_pred];

% Define function f(t) = lat(t) - targetLat using local interpolation
latAll_spline = spline(tAll, latAll);
f_lat = @(t) ppval(latAll_spline, t) - targetLat;

% Bisection on a bracket around the last time window if sign change exists
tA = Tpast(end)-6; tB = Tfuture(end);
if f_lat(tA)*f_lat(tB) > 0
    tA = Tpast(1); % fallback wider search
end

% Simple bisection
a = tA; b = tB;
for iter = 1:60
    m = 0.5*(a+b);
    if f_lat(a)*f_lat(m) <= 0
        b = m;
    else
        a = m;
    end
end
t_cross_lat = 0.5*(a+b);

% (b) When will intensity hit a threshold (e.g., I = 30 m/s)?
IAll = [I_fit; I_pred];
IAll_spline = spline(tAll, IAll);
gI = @(t) ppval(IAll_spline, t) - 30;

% Newton method from last known time
t0 = Tpast(end);
for kNewt = 1:20
    % numerical derivative
    h = 1e-3;
    gt = gI(t0);
    dgt = (gI(t0+h)-gI(t0-h))/(2*h);
    t0 = t0 - gt/(dgt + 1e-9);
end
t_I30 = t0;

%% 6) 3D visualization (time vs lat vs lon), colored by intensity
t_show = tAll;
lat_show = latAll;
lon_show = latAll*0; % placeholder to preallocate

% But we do have lonAll:
lonAll = [lon_fit; lon_pred];
I_show = IAll;

figure('Color','w'); hold on; grid on;
% Past track
plot3(tFinePast, lat_fit, lon_fit, 'o-', 'LineWidth', 1.5, 'MarkerSize', 3);
% Future track
plot3(Tfuture, lat_pred, lon_pred, '*-', 'LineWidth', 1.5, 'MarkerSize', 6);

% Color by intensity along combined path
% Build a colored scatter for clarity
scatter3(tAll, [lat_fit; lat_pred], [lon_fit; lon_pred], 28, IAll, 'filled');
cb = colorbar; cb.Label.String = 'Intensity (arb., e.g., m/s)';

xlabel('Time (h)');
ylabel('Latitude (deg)');
zlabel('Longitude (deg)');
title('Cyclone Path Forecast (Interpolation + Env-Corrected Motion)');
view(45,25);
legend({'Past (spline)','Future (env-corrected)','Intensity-colored'}, 'Location','best');

%% 7) Print key results
fprintf('Target latitude %.2f° crossed at t = %.2f h\n', targetLat, t_cross_lat);
fprintf('Intensity I=30 reached near t = %.2f h\n', t_I30);

%% 8) Notes for real data
% - Replace synthetic Tpast, lat_obs, lon_obs, I_obs with IMD/NOAA track.
% - Replace Temp/Hum/Pres with reanalysis or station data along the track.
% - Consider geodesic corrections (deg->km varies with latitude).
% - Swap linear regression with ridge/GPR if you want smoother corrections.
% - Use a Kalman Filter to fuse (lat,lon) and env-driven motion for robustness.
