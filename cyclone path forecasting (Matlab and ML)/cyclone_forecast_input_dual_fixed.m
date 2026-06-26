%% cyclone_predict_fullpath_map_integrated.m
% Rolling predictions on fixed map region [-10,30] lat and [60,130] lon
% Land = green, Ocean = blue. Uses full observed path and predicts Hf_future ahead each step.
clear; clc; close all;

%% -------------------- User settings --------------------
env_file  = 'cyclone_data.csv';
path_file = 'cyclone_sidr_path_with_time.csv';
Hf_future = 3;          % predicted points at each step (change as needed)
minPointsToFit = 10;     % minimum observed points to fit the delta model
pause_sec = 0.12;       % pause between steps for animation
predict_mode = 'delta'; % only delta mode implemented here (recommended)
latLim = [10 25];      % [south north]
lonLim = [85 95];      % [west east]
landColor  = [0.4 0.8 0.4];   % green
oceanColor = [0.6 0.85 1];    % blue
% -------------------------------------------------------

%% ---- 0) load files ----
if ~isfile(env_file), error('Env file not found: %s', env_file); end
if ~isfile(path_file), error('Path file not found: %s', path_file); end
T_env  = readtable(env_file);
T_path = readtable(path_file);


% Remove rows where lat or lon is missing
latCols = {'lat','latitude','Lat','Latitude'};
lonCols = {'lon','longitude','Lon','Longitude'};

tmpLat = intersect(latCols, T_path.Properties.VariableNames);
tmpLon = intersect(lonCols, T_path.Properties.VariableNames);

if isempty(tmpLat) || isempty(tmpLon)
    error('Path CSV must contain lat/lon columns.');
end

% Drop rows with missing coordinates
T_path = rmmissing(T_path, 'DataVariables', {tmpLat{1}, tmpLon{1}});
fprintf('After cleaning: path rows=%d\n', height(T_path));





fprintf('Loaded: env rows=%d, path rows=%d\n', height(T_env), height(T_path));

%% ---- 1) datetime detection ----
dt_candidates = {'datetime','date','Date','time','Time'};
dtEnvCols = intersect(dt_candidates, T_env.Properties.VariableNames);
dtPathCols = intersect(dt_candidates, T_path.Properties.VariableNames);
hasEnvTime  = ~isempty(dtEnvCols);
hasPathTime = ~isempty(dtPathCols);

if hasEnvTime
    try envTimes = datetime(T_env.(dtEnvCols{1})); catch; envTimes = []; hasEnvTime=false; warning('Failed to parse env datetime.'); end
else envTimes = []; end
if hasPathTime
    try pathTimes = datetime(T_path.(dtPathCols{1})); catch; pathTimes = []; hasPathTime=false; warning('Failed to parse path datetime.'); end
else pathTimes = []; end

%% ---- 2) map each path row to nearest env row (so we iterate full path) ----
nPath = height(T_path);
nEnv  = height(T_env);
mapEnvIdx = zeros(nPath,1);
if hasEnvTime && hasPathTime
    for i = 1:nPath
        [~, mapEnvIdx(i)] = min(abs(envTimes - pathTimes(i)));
    end
elseif hasEnvTime && ~hasPathTime
    for i = 1:nPath, mapEnvIdx(i) = min(i, nEnv); end
else
    % no env times or no path times -> map by index (clamped)
    for i = 1:nPath, mapEnvIdx(i) = min(i, nEnv); end
end

%% ---- 3) extract env features robustly ----
firstCol = @(cands,names) intersect(cands,names);
% Temperature
tmp = firstCol({'temp','temperature','Temp','TempMax','TempMin'}, T_env.Properties.VariableNames);
if ~isempty(tmp), Temp_env = double(T_env.(tmp{1})); else Temp_env = nan(nEnv,1); end
% Humidity
tmp = firstCol({'humidity','Hum','Humidity'}, T_env.Properties.VariableNames);
if ~isempty(tmp), Hum_env = double(T_env.(tmp{1})); else Hum_env = nan(nEnv,1); end
% Pressure
tmp = firstCol({'sealevelpressure','pressure','Pressure'}, T_env.Properties.VariableNames);
if ~isempty(tmp), Pres_env = double(T_env.(tmp{1})); else Pres_env = nan(nEnv,1); end
% Wind speed
tmp = firstCol({'windspeed','wind_speed','WindSpeed','windspd'}, T_env.Properties.VariableNames);
if ~isempty(tmp), WindSpd_env = double(T_env.(tmp{1})); else WindSpd_env = nan(nEnv,1); end
% Wind dir
tmp = firstCol({'winddir','wind_dir','windDirection','windbearing'}, T_env.Properties.VariableNames);
if ~isempty(tmp), WindDir_env = double(T_env.(tmp{1})); else WindDir_env = nan(nEnv,1); end

Temp_env = Temp_env(:); Hum_env = Hum_env(:); Pres_env = Pres_env(:);
WindSpd_env = WindSpd_env(:); WindDir_env = WindDir_env(:);

%% ---- 4) extract full path coords ----
tmpLat = firstCol({'lat','latitude','Lat','Latitude'}, T_path.Properties.VariableNames);
tmpLon = firstCol({'lon','longitude','Lon','Longitude'}, T_path.Properties.VariableNames);
if isempty(tmpLat) || isempty(tmpLon), error('Path CSV must contain lat/lon columns.'); end
Lat_all = double(T_path.(tmpLat{1})); Lon_all = double(T_path.(tmpLon{1}));
Lat_all = Lat_all(:); Lon_all = Lon_all(:);



% % --- Fix mismatched NaNs in Lat/Lon for geoshow ---
% latNaN = isnan(Lat_all);
% lonNaN = isnan(Lon_all);
% 
% % If one is NaN and the other isn’t, make both NaN
% mismatch = xor(latNaN, lonNaN);
% if any(mismatch)
%     warning('Fixing %d mismatched NaN(s) in Lat/Lon', sum(mismatch));
%     Lat_all(mismatch) = NaN;
%     Lon_all(mismatch) = NaN;
% end
















figure('Name','Cyclone Sidr Observed Path','Color','w','Visible','on');
latLim = [-10 30];  % same as before
lonLim = [60 130];

% If Mapping Toolbox is available
if exist('worldmap','file') == 2
    worldmap(latLim, lonLim);
    load coastlines;
    geoshow(coastlat, coastlon, 'DisplayType','polygon', 'FaceColor',[0.4 0.8 0.4]); % land = green
    setm(gca,'FFaceColor',[0.6 0.85 1]);  % ocean = blue
    hold on;
    geoshow(Lat_all, Lon_all, 'DisplayType','line','Color','k','LineWidth',2);
    scatterm(Lat_all, Lon_all, 50, 'r', 'filled'); % points
else
    % fallback if Mapping Toolbox is not available
    plot(Lon_all, Lat_all, '-k','LineWidth',2); hold on;
    scatter(Lon_all, Lat_all, 50,'r','filled');
    set(gca,'Color',[0.6 0.85 1]); % ocean
    axis([lonLim latLim]); axis equal;
    xlabel('Longitude'); ylabel('Latitude');
    % draw land manually (rough approximation)
    patch([-180 180 180 -180],[-90 -90 90 90],[0.4 0.8 0.4],'EdgeColor','none');
end

title('Cyclone Sidr Observed Path');
legend('Observed Path','Observed Points','Location','best');
grid on;












% Check for NaNs quickly
if all(isnan(Lat_all)) || all(isnan(Lon_all)), error('All lat or lon are NaN — cannot plot.'); end
fprintf('First 3 observed coords: \n'); disp([Lat_all(1:min(3,end)) Lon_all(1:min(3,end))]);

%% ---- 5) Prepare figure & map (Mapping toolbox if available) ----
hasMapping = (exist('worldmap','file') == 2) && (exist('geoshow','file') == 2);
figure('Name','Cyclone rolling prediction on fixed region','Color','w','Visible','on','Units','normalized','Position',[0.06 0.06 0.88 0.86]);

if hasMapping
    try
        worldmap(latLim, lonLim);
        % draw land polygon (coastlines dataset has coastlat/coastlon)
        try
            load coastlines; % coastlat, coastlon
            geoshow(coastlat, coastlon, 'DisplayType','polygon', 'FaceColor', landColor, 'EdgeColor','none');
        catch
            % if coastlines missing, fill map background as ocean
        end
        setm(gca,'FFaceColor', oceanColor);
        hold on;
    catch ME
        warning('Mapping toolbox worldmap failed: %s. Falling back to plain axes.', ME.message);
        hasMapping = false;
    end
end

if ~hasMapping
    ax = axes('Parent', gcf); hold(ax,'on'); box(ax,'on');
    % attempt to load alternate coastline variable if present
    try
        load coast; % older files
        patch(long, lat, landColor, 'EdgeColor','none'); % if long/lat exist
    catch
        % set ocean background
        set(gca,'Color', oceanColor);
    end
    % axes limits
    xlim(lonLim); ylim(latLim);
    axis equal;
    xlabel('Longitude (°E)'); ylabel('Latitude (°N)');
end

% plot full observed path faint (reference)
if hasMapping
    hObsAll = geoshow(Lat_all, Lon_all, 'DisplayType','line', 'LineStyle', ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
else
    hObsAll = plot(Lon_all, Lat_all, ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
end
hold on;

% placeholders for dynamic handles
hObsLine = []; hObsPts = []; hPredLine = []; hPredPts = []; hInfo = [];

% place info text
if hasMapping
    hInfo = textm(latLim(2)-0.8, lonLim(1)+1, '', 'FontSize',12, 'FontWeight', 'bold');
else
    hInfo = text(lonLim(1)+1, latLim(2)-0.8, '', 'FontSize',12, 'FontWeight', 'bold');
end

drawnow;

%% ---- 6) Rolling loop through full observed path ----
N = numel(Lat_all);
for i = 1:N
    % Observed up to current i
    obsLat = Lat_all(1:i);
    obsLon = Lon_all(1:i);
    obsT   = (1:i)'; % use simple index-based time for fitting purposes
    
    % Map env rows for those observations
    envIdx_for_obs = mapEnvIdx(1:i);
    % Build X_cols: time + available env features (in same order)
    X_cols = obsT; fNames = {'Time'};
    if ~all(isnan(Temp_env(envIdx_for_obs))), X_cols = [X_cols, Temp_env(envIdx_for_obs)]; fNames{end+1}='Temp'; end
    if ~all(isnan(Hum_env(envIdx_for_obs))),  X_cols = [X_cols, Hum_env(envIdx_for_obs)];  fNames{end+1}='Hum'; end
    if ~all(isnan(Pres_env(envIdx_for_obs))), X_cols = [X_cols, Pres_env(envIdx_for_obs)]; fNames{end+1}='Pres'; end
    if ~all(isnan(WindSpd_env(envIdx_for_obs))), X_cols = [X_cols, WindSpd_env(envIdx_for_obs)]; fNames{end+1}='WindSpd'; end
    if ~all(isnan(WindDir_env(envIdx_for_obs))), X_cols = [X_cols, WindDir_env(envIdx_for_obs)]; fNames{end+1}='WindDir'; end
    
    % Attempt to fit delta model if enough data
    doPredict = false; Lat_pred = []; Lon_pred = [];
    if i >= minPointsToFit
        X_train = X_cols(1:end-1,:); % features at time t to predict delta t->t+1
        dLat = diff(obsLat); dLon = diff(obsLon);
        if size(X_train,1) == numel(dLat) && numel(dLat) >= (minPointsToFit-1)
            muX = mean(X_train,1); sigX = std(X_train,[],1); sigX(sigX==0)=1;
            Ztrain = (X_train - muX) ./ sigX;
            Ztrain = [ones(size(Ztrain,1),1), Ztrain];
            beta_dlat = (Ztrain' * Ztrain) \ (Ztrain' * dLat);
            beta_dlon = (Ztrain' * Ztrain) \ (Ztrain' * dLon);
            doPredict = true;
        end
    end
    
    % If we can predict, extrapolate features for Hf_future steps and predict dLat/dLon
    if doPredict
        Tfuture = (obsT(end)+1 : obsT(end)+Hf_future)'; % simple index future steps
        m = size(X_cols,2);
        X_future = zeros(Hf_future, m);




for j = 1:m
    col_j = X_cols(:,j);
    if numel(unique(obsT))>1 && sum(~isnan(col_j))>=2
        % Replace straight line with spline interpolation
X_future(:,j) = interp1(obsT, col_j, Tfuture, 'linear', 'extrap');
    else
        X_future(:,j) = repmat(col_j(end), Hf_future,1);
    end
end







        Zf = (X_future - muX) ./ sigX; Zf = [ones(size(Zf,1),1), Zf];
       
        
        dLat_pred = Zf * beta_dlat; dLon_pred = Zf * beta_dlon;
        
        
        %Lat_pred = obsLat(end) + cumsum(dLat_pred);
        %Lon_pred = obsLon(end) + cumsum(dLon_pred);
   





    % obsLat, obsLon are the observed coordinates collected so far
Lat_pred = [obsLat(end); obsLat(end) + cumsum(dLat_pred)];
Lon_pred = [obsLon(end); obsLon(end) + cumsum(dLon_pred)];

    
    
    
    
    end
    
    %% ---- 7) Update plot: remove previous dynamic handles and redraw ----
    try
        if ~isempty(hObsLine), delete(hObsLine); end
        if ~isempty(hObsPts),  delete(hObsPts);  end
        if ~isempty(hPredLine), delete(hPredLine); end
        if ~isempty(hPredPts),  delete(hPredPts);  end
    catch
        % ignore deletion errors
    end
    










% observed-so-far (smooth black)
numSmooth = 200; % number of points for smooth curve
validIdx = find(~isnan(obsLat) & ~isnan(obsLon));
if numel(validIdx) >= 2
    t  = 1:numel(validIdx);
    tt = linspace(1, numel(validIdx), numSmooth);
    smoothLat = spline(t, obsLat(validIdx), tt);
    smoothLon = spline(t, obsLon(validIdx), tt);
else
    smoothLat = obsLat;
    smoothLon = obsLon;
end

if hasMapping
    hObsLine = plotm(smoothLat, smoothLon, '-k', 'LineWidth', 2.2);
    hObsPts  = scatterm(obsLat, obsLon, 40, 'k', 'filled');
else
    hObsLine = plot(smoothLon, smoothLat, '-k', 'LineWidth', 2.2); hold on;
    hObsPts  = scatter(obsLon, obsLat, 40, 'k', 'filled');
    axis equal;
    xlim(lonLim); ylim(latLim);
end



    

    
    % predicted future (red)
    if doPredict && ~isempty(Lat_pred)
        if hasMapping
            hPredLine = plotm(Lat_pred, Lon_pred, '-.', 'Color', [1 0 0], 'LineWidth', 2);
            hPredPts  = scatterm(Lat_pred, Lon_pred, 60, 'r', 'filled');
        else
            hPredLine = plot(Lon_pred, Lat_pred, '-.', 'Color', [1 0 0], 'LineWidth', 2);
            hPredPts  = scatter(Lon_pred, Lat_pred, 60, 'r', 'filled');
        end
    end
    
    % update info text
    infoStr = sprintf('Step %d / %d', i, N);
    if doPredict
        infoStr = [infoStr sprintf('  (predicted %d ahead)', Hf_future)];
    else
        infoStr = [infoStr '  (insufficient data to predict)'];
    end
    if hasMapping
        set(hInfo,'String',infoStr);
    else
        set(hInfo,'String',infoStr,'Position',[lonLim(1)+0.5, latLim(2)-0.8]);
    end
    
    drawnow;
    pause(pause_sec);
end

fprintf('Done rolling predictions across entire observed path.\n');