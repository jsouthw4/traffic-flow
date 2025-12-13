%% Ring-road baseline (IDM) + jam emergence
clear; close all; clc;

%% Parameters
L     = 260;        % ring length (m)
ell   = 4.5;        % vehicle length (m)
dt    = 0.05;       % time step (s)
Tend  = 300;        % total time (s)
Nt    = round(Tend/dt) + 1;
t     = (0:Nt-1)*dt;

% IDM parameters (tweakable)
p.v0    = 30;       % desired speed (m/s)
p.a     = 1.2;      % max acceleration (m/s^2)
p.b     = 3.0;      % 2.0 is comfortable braking (m/s^2). Greater is strong breaking
p.T     = 1.2;      % 1.2 approx desired time headway (s)
p.s0    = 1.0;      % minimum spacing (m)
p.delta = 4;        % acceleration exponent

% Experiments (densities): choose car counts
N_list = [24, 28, 32, 60];

% For fundamental diagram measurement window (avoid transient)
t_meas_start = 150;     % seconds
i_meas0 = find(t >= t_meas_start, 1, 'first');

%% Storage for FD results
rho_list = zeros(size(N_list));
q_list   = zeros(size(N_list));
vbar_list= zeros(size(N_list));

%% Loop over densities
for run = 1:numel(N_list)
    N = N_list(run);

    % Initial positions equally spaced
    s = (0:N-1)' * (L/N);

    % Initial speeds near desired with small noise
    v = p.v0 * ones(N,1) + 0.5*randn(N,1);
    v(v<0) = 0;

    % Small perturbation to seed waves
    v(1) = 0.5*p.v0;

    % Logs for plots (space-time and speed heatmap)
    s_log = zeros(N,Nt);
    v_log = zeros(N,Nt);

    s_log(:,1) = s;
    v_log(:,1) = v;

    % For flow estimate: count laps (when car crosses s=0)
    lap_count = zeros(N,1);

    % Sim loop (semi-implicit Euler: update v then s)
    for k = 2:Nt
        % Leader indices (circular)
        ip1 = [2:N 1]';

        % Gaps and relative speeds
        ds_raw = s(ip1) - s;
        ds_mod = mod(ds_raw, L);          % wrap
        gap    = ds_mod - ell;
        gap(gap < 0.1) = 0.1;             % avoid blowups

        dv     = v - v(ip1);

        % IDM desired gap
        sstar = p.s0 + v*p.T + (v.*dv)./(2*sqrt(p.a*p.b));
        sstar(sstar < p.s0) = p.s0;

        % IDM acceleration
        acc = p.a*(1 - (v/p.v0).^p.delta - (sstar./gap).^2);

        % Integrate
        v = v + dt*acc;
        v(v<0) = 0;

        s_prev = s;
        s = s + dt*v;

        % Wrap positions and count laps
        crossed = (s_prev < L) & (s >= L);
        lap_count = lap_count + crossed;
        s = mod(s, L);

        % Log
        s_log(:,k) = s;
        v_log(:,k) = v;
    end

    %% Fundamental diagram metrics for this run
    rho = N / L;  % veh/m
    rho_list(run) = rho;

    % Average speed over measurement window
    vbar = mean(mean(v_log(:,i_meas0:end), 2));
    vbar_list(run) = vbar;

    % Flow estimate: total laps completed after measurement start
    laps_meas = lap_count;  % laps over whole sim, but we want only after i_meas0
    % Better: recompute laps only after i_meas0 from s_log
    laps2 = 0;
    for i = 1:N
        si = s_log(i,i_meas0:end);
        
        % Count wrap-around events (vehicle passes reference point once per lap)
        laps2 = laps2 + sum(diff(si) < -0.5*L); % wrap detected by big negative jump
    end
    time_meas = Tend - t_meas_start;
    q = (laps2 / time_meas);     % vehicles per second crossing s=0
    q_list(run) = q;

    %% Plots for a couple representative runs
figure('Color','w');
s_unwrap = unwrapPositions(s_log, L);
imagesc(t, 1:N, s_unwrap);
axis xy; xlabel('time (s)'); ylabel('car index');
title(sprintf('Space-time diagram (N=%d, density=%.4f veh/m)', N, rho));
colorbar;

figure('Color','w');
imagesc(t, 1:N, v_log);
axis xy; xlabel('time (s)'); ylabel('car index');
title(sprintf('Speed heatmap (N=%d)', N));
colorbar;

end

%% Fundamental diagram plot
figure('Color','w');
plot(rho_list, q_list, 'o-', 'LineWidth', 1.5);
xlabel('density \rho (veh/m)'); ylabel('flow q (veh/s)');
title('Fundamental diagram from ring simulation');

figure('Color','w');
plot(rho_list, vbar_list, 's-', 'LineWidth', 1.5);
xlabel('density \rho (veh/m)'); ylabel('mean speed (m/s)');
title('Mean speed vs density');

%% Helper: unwrap positions for space-time clarity
function s_unwrap = unwrapPositions(s_log, L)
    [N,Nt] = size(s_log);
    s_unwrap = zeros(N,Nt);
    for i = 1:N
        s_unwrap(i,1) = s_log(i,1);
        for k = 2:Nt
            ds = s_log(i,k) - s_log(i,k-1);
            if ds < -0.5*L
                ds = ds + L;
            elseif ds > 0.5*L
                ds = ds - L;
            end
            s_unwrap(i,k) = s_unwrap(i,k-1) + ds;
        end
    end
end
