%% Assignment #1: Quasi-static HEV model and Rule Based Control

%% Group information
% Group number: 47
% Students:
%   Student 1, sXXXXXX
%   Student 2, sXXXXXX
%   Student 3, sXXXXXX

%% Load the cycle and vehicle data
% Clear the workspace and add all required folders to the MATLAB path.
clear
addpath("models")
addpath("data")
addpath("utilities")

% Load the WLTP driving cycle data.
mis    = load(fullfile("data", "WLTP.mat"));
time   = mis.time_s;               % [s]
vehSpd = mis.speed_km_h ./ 3.6;   % convert to [m/s]
vehAcc = mis.acceleration_m_s2;   % [m/s^2]

% Load and scale vehicle data for Group 47:
%   Engine power    = 118 000 W
%   E-Machine power =  27 000 W
%   Battery energy  =   1 900 Wh
veh = load(fullfile("data", "vehData.mat"));
veh = scaleVehData(veh, 118e3, 27e3, 1.9e3);

% Load gear-shift schedule for the transmission controller.
load(fullfile("data", "transmControlData.mat"))

%% Simulation loop
% Initialize SOC and starting gear.
x(1) = 0.6;   % initial SOC = 60%
GN0  = 1;     % initial gear number

for n = 1:length(time)
    % Transmission controller (provided): selects gear from vehicle speed.
    GN(n) = gearControl(vehSpd(n), GN0, upSpd, downSpd);

    % Powerflow controller (our RBC): computes engine torque-split factor.
    engAlpha(n) = powerflowControl(x(n), GN(n), vehSpd(n), vehAcc(n), veh);

    % Advance the HEV model one timestep.
    [x(n+1), stageCost, unfeas, engPrf(n), emPrf(n), battPrf(n), vehPrf(n)] = ...
        hev_model(x(n), [GN(n), engAlpha(n)], [vehSpd(n), vehAcc(n)], veh);

    GN0 = GN(n);
end

%% Results analysis
% Convert non-scalar struct arrays into scalar structs (Nx1 arrays).
engPrf  = structArray2struct(engPrf);
emPrf   = structArray2struct(emPrf);
battPrf = structArray2struct(battPrf);
vehPrf  = structArray2struct(vehPrf);

% Pack all profiles into a single parent structure.
prof.engPrf  = engPrf;
prof.emPrf   = emPrf;
prof.battPrf = battPrf;
prof.vehPrf  = vehPrf;

% Compute fuel consumption, fuel economy and final SOC.
fuelConsumption = trapz(time, engPrf.fuelFlwRate) .* 1e-3;              % [kg]
dist            = trapz(time, vehSpd) * 1e-3;                           % [km]
fuelEconomy     = fuelConsumption / (veh.eng.fuelDensity * dist) * 1e2; % [l/100km]
finalSOC        = battPrf.battSOC(end);

% Print results.
fprintf("Fuel consumption: %.2f kg\n",      fuelConsumption)
fprintf("Fuel economy:     %.2f l/100km\n", fuelEconomy)
fprintf("Final SOC:        %.3f\n",         finalSOC)

% Plot time profiles: speed, SOC, gear, alpha, cumulative fuel.
mainProfiles(prof);

% Plot component power profiles.
powerProfiles(prof);

% Plot engine operating points on the BSFC map.
engMapWithPF(veh.eng, prof, 'bsfc');

% Plot e-machine operating points on the efficiency map.
emMapWithPF(veh.em, prof);

%% Save results
% Transform and pack profiles (already done above).
% Store results
save("results.mat", "prof", "fuelConsumption", "fuelEconomy", "finalSOC")

%% (Optional) Improving the torque-split controller
% Re-run the simulation using powerflowControlExtra, which adds hard SOC
% thresholds (soc_min = 0.40, soc_max = 0.80) on top of the base logic.

x_ex(1) = 0.6;
GN0_ex  = 1;

for n = 1:length(time)
    GN_ex(n)       = gearControl(vehSpd(n), GN0_ex, upSpd, downSpd);
    engAlpha_ex(n) = powerflowControlExtra(x_ex(n), GN_ex(n), vehSpd(n), vehAcc(n), veh);
    [x_ex(n+1), ~, ~, engPrf_ex(n), emPrf_ex(n), battPrf_ex(n), vehPrf_ex(n)] = ...
        hev_model(x_ex(n), [GN_ex(n), engAlpha_ex(n)], [vehSpd(n), vehAcc(n)], veh);
    GN0_ex = GN_ex(n);
end

engPrf_ex  = structArray2struct(engPrf_ex);
emPrf_ex   = structArray2struct(emPrf_ex);
battPrf_ex = structArray2struct(battPrf_ex);
vehPrf_ex  = structArray2struct(vehPrf_ex);

prof.engPrf  = engPrf_ex;
prof.emPrf   = emPrf_ex;
prof.battPrf = battPrf_ex;
prof.vehPrf  = vehPrf_ex;

fuelConsumption = trapz(time, engPrf_ex.fuelFlwRate) .* 1e-3;
fuelEconomy     = fuelConsumption / (veh.eng.fuelDensity * dist) * 1e2;
finalSOC        = battPrf_ex.battSOC(end);

fprintf("--- Extra Feature Results ---\n")
fprintf("Fuel consumption: %.2f kg\n",      fuelConsumption)
fprintf("Fuel economy:     %.2f l/100km\n", fuelEconomy)
fprintf("Final SOC:        %.3f\n",         finalSOC)

mainProfiles(prof);

% Store results
save("results_extra.mat", "prof", "fuelConsumption", "fuelEconomy", "finalSOC")

%% The torque-split controller
% -----------------------------------------------------------------------
% FUNCTION: powerflowControl
%
% Rule-based torque-split controller for a p2 parallel HEV.
% The engine is kept as close as possible to its Optimal Operating Line
% (OOL). The EM covers the remaining torque or charges the battery.
%
% INPUTS:
%   soc     - battery state of charge [-]  (needed: charge/discharge decision)
%   gearNum - current gear [1-5]           (needed: shaft speed and T_dem)
%   spd     - vehicle speed [m/s]          (needed: kinematics, EV threshold)
%   acc     - vehicle acceleration [m/s^2] (needed: traction force)
%   veh     - vehicle data struct          (needed: all maps and parameters)
%
% OUTPUT:
%   alpha   - engine torque-split factor = T_eng / T_dem [-]
%             0 = full electric  |  1 = full thermal
%
% KEY PARAMETERS:
%   soc_tgt = 0.60  SOC charge-sustaining target
%   v_pe    = 5.0   speed threshold for pure-electric mode [m/s]
% -----------------------------------------------------------------------
function alpha = powerflowControl(soc, gearNum, spd, acc, veh)

    % Tunable parameters
    soc_tgt = 0.60;
    v_pe    = 5.0;

    % Driveline kinematics
    omega_wh  = spd / veh.wh.r;
    F_veh     = veh.body.f0 + veh.body.f1*spd + veh.body.f2*spd^2 + veh.body.m*acc;
    T_wh      = F_veh * veh.wh.r;
    tau_gb    = veh.gb.tau(gearNum);
    omega_sh  = omega_wh * veh.fd.tau * tau_gb;
    T_dem     = T_wh   / (veh.fd.tau * tau_gb);

    % Engine OOL torque and max torque at this shaft speed
    T_ool    = veh.eng.ool(omega_sh);
    T_engMax = veh.eng.maxTrq(omega_sh);

    % EM torque limits converted to equivalent torque at gearbox input
    omega_em  = omega_sh * veh.tc.tau;
    Tgb_emMax =  veh.tc.tau * veh.em.maxTrq(omega_em);
    Tgb_emMin =  veh.tc.tau * veh.em.minTrq(omega_em);

    % Decision tree
    if T_dem <= 0
        % Braking / coasting: engine off, EM regenerates
        alpha = 0;

    elseif spd < v_pe
        % Low speed: pure electric
        alpha = 0;

    elseif T_dem > T_ool && soc >= soc_tgt
        % High demand, battery OK: power-split, engine near OOL
        T_eng = max(T_ool, T_dem - Tgb_emMax);
        T_eng = min(T_eng, T_engMax);
        T_eng = max(T_eng, 0);
        alpha = T_eng / T_dem;

    elseif soc < soc_tgt
        % SOC low: battery charging — engine above OOL, charges via EM
        T_eng = min(T_ool, T_dem - Tgb_emMin);
        T_eng = min(T_eng, T_engMax);
        T_eng = max(T_eng, 0);
        alpha = T_eng / T_dem;

    else
        % Moderate demand, battery OK: engine on OOL, EM covers rest
        T_eng = T_ool;
        T_eng = min(T_eng, T_engMax);
        T_eng = max(T_eng, 0);
        alpha = T_eng / T_dem;
    end

    % Clamp to physical bounds
    alpha = max(0, min(1, alpha));
end

%% (Optional) A torque-split controller, with an extra feature
% -----------------------------------------------------------------------
% FUNCTION: powerflowControlExtra
%
% Extended RBC with hard SOC threshold enforcement (Extra Feature #1).
% Same inputs and output as powerflowControl.
%
% EXTRA PARAMETERS:
%   soc_min = 0.40  lower bound: forces engine-dominant mode to stop discharge
%   soc_max = 0.80  upper bound: forces maximum EM use to stop overcharge
% -----------------------------------------------------------------------
function alpha = powerflowControlExtra(soc, gearNum, spd, acc, veh)

    % Tunable parameters
    soc_tgt = 0.60;
    soc_min = 0.40;
    soc_max = 0.80;
    v_pe    = 5.0;

    % Driveline kinematics (same as base controller)
    omega_wh  = spd / veh.wh.r;
    F_veh     = veh.body.f0 + veh.body.f1*spd + veh.body.f2*spd^2 + veh.body.m*acc;
    T_wh      = F_veh * veh.wh.r;
    tau_gb    = veh.gb.tau(gearNum);
    omega_sh  = omega_wh * veh.fd.tau * tau_gb;
    T_dem     = T_wh   / (veh.fd.tau * tau_gb);

    T_ool    = veh.eng.ool(omega_sh);
    T_engMax = veh.eng.maxTrq(omega_sh);
    omega_em  = omega_sh * veh.tc.tau;
    T_emMax   = veh.em.maxTrq(omega_em);
    T_emMin   = veh.em.minTrq(omega_em);
    Tgb_emMax =  veh.tc.tau * T_emMax;
    Tgb_emMin =  veh.tc.tau * T_emMin;

    % Extended decision tree
    if T_dem <= 0
        alpha = 0;

    elseif soc < soc_min
        % Battery critically low: engine dominant, EM at minimum
        T_em  = max(T_emMin, 0);         % do not discharge further
        T_eng = T_dem - veh.tc.tau * T_em;
        T_eng = min(T_eng, T_engMax);
        T_eng = max(T_eng, 0);
        alpha = T_eng / T_dem;

    elseif soc > soc_max
        % Battery too full: maximise EM discharge (Max EM mode)
        if (T_dem / veh.tc.tau) <= T_emMax
            alpha = 0;                   % pure electric suffices
        else
            T_eng = T_dem - Tgb_emMax;
            T_eng = min(T_eng, T_engMax);
            T_eng = max(T_eng, 0);
            alpha = T_eng / T_dem;
        end

    elseif spd < v_pe
        alpha = 0;

    elseif T_dem > T_ool && soc >= soc_tgt
        T_eng = max(T_ool, T_dem - Tgb_emMax);
        T_eng = min(T_eng, T_engMax);
        T_eng = max(T_eng, 0);
        alpha = T_eng / T_dem;

    elseif soc < soc_tgt
        T_eng = min(T_ool, T_dem - Tgb_emMin);
        T_eng = min(T_eng, T_engMax);
        T_eng = max(T_eng, 0);
        alpha = T_eng / T_dem;

    else
        T_eng = T_ool;
        T_eng = min(T_eng, T_engMax);
        T_eng = max(T_eng, 0);
        alpha = T_eng / T_dem;
    end

    % Clamp to physical bounds
    alpha = max(0, min(1, alpha));
end
