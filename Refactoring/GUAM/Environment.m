classdef Environment < handle
    properties
        WindSpd_kts     = 0;    % WindSpd (kts)          Wind speed     [knots]
        WindDir_deg     = 0;    % WindDir (True, deg)    Wind direction [deg]
        VerticalWind    = 0;    % VerticalWind           수직 바람   [ft/s]

        TurbulenceLevel = 0;    % Turbulence Level       난류 강도 레벨 (0=None,1=Light,2=Moderate,3=Severe)
        TurbSwitch      = 1;    % Level Select / Switch  난류 on/off

        dT              = 0;    % atmos76 dt 입력        표준대기 온도편차 [deg R]

        Ts              = 0.01; % 고정 스텝 샘플 시간     [s] (난류 이산화에 사용)

        turbulence_level = [...
            0,   0,      0;...      % turb_Ku, turb_Kv, turb_Kw  (Level 0: None)
            0.2, 0.2,    0.17;...   % turb_Ku, turb_Kv, turb_Kw  (Level 1: Light)
            0.4, 0.4,    0.34;...   % turb_Ku, turb_Kv, turb_Kw  (Level 2: Moderate)
            0.4, 1.152,  0.34];     % turb_Ku, turb_Kv, turb_Kw  (Level 3: Severe)

        RandomSeeds     = [23341, 23342, 23343, 23344]; % [u_turb, v_turb, w_turb, p_gust] 시드

        dt
    end

    properties (Constant)
        KTS2FPS = 1.6878099;    % -kts2fps   knots  -> ft/s
        D2R     = pi/180;       % d2r        deg    -> rad
        FPS2KTS = 0.5924838;    % fpstoknots ft/s   -> knots
        M2FT    = 3.2808399;    % m2ft 3     m      -> ft
    end

    %% ------------------------------------------------------------------
    %  내부 상태 (Dryden 필터의 이산 상태 등, 적분기 블록 대응)
    %  ------------------------------------------------------------------
    properties (Access = private)
        turbState   = [];       % Dryden 난류 필터 이산 상태 (struct: xu,xv,xw,xp,xq,xr)
        noiseStream = [];       % 난류 백색잡음용 RandStream (struct: u,v,w,p) - 재현성
    end

    %% ==================================================================
    methods
        % --------------------------------------------------------------
        function obj = Environment(params)
            % if nargin >= 1 && ~isempty(params)
            %     f = fieldnames(params);
            %     for k = 1:numel(f)
            %         if isprop(obj, f{k})
            %             obj.(f{k}) = params.(f{k});
            %         end
            %     end
            % end
            % obj.reset();

            obj.dt = params.dt;
        end
    end     
    methods (Static)
        function band_limited_white_noise(obj)
            % Generate band-limited white noise (zero-mean) using an FIR lowpass filter.
            % Usage: y = obj.bandLimitedWhiteNoise(nSamples, fs, fc)
            % Inputs:
            %   nSamples - number of output samples
            %   fs       - sampling frequency (Hz)
            %   fc       - cutoff frequency (Hz) for lowpass (fc < fs/2)
            % Output:
            %   y        - column vector of length nSamples containing band-limited noise
            %
            % Notes:
            % - Uses a windowed-sinc FIR lowpass filter (Hann window).
            % - Normalizes output to unit standard deviation.

            fs = 1/obj.dt;
            fc = obj.fc;
            if nargin < 4
                error('bandLimitedWhiteNoise requires (obj, nSamples, fs, fc).');
            end
            if fc <= 0 || fc >= fs/2
                error('Cutoff fc must satisfy 0 < fc < fs/2.');
            end
            % Generate white Gaussian noise
            x = randn();
            % Design FIR lowpass: choose filter length proportional to fs/fc
            % but keep it reasonable and odd.
            % Transition width ~ fc/4 => N ~= 4*fs/fc
            N = round(4 * fs / fc);
            N = max(N, 11);            % minimum length
            if mod(N,2)==0, N = N+1; end
            M = (N-1)/2;
            % Ideal sinc lowpass (normalized frequency)
            fcn = fc / fs; % normalized cutoff (0..0.5)
            n = (-M:M).';
            h = 2 * fcn * sinc(2 * fcn * n); % sinc(x) = sin(pi x)/(pi x) in MATLAB's sinc
            % Apply Hann window
            w = (0.5 + 0.5*cos(2*pi*n/(N-1))).';
            h = h .* w;
            % Normalize filter to unity RMS gain for white noise -> normalize so output std = 1
            % For white noise input with unit variance, output variance = sum(h.^2)
            h = h / sqrt(sum(h.^2));
            % Filter with zero-phase capability using filtfilt if signal long enough
            if nSamples > 3*N
                y = filtfilt(h,1,x);
            else
                y = filter(h,1,x);
            end
            % Ensure zero mean and unit std
            y = y - mean(y);
            y = y / std(y);
        end
    end
end