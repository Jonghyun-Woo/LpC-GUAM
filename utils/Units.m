classdef Units
    % Unit conversion constants — instance mirrors setUnits(distUnit, massUnit).
    %
    % Wraps lib/utilities/setUnits.m as a class so unit systems are
    % selectable at construction time while keeping the same field names.
    %
    % Allowable distUnit : 'm', 'cm', 'ft', 'in', 'SI'  (SI → 'm'+'kg')
    % Allowable massUnit : 'kg', 'gram', 'slug', 'slinch', 'lbm'
    %   (metric and English units cannot be mixed)
    %
    % Usage:
    %   u = Units('ft', 'slug')   % English — default for GUAM simulation
    %   u = Units('m',  'kg')     % SI
    %   u = Units()               % defaults to 'ft', 'slug'
    %
    %   wingspan = 10 * u.m       % 10 m expressed in base-unit (ft)
    %   V_kts    = V_fps / u.knot % ft/s → knots

    properties
        % --- Base: distance ---
        ft;  in;  m;  cm

        % --- Base: mass ---
        slug;  slinch;  lbm;  kg;  gram

        % --- Base: time / angle ---
        s;  rad

        % --- Temperature ---
        K;  degR;  zeroC;  zeroF

        % --- Derived: distance ---
        km;  mm;  nmile;  smile

        % --- Derived: time ---
        hour

        % --- Derived: speed ---
        knot

        % --- Derived: acceleration ---
        g0

        % --- Derived: area ---
        m2;  ft2;  in2

        % --- Derived: force ---
        lbf;  N;  dyne;  poundal

        % --- Derived: pressure ---
        Pa;  mbar;  psf;  psi;  atmos;  mmHg;  inHg

        % --- Derived: energy ---
        J;  erg;  ftlbf

        % --- Derived: power ---
        W;  hp

        % --- Derived: angle ---
        deg;  arcmin;  arcsec
    end

    methods
        function obj = Units(distUnit, massUnit)
            % Mirrors setUnits(distUnit, massUnit) — see lib/utilities/setUnits.m.
            % Default: Units('ft', 'slug')

            % --- defining constants (NIST SP811) ---
            defineFoot   = 0.3048;       % foot = 0.3048 m
            definePound  = 0.45359237;   % lbm  = 0.45359237 kg
            defineG      = 9.80665;      % g0   = 9.80665 m/s²
            defineKRankine = 1.8;        % 1 K  = 1.8 °R
            zeroCelsius    = 273.15;
            zeroFahrenheit = zeroCelsius * defineKRankine - 32.0;
            defineNM     = 1852.0;       % nautical mile = 1852 m
            defineAtmos  = 101325.0;     % standard atmosphere [Pa]
            defineMmHg   = 133.3224;     % mm Hg [Pa]

            % --- defaults ---
            if nargin < 1,  distUnit = 'ft';    end
            if nargin < 2,  massUnit = 'slug';  end
            if strcmpi(distUnit, 'SI')
                distUnit = 'm';  massUnit = 'kg';
            end

            % --- distance base unit ---
            if strcmpi(distUnit, 'm')
                metric = true;  m = 1.0;       cm = 0.01;
            elseif strcmpi(distUnit, 'cm')
                metric = true;  m = 100.0;     cm = 1.0;
            elseif strcmpi(distUnit, 'ft')
                metric = false; ft = 1.0;      in = 1.0/12.0;
            elseif strcmpi(distUnit, 'in')
                metric = false; ft = 12.0;     in = 1.0;
            else
                error('Units: invalid distUnit "%s".', distUnit);
            end

            if metric
                ft = defineFoot * m;    in = ft / 12.0;
            else
                m  = ft / defineFoot;   cm = m / 100.0;
            end

            g0 = defineG * m;   % standard free fall in base distance unit

            % --- mass base unit ---
            if strcmpi(massUnit, 'kg')
                if ~metric, error('Units: cannot mix "kg" with "%s".', distUnit); end
                kg = 1.0;  gram = 0.001;
            elseif strcmpi(massUnit, 'gram')
                if ~metric, error('Units: cannot mix "gram" with "%s".', distUnit); end
                kg = 1000.0;  gram = 1.0;
            elseif strcmpi(massUnit, 'slug')
                if metric,  error('Units: cannot mix "slug" with "%s".', distUnit); end
                slug = 1.0;  slinch = 12.0;  lbm = slug * (ft/g0);
            elseif strcmpi(massUnit, 'slinch')
                if metric,  error('Units: cannot mix "slinch" with "%s".', distUnit); end
                slug = 1/12.0;  slinch = 1.0;  lbm = slug * (ft/g0);
            elseif strcmpi(massUnit, 'lbm')
                if metric,  error('Units: cannot mix "lbm" with "%s".', distUnit); end
                lbm = 1.0;  slug = lbm * (g0/ft);  slinch = 12.0 * slug;
            else
                error('Units: invalid massUnit "%s".', massUnit);
            end

            if metric
                lbm    = definePound * kg;
                slug   = lbm * (g0/ft);
                slinch = 12.0 * slug;
            else
                kg   = lbm / definePound;
                gram = 0.001 * kg;
            end

            % --- base units ---
            obj.m  = m;    obj.cm = cm;
            obj.ft = ft;   obj.in = in;

            obj.kg      = kg;    obj.gram   = gram;
            obj.lbm     = lbm;   obj.slug   = slug;
            obj.slinch  = slinch;

            obj.s   = 1.0;
            obj.rad = 1.0;

            if metric
                obj.K    = 1.0;                obj.degR = 1.0 / defineKRankine;
            else
                obj.K    = defineKRankine;     obj.degR = 1.0;
            end
            obj.zeroC = zeroCelsius;
            obj.zeroF = zeroFahrenheit;

            % --- derived: distance ---
            obj.km    = 1000.0    * obj.m;
            obj.mm    = 0.001     * obj.m;
            obj.nmile = defineNM  * obj.m;
            obj.smile = 5280.0    * obj.ft;

            % --- derived: time ---
            obj.hour = 3600.0 * obj.s;

            % --- derived: speed ---
            obj.knot = obj.nmile / obj.hour;

            % --- derived: acceleration ---
            obj.g0 = g0;

            % --- derived: area ---
            obj.m2  = obj.m  * obj.m;
            obj.ft2 = obj.ft * obj.ft;
            obj.in2 = obj.in * obj.in;

            % --- derived: force ---
            obj.lbf     = obj.lbm  * obj.g0;
            obj.N       = obj.kg   * obj.m;
            obj.dyne    = obj.gram  * obj.cm;
            obj.poundal = obj.lbm  * obj.ft;

            % --- derived: pressure ---
            obj.Pa    = obj.N   / obj.m2;
            obj.mbar  = 100.0   * obj.Pa;
            obj.psf   = obj.lbf / obj.ft2;
            obj.psi   = obj.lbf / obj.in2;
            obj.atmos = defineAtmos * obj.Pa;
            obj.mmHg  = defineMmHg  * obj.Pa;
            obj.inHg  = obj.mmHg * obj.in / obj.mm;

            % --- derived: energy ---
            obj.J     = obj.N    * m;          % m is the local variable (= base dist)
            obj.erg   = obj.dyne * obj.cm;
            obj.ftlbf = obj.lbf  * obj.ft;

            % --- derived: power ---
            obj.W  = obj.J / obj.s;
            obj.hp = 550.0 * obj.ftlbf / obj.s;

            % --- derived: angle ---
            obj.deg    = pi / 180.0;
            obj.arcmin = obj.deg    / 60.0;
            obj.arcsec = obj.arcmin / 60.0;
        end
    end
end
