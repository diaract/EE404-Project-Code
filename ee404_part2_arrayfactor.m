function ee404_part2_arrayfactor
%% ========================================================================
%  EE 404 Antennas - Project Part 2
%  Array factor and total radiation pattern of an N-element uniform linear
%  array, evaluated for a USER-SUPPLIED element pattern (passed as a
%  function handle, exactly as the assignment states: "for a given antenna
%  radiation pattern ...").
%
%  Geometry: array along x, half-wave dipoles along z, spacing d = lambda/4.
%  The angle 'ang' is measured FROM THE ARRAY AXIS (x):
%       ang = 90 deg  -> broadside ;  ang = 0/180 deg -> end-fire.
%  In the CST "Constant Theta, Cut Angle = 90" cut this equals CST's phi.
%
%  Structure:
%     arrayFactor()  - closed-form |AF| with the psi->0 limit handled
%     dipoleFigure8() - half-wave dipole element factor (x-z cut)
%     beamMetrics()  - main-beam angle, HPBW and side-lobe level (auto)
%     rectPlot()/polarPlot() - plotting helpers (polar floor handled by RLim)
% =========================================================================

    cfg.N       = 5;
    cfg.f       = 1.3e9;
    cfg.c       = 299792458;
    cfg.lambda  = cfg.c / cfg.f;
    cfg.d       = cfg.lambda / 4;          % 57.69 mm
    cfg.phases  = [0 15 30 45];            % progressive phase shifts [deg]
    cfg.floorDB = -40;

    angDeg = linspace(0, 360, 4001);
    angRad = deg2rad(angDeg);

    % element models: {title suffix, x-axis label, element function handle}
    % NOTE: the x-y (CST theta=90) cut sweeps CST's azimuth phi; the x-z cut
    % sweeps a DIFFERENT in-plane angle (alpha) from the array axis.
    afXlab = '\phi from array axis (x-y / CST \theta=90 cut) (deg)';
    models = { 'Omni (x-y / CST \theta=90)', afXlab,                                    @(a) ones(size(a)) ; ...
               'Dipole element (x-z cut)',   '\alpha from array axis (x-z cut) (deg)',  @dipoleFigure8 };
    nCol = size(models,1) + 1;

    fprintf('Phase |  AF main-beam | HPBW(AF) | SLL(AF)\n');
    fprintf('------+---------------+----------+--------\n');
    for p = cfg.phases
        AF = arrayFactor(cfg.N, cfg.lambda, cfg.d, angRad, deg2rad(p));
        m  = beamMetrics(angDeg, AF);
        fprintf('%4d  |   %6.1f deg  | %5.1f deg | %5.1f dB\n', ...
                 p, m.peakAngle, m.hpbw, m.sll_dB);

        figure('Name',sprintf('Phase shift = %d deg',p),'Color','w', ...
               'Position',[60 60 1200 700]);

        % column 1 : array factor (x-y / CST cut)
        rectPlot (subplot(2,nCol,1),      angDeg, AF, 'Array Factor', afXlab, cfg.floorDB);
        polarPlot(subplot(2,nCol,nCol+1), angRad, AF, 'Array Factor', cfg.floorDB);

        % remaining columns : total pattern = element x AF
        for e = 1:size(models,1)
            EF = models{e,3}(angRad);  EF = EF / max(EF);
            T  = EF .* AF;
            rectPlot (subplot(2,nCol,1+e),        angDeg, T, ...
                      ['Total: ' models{e,1}], models{e,2}, cfg.floorDB);
            polarPlot(subplot(2,nCol,nCol+1+e),   angRad, T, ...
                      ['Total: ' models{e,1}], cfg.floorDB);
        end
    end

    % combined beam-steering overlay (array factor, CST cut)
    figure('Name','AF beam steering (all phases)','Color','w'); hold on; grid on;
    for p = cfg.phases
        AF = arrayFactor(cfg.N, cfg.lambda, cfg.d, angRad, deg2rad(p));
        plot(angDeg, max(20*log10(AF+eps), cfg.floorDB), 'LineWidth',1.6, ...
             'DisplayName', sprintf('\\Phi = %d\\circ', p));
    end
    xlim([0 360]); ylim([cfg.floorDB 0]);
    xlabel('\phi from array axis (deg)'); ylabel('Normalized AF [dB]');
    title('Array-factor beam steering (x-y / CST \theta=90 cut)');
    legend('Location','best');

    fprintf(['\nAbsolute level: read the gain from CST (~7.26 dBi at 0 deg).\n' ...
             '2.15 + 10log10(N) is only a rough rule and overestimates at d = lambda/4.\n']);
end

%% ------------------------------------------------------------------------
function AF = arrayFactor(N, lambda, d, angRad, betaRad)
% Closed-form magnitude of a uniform linear array factor, normalized.
    k   = 2*pi / lambda;
    psi = k*d*cos(angRad) + betaRad;        % inter-element phase progression
    s   = sin(psi/2);
    AF  = sin(N*psi/2) ./ s;                % |AF| = |sin(N psi/2)/sin(psi/2)|
    AF(abs(s) < 1e-9) = N;                   % limit value where psi -> 2*pi*m
    AF  = abs(AF);
    AF  = AF / max(AF);
end

%% ------------------------------------------------------------------------
function EF = dipoleFigure8(angRad)
% Half-wave dipole along z, observed in the x-z plane; angle from array axis.
    EF = abs( cos((pi/2).*sin(angRad)) ./ cos(angRad) );
    EF(abs(cos(angRad)) < 1e-6) = 0;         % enforce null along dipole axis
    EF(~isfinite(EF)) = 0;                    % (cos(pi/2) is ~6e-17, not 0)
end

%% ------------------------------------------------------------------------
function m = beamMetrics(angDeg, patLin)
% Main-beam angle, 3 dB beamwidth and worst side-lobe level (all approximate).
    patDB   = 20*log10(patLin/max(patLin) + eps);
    [~, i0] = max(patLin);
    m.peakAngle = angDeg(i0);
    m.hpbw      = halfPowerWidth(angDeg, patDB, i0);
    m.sll_dB    = sideLobeLevel(angDeg, patDB);
end

function w = halfPowerWidth(angDeg, patDB, i0)
    n = numel(patDB);  thr = -3;  iR = i0;  iL = i0;
    for s = 1:n-1
        j = i0 + s;  if j > n, j = j - n; end
        if patDB(j) <= thr, iR = j; break; end
    end
    for s = 1:n-1
        j = i0 - s;  if j < 1, j = j + n; end
        if patDB(j) <= thr, iL = j; break; end
    end
    w = angDeg(iR) - angDeg(iL);
    if w < 0, w = w + 360; end
end

function sll = sideLobeLevel(angDeg, patDB)
% Highest side-lobe level over the unique 0..180 deg half (avoids counting
% the symmetric mirror main beam). Main lobe = region around the peak down
% to its first minima; the peak of everything outside it is the SLL.
    idx = find(angDeg <= 180);
    p   = patDB(idx);
    [~, ip] = max(p);
    iL = ip;  while iL > 1        && p(iL-1) < p(iL),  iL = iL - 1;  end
    iR = ip;  while iR < numel(p) && p(iR+1) < p(iR),  iR = iR + 1;  end
    mask = true(size(p));  mask(iL:iR) = false;     % drop the main lobe
    if any(mask), sll = max(p(mask)); else, sll = NaN; end
end

%% ------------------------------------------------------------------------
function rectPlot(ax, angDeg, patLin, ttl, xlab, floorDB)
    dB = max( 20*log10(patLin/max(patLin) + eps), floorDB );
    plot(ax, angDeg, dB, 'LineWidth', 1.6);  grid(ax, 'on');
    xlim(ax, [0 360]);  ylim(ax, [floorDB 0]);
    xlabel(ax, xlab);  ylabel(ax, 'Normalized [dB]');
    title(ax, ttl);
end

function polarPlot(ax, angRad, patLin, ttl, floorDB)
    dB  = max( 20*log10(patLin/max(patLin) + eps), floorDB );
    pos = ax.Position;  delete(ax);            % swap cartesian tile -> polar
    pax = polaraxes('Position', pos);
    polarplot(pax, angRad, dB, 'LineWidth', 1.6);
    pax.ThetaZeroLocation = 'top';  pax.ThetaDir = 'clockwise';
    pax.RLim = [floorDB 0];
    title(pax, ttl);
end
