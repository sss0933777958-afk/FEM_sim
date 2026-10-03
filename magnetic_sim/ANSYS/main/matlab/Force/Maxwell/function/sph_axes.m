function out = sph_axes(model, geom, variant, Rum, L, opt)
%SPH_AXES  b.b and its gradient on the three actuator axes, from a solid-harmonic
%          model of given range and degree.
%
%   out = SPH_AXES(model, geom, variant, Rum, L)
%   out = SPH_AXES(model, geom, variant, Rum, L, opt)
%
%   The model is fitted on every node with |p| <= Rum at degree L, then
%   evaluated on NQ equidistant points over |s| <= Rum on each axis, with the
%   pole that axis points at driven at 1 A:  x_a <- P1, y_a <- P3, z_a <- P6
%   (the pairing of the smoothness gate).
%
%   opt.NQ    points per axis                         (default 2001)
%   opt.save  write data/<model>/.mat/sphax_R<R>_L<L>.mat   (default true)
%
%   Fields of the record, named as plot/<model>/sph_bb.m reads them:
%       s      NQ x 1   along-axis coordinate [um]
%       GS     NQ x 3   d(b.b)/ds on x_a, y_a, z_a      [mT^2/um], times SG
%       D3     NQ x 3   d3(b.b)/ds^3, signed as the gate signs it
%       bbS    NQ x 1   b.b on x_a                      [mT^2]
%       POLE, SPAN, LDEG, K, NMAE, nsc (sign changes of D3 per axis), spos
%
%   See also SPH_FIELD.

    if nargin < 6, opt = struct(); end
    NQ = 2001;   if isfield(opt,'NQ')   && ~isempty(opt.NQ),   NQ = opt.NQ;       end
    dosave = true; if isfield(opt,'save') && ~isempty(opt.save), dosave = opt.save; end

    here = fileparts(mfilename('fullpath'));          % .../Force/Maxwell/utils
    FMX  = fileparts(here);
    ROOT = fileparts(fileparts(FMX));                 % .../main/matlab
    CALR = fullfile(ROOT, 'Flux', 'Maxwell');
    addpath(fullfile(CALR,'function'), fullfile(CALR,'utils'), fullfile(CALR,'common_path'), ...
            fullfile(FMX,'function'), fullfile(FMX,'function'));

    cfg = model_config(model, geom);
    raw = extract_maxwell_data(cfg, 'all', variant);
    ad  = build_actuator_data(raw, cfg);
    inb = find(ad.r2 <= (Rum*1e-6)^2);
    [bq, gq, si] = sph_field(Rum, struct('P',ad.Pa(inb,:)*1e6,'B',ad.Ba(inb,:,:)), ...
                             struct('L',L,'Rn',Rum,'quiet',true,'condmax',1500));

    POLE = [1 3 6];   SG = [1 1 -1];
    s  = linspace(-Rum, Rum, NQ).';
    GS = zeros(NQ,3);   D3 = GS;   nsc = zeros(1,3);   spos = cell(1,3);
    for a = 1:3
        I  = zeros(6,1);   I(POLE(a)) = 1;
        P  = zeros(NQ,3);  P(:,a) = s;
        % signed like the gate (SG), as x150.mat stores it: z_a is walked towards
        % P6, which sits at -z_a, so the raw derivative there has the other sign
        g  = gq(P, I, 'du');    GS(:,a) = SG(a)*g(:,a);
        d3 = gq(P, I, 'd3');    D3(:,a) = SG(a)*d3(:,a);
        v  = D3(:,a);   k = find(diff(sign(v(v ~= 0))) ~= 0);
        sv = s(v ~= 0);   nsc(a) = numel(k);   spos{a} = sv(k).';
    end
    I  = zeros(6,1);   I(1) = 1;
    b  = bq([s zeros(NQ,2)], I);
    out = struct('model',model, 'variant',variant, 's',s, 'GS',GS, 'D3',D3, ...
                 'bbS',sum(b.^2,2), 'POLE',POLE, 'SG',SG, 'SPAN',Rum, 'R',Rum, ...
                 'LDEG',L, 'K',si.K, 'NMAE',si.NMAE_all, 'Nnodes',numel(inb), ...
                 'nsc',nsc, 'spos',{spos}, 'made',datestr(now,'yyyy-mm-dd HH:MM')); %#ok<TNOW1,DATST>
    fprintf('[sph_axes] R <= %g um, L = %d (K = %d), %d nodes, field NMAE %.4f %%%s', ...
            Rum, L, si.K, numel(inb), si.NMAE_all, newline);
    nm = {'x_a','y_a','z_a'};
    for a = 1:3
        fprintf('  %s (P%d): gradient %.4g .. %.4g mT^2/um, d3 sign changes %d at s = %s um%s', ...
                nm{a}, POLE(a), min(GS(:,a)), max(GS(:,a)), nsc(a), mat2str(round(spos{a})), newline);
    end
    if dosave
        fo = fullfile(FMX,'utils','data', sprintf('sphax_R%d_L%d.mat', round(Rum), L));
        save(fo, '-struct', 'out');
        fprintf('saved %s%s', fo, newline);
    end
end
