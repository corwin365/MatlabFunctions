function estSec = predict_fig2plotly_time(fig)
if nargin<1 || ~ishandle(fig), fig = gcf; end

% Count primitives
lines   = findall(fig,'Type','line');
nLine   = sum(arrayfun(@(h) numel(h.XData), lines));

sc      = findall(fig,'-isa','matlab.graphics.chart.primitive.Scatter');
nScatt  = sum(arrayfun(@(h) numel(h.XData), sc));

surfs   = findall(fig,'-isa','matlab.graphics.chart.primitive.Surface');
nSurf   = sum(arrayfun(@(h) numel(h.ZData), surfs));   % cells on grid

patches = findall(fig,'-isa','matlab.graphics.primitive.Patch');
nFaces  = sum(arrayfun(@(h) size(h.Faces,1), patches), 'omitnan');

imgs    = findall(fig,'Type','image');
nPix    = sum(arrayfun(@(h) numel(h.CData), imgs));

txt     = findall(fig,'Type','text');  nText = numel(txt);
axs     = findall(fig,'Type','axes');  nAxes = numel(axs);

% Rough weights (seconds per unit) — adjust to your machine if desired
w.line  = 1.0e-6;     % per line point
w.scatt = 1.2e-6;     % per scatter point
w.surf  = 3.0e-6;     % per surface grid cell
w.face  = 6.0e-6;     % per patch face
w.pix   = 2.0e-7;     % per image pixel (SVG images are heavy)
w.text  = 1.0e-3;     % per text object
w.axes  = 0.15;       % per axes (layout/annotations)
base    = 0.4;        % fixed overhead

estSec = base + ...
         w.line*nLine + w.scatt*nScatt + w.surf*nSurf + ...
         w.face*nFaces + w.pix*nPix + w.text*nText + w.axes*nAxes;
end