function hAxes = dscatter(X, varargin)
% DSCATTER creates a scatter plot coloured by density.
%
%   DSCATTER(X,Y) creates a scatterplot of X and Y at the locations
%   specified by the vectors X and Y (which must be the same size), colored
%   by the density of the points.
%
%   DSCATTER(...,'MARKER',M) allows you to set the marker for the
%   scatter plot. Default is 's', square.
%
%   DSCATTER(...,'MSIZE',MS) allows you to set the marker size for the
%   scatter plot. Default is 10.
%
%   DSCATTER(...,'FILLED',false) sets the markers in the scatter plot to be
%   outline. The default is to use filled markers.
%
%   DSCATTER(...,'PLOTTYPE',TYPE) allows you to create other ways of
%   plotting the scatter data. Options are "surf','mesh' and 'contour'.
%   These create surf, mesh and contour plots colored by density of the
%   scatter data.
%
%   DSCATTER(...,'BINS',[NX,NY]) allows you to set the number of bins used
%   for the 2D histogram used to estimate the density. The default is to
%   use the number of unique values in X and Y up to a maximum of 200.
%
%   DSCATTER(...,'SMOOTHING',LAMBDA) allows you to set the smoothing factor
%   used by the density estimator. The default value is 20 which roughly
%   means that the smoothing is over 20 bins around a given point.
%
%   DSCATTER(...,'LOGY',true) uses a log scale for the yaxis.
%
%   Examples:
%
%       [data, params] = fcsread('SampleFACS');
%       dscatter(data(:,1),10.^(data(:,2)/256),'log',1)
%       % Add contours
%       hold on
%       dscatter(data(:,1),10.^(data(:,2)/256),'log',1,'plottype','contour')
%       hold off
%       xlabel(params(1).LongName); ylabel(params(2).LongName);
%       
%   See also FCSREAD, SCATTER.

% Copyright 2003-2004 The MathWorks, Inc.
% $Revision:  $   $Date:  $

% Reference:
% Paul H. C. Eilers and Jelle J. Goeman
% Enhancing scatterplots with smoothed densities
% Bioinformatics, Mar 2004; 20: 623 - 628.

% Modified by Ethan Blackwood in 2020 to accept a third data argument, Z,
% to create a 3D scatterplot.
% First input can also be a n x p matrix, where p = 2 or 3, instead of
% inputting each dimension separately.

lambda = [];
nbins = [];
plottype = 'scatter';
msize = 10;
marker = 's';
logy = false;
filled = true;
is3d = false;
ndims = 2;

if isvector(X)
    X = X(:);
    
    % combine Y and possibly Z vectors into one matrix
    assert(nargin >= 2 && ~ischar(varargin{1}), 'Not enough input variables (less than 2)!');
    assert(numel(varargin{1}) == size(X, 1), 'Input size mismatch');
    X(:, 2) = varargin{1}(:);
    
    varargin = varargin(2:end);
    
    if nargin >= 3 && ~ischar(varargin{1})
        is3d = true;
        ndims = 3;
        assert(numel(varargin{1}) == size(X, 1), 'Input size mismatch');
        X(:, 3) = varargin{1}(:);
        varargin = varargin(2:end);
        
        assert(ischar(varargin{1}), 'Too many input variables (more than 3)!');
    end
else
    ndims = size(X, 2);
    assert(ndims == 2 || ndims == 3, 'Too many input variables (more than 3)!');
    is3d = ndims == 3;
end

% remove rows with nans
X = X(all(~isnan(X), 2), :);

n = size(X, 1);
    
nkwargs = length(varargin);

if nkwargs > 0
    if rem(nkwargs, 2) == 1
        error('Bioinfo:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'smoothing','bins','plottype','logy','marker','msize','filled'};
    for j=1:2:nkwargs
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs); %#ok
        if isempty(k)
            error('Bioinfo:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % smoothing factor
                    if isnumeric(pval)
                        lambda = pval;
                    else
                        error('Bioinfo:InvalidScoringMatrix','Invalid smoothing parameter.');
                    end
                case 2
                    if isscalar(pval)
                        nbins = repmat(pval, 1, ndims);
                    else
                        nbins = pval;
                    end
                case 3
                    plottype = pval;
                    if is3d
                        error('Cannot specify plot type for 3D input');
                    end
                case 4
                    logy = pval;
                    X(:, 2) = log10(X(:, 2));
                case 5
                    marker = pval;
                case 6
                    msize = pval;
                case 7
                    filled = pval;
            end
        end
    end
end

mins = min(X);
maxes = max(X);

if isempty(nbins)
    nbins = arrayfun(@(dim) min(numel(unique(X(:, dim))), 200), 1:ndims);
end

if isempty(lambda)
    lambda = 20;
end

edges = cell(ndims, 1);
ctrs = cell(ndims, 1);
bin = zeros(n,ndims);
% Reverse the columns to put the first column of X along the horizontal
% axis, the second along the vertical.
dimorder = 1:ndims;
dimorder(1) = 2; dimorder(2) = 1;

for kD = 1:ndims
    edges{kD} = linspace(mins(kD), maxes(kD), nbins(kD)+1);
    ctrs{kD} = edges{kD}(1:end-1) + .5*diff(edges{kD});
    edges{kD} = [-Inf edges{kD}(2:end-1) Inf];
    
    bin(:,dimorder(kD)) = discretize(X(:, kD), edges{kD});
end

density = accumarray(bin,1,nbins(dimorder)) ./ n;

for kD = 1:ndims
    density = shiftdim(smooth1D(density, nbins(kD) / lambda), 1);
end
% = filter2D(H,lambda);

if logy
    ctrs{2} = 10.^ctrs{2};
    X(:, 2) = 10.^X(:, 2);
end

if is3d
    density = density ./ max(density(:));
    ind = sub2ind(size(density), bin(:, 1), bin(:, 2), bin(:, 3));
    col = density(ind);
    
    if filled
        h = scatter3(X(:, 1), X(:, 2), X(:, 3), msize, col, marker, 'filled');
    else
        h = scatter3(X(:, 1), X(:, 2), X(:, 3), msize, col, marker);
    end
    
else
    okTypes = {'surf','mesh','contour','image','scatter'};
    k = strmatch(lower(plottype), okTypes); %#ok
    if isempty(k)
        error('dscatter:UnknownPlotType',...
            'Unknown plot type: %s.',plottype);
    elseif length(k)>1
        error('dscatter:AmbiguousPlotType',...
            'Ambiguous plot type: %s.',plottype);
    else
        switch(k)
            
            case 1 %'surf'
                h = surf(ctrs{1},ctrs{2},density,'edgealpha',0);
            case 2 % 'mesh'
                h = mesh(ctrs{1},ctrs{2},density);
            case 3 %'contour'
                [~, h] =contour(ctrs{1},ctrs{2},density);
            case 4 %'image'
                nc = 256;
                density = density./max(density(:));
                colormap(repmat(linspace(1,0,nc)',1,3));
                h =image(ctrs{1},ctrs{2},floor(nc.*density) + 1);
            case 5 %'scatter'
                density = density./max(density(:));
                ind = sub2ind(size(density),bin(:,1),bin(:,2));
                col = density(ind);
                if filled
                    h = scatter(X(:, 1), X(:, 2),msize,col,marker,'filled');
                else
                    h = scatter(X(:, 1), X(:, 2),msize,col,marker);
                end
        end        
    end
end

if logy
    set(gca,'yscale','log');
end
if nargout > 0
    hAxes = get(h,'parent');
end
%%%% This method is quicker for symmetric data.
% function Z = filter2D(Y,bw)
% z = -1:(1/bw):1;
% k = .75 * (1 - z.^2);
% k = k ./ sum(k);
% Z = filter2(k'*k,Y);
end


function Z = smooth1D(Y,lambda)
m = size(Y, 1);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1; 
Z = reshape((E + P) \ reshape(Y, m, []), size(Y));
end


