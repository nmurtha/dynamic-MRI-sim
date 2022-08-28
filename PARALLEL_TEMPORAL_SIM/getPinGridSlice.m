function [ pinGrids ] = getPinGridSlice( nr,nc )
%getPinGridSlice Generate pin grid spatial features.
%
%   Generates a 2D slice that can be put into a 3D phantom array. Contains
%   high resolution pin grid features, to test spatial resolution
%   recovery in frequency and phase encode directions.
%
%   This feature was originally designed for a phantom with slices
%   measuring 192x192, but this code now has the capability for general nr
%   x nc slice size. However, scaling capability is limited. The minimum is
%   175x175, beyond this the features will get bunched up. There's really
%   no maximum, but for very large sizes there will be lots of empty space
%   because the pin widths are fixed.
%
%    INPUTS
%     nr          -> Number of rows in the output image.
%     nc          -> Number of cols in the output image.
%
%    OUTPUTS
%     pinGrids    -> pin pair features. Measures (nr x nc). Output is a
%                    mask of 0's and 1's.
%
%    EXTERNALS
%     [none]

% Warn if array is too small
if ((nr < 175) || (nc < 175))
    disp('GETPINGRIDSLICE: Dimensions less than 175 will produce bunched features.')
end

% Allocate memory for pin pair mask.
pinGrids = zeros(nr,nc); 

% Set value to which pin pair mask features will be set to. 
val = 1;

% Set margin width for air outside phantom.
margin_r = 0.05 .* nr;
margin_c = 0.05 .* nc;

% Set the number of pins per side in the grid features.
npins = 4;

% Set number of features in the top, middle and bottom "rows"
ntop = 4;
nmid = 3;
nbot = 2;

% Set initial pin seperation and thickness, as well as the increase in
% seperation thickness that will be seen. These values will be increased as
% the features are built. Additionally, calculate the "unit width" of the
% smallest feature (i.e. sum of all pin widths and gaps between them). The
% larger features will be integer multiples of this unit width.
pin_sep = 1;
pin_thick = 1;
sep_incr = 1; 
w_unit = npins*pin_thick + (npins-1)*pin_sep;

%  ------------------------------------------------------------------------
%  ---------------------   TOP ROW OF FEATURES   --------------------------
%  ------------------------------------------------------------------------

% Calculate step length between features. The top features are
% designed such that they're all equally spaced by their own pin length.
steplength_top = round((nc-2*margin_c-max(cumsum(1:ntop))*w_unit)/(ntop+1));

% Get initial column to start building top features from.
col = round(margin_c + steplength_top);

% Construct the top set of pin grid pin features, resetting the row we
% begin construction at so that the top of the features are level.
for c = 1:ntop
    % Reset row to ensure features are level
    row = round(margin_r + w_unit);
    for t = 1:npins
        % Fill pins in current row
        for j = 0:2:2*npins-2
            pinGrids(row:row+pin_thick-1,col+j*pin_thick:col+(j+1)*pin_thick-1) = val;
        end
        
        % Increment row by appropriate gap for the next pin
        row = row + pin_thick + pin_sep;
    end
    % Step the column across the feature width and the step length to make
    % the next feature.
    col = col + c*w_unit + steplength_top;
    % Advance pin seperation and pin thickness for next feature.
    pin_sep = pin_sep + sep_incr;
    pin_thick = pin_thick + sep_incr;
end

% Get bottom row of last feature by "undoing" the row incrementation in the
% above loop.
row_mid = row - (pin_thick-sep_incr);



%  ------------------------------------------------------------------------
%  ---------------------   MID ROW OF FEATURES   --------------------------
%  ------------------------------------------------------------------------

% Calculate step length between features. The top features are
% designed such that they're all equally spaced by their own pin length.
steplength_mid = round((nc-2*margin_c-max(cumsum(ntop+1:nmid+ntop))*w_unit)/(nmid+1));

% Get vertical seperation between top and mid features. This seperation is
% calculated such that the seperation between the middle and bottom rows is
% equivalent.
delta = ((nr-2*margin_r)-w_unit*(3*ntop+2*nmid+nbot+2))/2;

% Get initial column to start building top features from.
col = round(margin_c + steplength_mid);

% Construct the mid set of pin grid pin features, resetting the row we
% begin construction at so that the top of the features are level.
for c = 1:nmid
    % Reset row to ensure features are level
    row = round(row_mid + delta);
    for t = 1:npins
        % Fill pins in current row
        for j = 0:2:2*npins-2
            pinGrids(row:row+pin_thick-1,col+j*pin_thick:col+(j+1)*pin_thick-1) = val;
        end
        
        % Increment row by appropriate gap for the next pin
        row = row + pin_thick + pin_sep;
    end
    % Step the column across the feature width and the step length to make
    % the next feature.
    col = col + (ntop+c)*w_unit + steplength_mid;
    % Advance pin seperation and pin thickness for next feature.
    pin_sep = pin_sep + sep_incr;
    pin_thick = pin_thick + sep_incr;
end

% Get bottom row of last feature by "undoing" the row incrementation in the
% above loop.
row_bot = row - (pin_thick-sep_incr);



%  ------------------------------------------------------------------------
%  ---------------------   BOT ROW OF FEATURES   --------------------------
%  ------------------------------------------------------------------------

% Calculate step length between features. The top features are
% designed such that they're all equally spaced by their own pin length.
steplength_bot = round((nc-2*margin_c-max(cumsum(ntop+nmid+1:nbot+nmid+ntop))*w_unit)/(nbot+1));

% Get initial column to start building top features from.
col = round(margin_c + steplength_bot);

% Construct the bot set of pin grid pin features, resetting the row we
% begin construction at so that the top of the features are level.
for c = 1:nbot
    % Reset row to ensure features are level
    row = round(row_bot + delta);
    for t = 1:npins
        % Fill pins in current row
        for j = 0:2:2*npins-2
            pinGrids(row:row+pin_thick-1,col+j*pin_thick:col+(j+1)*pin_thick-1) = val;
        end
        
        % Increment row by appropriate gap for the next pin
        row = row + pin_thick + pin_sep;
    end
    % Step the column across the feature width and the step length to make
    % the next feature.
    col = col + (ntop+nmid+c)*w_unit + steplength_bot;
    % Advance pin seperation and pin thickness for next feature.
    pin_sep = pin_sep + sep_incr;
    pin_thick = pin_thick + sep_incr;
end


end

