

% cbrew.m NPH UoB 20170901

% Generate cbrewer colormap with simpler inputs. Also generate NPH custom
% colormaps.

% ONE INPUT:
% Specify either the name of the colormap or it's number below:
%
% cmap = cbrew(1);      or   cmap = cbrew('Spectral');
% cmap = cbrew(5);      or   cmap = cbrew('RdBu');
% cmap = cbrew(21);     or   cmap = cbrew('Oranges');
%
% default length is 64 element length colormaps.
% Name is case-insensitive.

% TWO INPUTS:
% Also specify the length of the colormap:
%
% cmap = cbrew(1,32);   or   cmap = cbrew('Spectral',32);
% cmap = cbrew(5,32);   or   cmap = cbrew('RdBu',32);
%
% Note that divergent colormaps are flipud()'ed as per the CBrewer Preview.

% EDIT: NOW ADDED NEW PIVOT CAPABILITY!!!!
% cbrew('RdBu',64,'pivot',0.75);


function cmap = cbrew(varargin)

% Same order as in the preview supplied with cbrewer:
cmapnames = {
    'Spectral'  , 'div'
    'RdYlGn'    , 'div'
    'RdYlBu'    , 'div'
    'RdGy'      , 'div'
    'RdBu'      , 'div'
    'PuOr'      , 'div'
    'PRGn'      , 'div'
    'PiYG'      , 'div'
    'BrBG'      , 'div'
    'YlOrRd'    , 'seq'
    'YlOrBr'    , 'seq'
    'YlGnBu'    , 'seq'
    'YlGn'      , 'seq'
    'Reds'      , 'seq'
    'RdPu'      , 'seq'
    'Purples'   , 'seq'
    'PuRd'      , 'seq'
    'PuBuGn'    , 'seq'
    'PuBu'      , 'seq'
    'OrRd'      , 'seq'
    'Oranges'   , 'seq'
    'Greys'     , 'seq'
    'Greens'    , 'seq'
    'GnBu'      , 'seq'
    'BuPu'      , 'seq'
    'BuGn'      , 'seq'
    'Blues'     , 'seq'
    'Set3'      , 'qual'
    'Set2'      , 'qual'
    'Set1'      , 'qual'
    'Pastel2'   , 'qual'
    'Pastel1'   , 'qual'
    'Paired'    , 'qual'
    'Dark2'     , 'qual'
    'Accent'    , 'qual'};


customcmaps = {...
    'nph_RdYlBu'
    'nph_Spectral'
    'nph_Parula'
    'nph_Rainbow'
    'nph_BuOr'
    'nph_BuOr2'
    'nph_BlueOrange'
    'nph_BuOrRd'
    'nph_alt_jet'
    'nph_RdYlBuGrey'
    'nph_modspectral'
    'nph_RdBuPastel'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLORMAP LENGTH
switch length(varargin)
    case 0 % NO INPUTS, default
        cmap = flipud(cbrewer('div','RdBu',64));
        return;
    case 1 % ONE INPUT
        len = 64; % default 64x3 element colormap vector
    otherwise % TWO OR MORE INPUTS
%         if isnumeric(varargin{2})
            len = varargin{2}; % user specified length       
%         else
%             len = 64;
%         end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLORMAP TYPE
switch isnumeric(varargin{1})
    
    case 1 % NUMERIC INPUT
        
        cname = cmapnames{varargin{1},1};
        ctype = cmapnames{varargin{1},2};
        
        cmap = cbrewer(ctype,cname,len);
        
    case 0 % STRING INPUT
        
        % parse input:
        switch any(strcmpi(varargin{1},cmapnames(:,1)))
            
            case 1 % if it is a recognised cbrewer one:
                
                % get correct letter case from list
                cname = cmapnames{strcmpi(varargin{1},cmapnames(:,1)),1};
                ctype = cmapnames{strcmpi(varargin{1},cmapnames(:,1)),2};
                
                cmap = cbrewer(ctype,cname,len);
            
            case 0
                % else try my custom colormaps:
                switch any(strcmpi(varargin{1},customcmaps))
                    case 1
                        cmap = custom_cmap(varargin{1},len);
                        ctype = 'custom';
                    otherwise % else not recognised, return
                        disp(['Unrecognised colormap: ' varargin{1}])
                    return
                end
                
            
        end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRECT FOR LENGTHS OF LESS THAN EIGHT (cbrewer doesn't like them, not full range)
switch ctype
    case {'qual' 'seq' 'div'}
        if len < 11
            fullcmap = cbrewer(ctype,cname,11);
            cmap = interp1(1:11,fullcmap,linspace(2,10,len));
        end
    otherwise
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLIP DIVERGING COLORMAPS:
% for some reason div colorbars are flipped compared to the cbrewer preview.
switch ctype
    case 'div'
        cmap = flipud(cmap);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIVOT
% Add the ability to change the pivot of the selected colormap. So for a
% diverging colormap, move the zero position.
% Pivot values are currently defined as between 0 and 1, where 0 means all
% the way to the bottom and 1 mean all the way to the top. I might change
% this at some point, but it's the simplest way without including actual
% clims into this function. By default pivot is 0.5, i.e. in the middle.

if any(strcmpi('pivot',varargin))
    
    pivot = varargin{find(strcmpi('pivot',varargin))+1};
    
    % Parse input
        if pivot < 0 || pivot > 1
            disp('Pivot value must be between 0 and 1.')
            pivot = 0.5;
        end
        
        
    if pivot ~= 0.5 % requires a change
        
        cmapi = griddedInterpolant({linspace(1,len,len),1:3},cmap,'linear');
        
        midpoint = floor(0.5*len) + 1;
        
        len1 = floor(len * pivot);
        len2 = len - len1;
        
%         if any(strcmpi('power',varargin))
%             pow = varargin{find(strcmpi('power',varargin))+1};
%             cmap = cmapi({[powerspace(1,midpoint,len1,pow) powerspace(midpoint,len,len2,pow)],1:3});
%         else
            cmap = cmapi({[linspace(1,midpoint,len1) linspace(midpoint,len,len2)],1:3});
%         end
        
        
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % POWER SPACE
% % tired of always having linearly progressing colormaps? Well now with
% % powerspace, you can have a colormap that ramps in any power you want!
% % Slow it down to be data = sqrt(data) or speed up to be data =
% % log10(data) or whatever.
% 
% % POW should be a 1 or 2 element vector for 'seq' and 'div' colormaps.
% 
% if any(strcmpi('power',varargin))
%     
%     % Extract power space value
%     pow = varargin{find(strcmpi('power',varargin))+1};
%     
%     % Pivot will have already had the correct power spacing applied above:
%     if any(strcmpi('pivot',varargin))
%         pow = 1;
%     end
%     
%     cmapi = griddedInterpolant({linspace(1,len,len),1:3},cmap,'linear');
%         
%     
% end








end % end function





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ccmapi = custom_cmap(ccmap_name,len)

switch lower(ccmap_name)
    case lower('nph_RdYlBu')
        ccmap = [...
            0.192156862745098         0.211764705882353         0.584313725490196
            0.270588235294118         0.458823529411765         0.705882352941177
            0.454901960784314         0.67843137254902          0.819607843137255
            0.670588235294118         0.850980392156863         0.913725490196078
            0.87843137254902          0.952941176470588         0.972549019607843
            1                         1                         1
            0.996078431372549         0.87843137254902          0.564705882352941
            0.992156862745098         0.682352941176471         0.380392156862745
            0.956862745098039         0.427450980392157         0.262745098039216
            0.843137254901961         0.188235294117647         0.152941176470588
            0.647058823529412                         0         0.149019607843137];
        
    case lower('nph_Spectral')
        ccmap = [...
            0.368627450980392         0.309803921568627         0.635294117647059
            0.196078431372549         0.533333333333333         0.741176470588235
                          0.4          0.76078431372549         0.647058823529412
            0.670588235294118         0.866666666666667         0.643137254901961
            0.901960784313726          0.96078431372549         0.596078431372549
                            1                         1                         1
            0.996078431372549          0.87843137254902         0.545098039215686
            0.992156862745098         0.682352941176471         0.380392156862745
            0.956862745098039         0.427450980392157         0.262745098039216
            0.835294117647059         0.243137254901961         0.309803921568627
            0.619607843137255       0.00392156862745098         0.258823529411765];
      case lower('nph_Parula')
        ccmap = [...
                         0        0.0745098039215686         0.352941176470588
         0.105098039215686         0.193725490196078         0.476078431372549
         0.210980392156863         0.290196078431373         0.544313725490196
         0.310588235294118         0.373333333333333          0.56078431372549
         0.404705882352941         0.443921568627451          0.56078431372549
         0.501960784313725         0.513725490196078         0.549019607843137
         0.603137254901961         0.581960784313726         0.512156862745098
         0.728627450980392         0.658039215686275         0.425098039215686
         0.864313725490196         0.729411764705882         0.296470588235294
         0.991372549019608         0.780392156862745          0.11843137254902
         0.984313725490196         0.819607843137255        0.0784313725490196];
     case lower('nph_Rainbow')
        ccmap = [...
                       0.2         0.431372549019608         0.690196078431373
         0.282352941176471         0.626274509803922          0.83921568627451
         0.507450980392157         0.784313725490196         0.854901960784314
         0.513725490196078         0.756862745098039          0.70156862745098
         0.509803921568627         0.709803921568627         0.349803921568627
         0.509803921568627         0.709803921568627         0.301960784313725
         0.621960784313725         0.749803921568627         0.337254901960784
         0.956862745098039         0.894509803921569         0.423529411764706
         0.972549019607843         0.803137254901961         0.405490196078431
         0.870588235294118         0.434117647058823         0.309803921568627
         0.882352941176471         0.545098039215686         0.513725490196078];
    case {lower('nph_BuOr'),lower('nph_BlueOrange')}
        ccmap = [...
         0.0313725490196078         0.188235294117647         0.419607843137255
         0.135455163587567         0.449725189198224         0.713326283812265
         0.290980392156863         0.594509803921569         0.789019607843137
         0.666666666666667         0.812156862745098         0.898823529411765
         0.843977272727273          0.90289592760181         0.959858696515043
                         1                         1                         1
         0.994494720965309         0.878849015590757         0.760544855110186
         0.992156862745098         0.722352941176471         0.484313725490196
         0.954509803921569                      0.44         0.106666666666667
         0.856555258467023         0.288205717830114       0.00589314194577355
         0.498039215686275         0.152941176470588        0.0156862745098039];
    case lower('nph_BuOr2')
        ccmap = [...
         0.0313725490196078         0.250980392156863         0.505882352941176
        0.0313725490196078         0.407843137254902         0.674509803921569
         0.168627450980392         0.549019607843137         0.745098039215686
          0.23503474718428         0.629359118440992          0.79966926529648
         0.305882352941176         0.701960784313725         0.827450980392157
         0.482352941176471                       0.8         0.768627450980392
         0.658823529411765         0.866666666666667         0.709803921568627
                       0.8          0.92156862745098         0.772549019607843
         0.841297613861865         0.938084181040837         0.814072637095465
          0.87843137254902         0.952941176470588         0.858823529411765
         0.968627450980392         0.988235294117647         0.941176470588235
                         1                         1         0.898039215686275
                         1         0.968627450980392         0.737254901960784
         0.998039215686275         0.936523524464454         0.657995937382857
         0.996078431372549         0.890196078431372         0.568627450980392
         0.996078431372549         0.768627450980392         0.309803921568627
         0.996078431372549                       0.6          0.16078431372549
         0.925490196078431          0.43921568627451        0.0784313725490196
         0.871378944461402         0.363483940693372        0.0302097583219334
                       0.8         0.298039215686275       0.00784313725490196
                       0.6         0.203921568627451        0.0156862745098039
                       0.4         0.145098039215686        0.0235294117647059];
     
    case lower('nph_BuOrRd')
        ccmap = [...
         0.0313725490196078         0.188235294117647         0.419607843137255
         0.135455163587567         0.449725189198224         0.713326283812265
         0.290980392156863         0.594509803921569         0.789019607843137
         0.666666666666667         0.812156862745098         0.898823529411765
         0.843977272727273          0.90289592760181         0.959858696515043
                         1                         1                         1
         0.994494720965309         0.887878787878788         0.732658819912354
         0.992156862745098         0.762745098039216         0.548235294117647
         0.947450980392157         0.427450980392157         0.295686274509804
          0.84848265776899          0.19888857901819         0.130279976588332
         0.498039215686275                         0                         0];
    case lower('nph_alt_jet')
        ccmap = [...
         0.529411764705882         0.850980392156863         0.992156862745098
         0.384313725490196         0.694117647058824         0.988235294117647
        0.0470588235294118          0.32156862745098          0.76078431372549
        0.0862745098039216         0.611764705882353         0.290196078431373
         0.113725490196078         0.741176470588235          0.32156862745098
         0.909803921568627         0.945098039215686         0.192156862745098
         0.996078431372549         0.925490196078431         0.203921568627451
         0.992156862745098         0.694117647058824         0.356862745098039
         0.988235294117647         0.317647058823529          0.12156862745098
         0.941176470588235        0.0470588235294118        0.0980392156862745
         0.788235294117647        0.0313725490196078        0.0745098039215686];
    case lower('nph_RdYlBuGrey')
        ccmap = [...
        0.985607843137255         0.985607843137255         0.985098039215686
         0.729411764705882         0.770196078431373         0.849411764705882
         0.564117647058824         0.648117647058824         0.765882352941176
         0.591372549019608         0.709019607843137         0.827372549019608
          0.76156862745098         0.843921568627451         0.856705882352941
         0.997647058823529         0.984313725490196         0.792941176470588
         0.987450980392157          0.80078431372549          0.43843137254902
         0.985882352941177         0.567058823529412         0.276862745098039
         0.947450980392157         0.257254901960784         0.184313725490196
         0.770196078431373       0.00549019607843137         0.167058823529412
          0.56078431372549       0.00235294117647059         0.171764705882353];
    case lower('nph_modspectral')
        ccmap = flipud([...
            153 55 69
            193 79 88
            222 127 120
            231 172 141
            240 212 169
            255 255 255
            215 234 235
            157 211 223
            91 191 218
            25 161 202
            14 107 162] ./ 255);
    case lower('nph_RdBuPastel')
        ccmap = flipud([...
            0.764705882352941         0.150708458565908                         0
            0.811764705882353         0.294117647058823         0.164705882352941
            0.854901960784314         0.466666666666667         0.366666666666667
            0.905882352941177         0.641176470588235         0.576470588235294
            0.949019607843137         0.811372549019608         0.776078431372549
            0.949019607843137         0.949019607843137         0.949019607843137
            0.803921568627451         0.892156862745098         0.921568627450981
            0.629411764705882         0.788235294117647         0.835294117647059
            0.443137254901961         0.678431372549019         0.749019607843137
            0.264705882352942         0.570588235294118         0.663529411764706
            0.0862745098039216         0.462745098039216         0.580392156862745]);
end

% interpolate to right length:

ccmapi = interp1(1:size(ccmap,1),ccmap,linspace(1,size(ccmap,1),len));

end% end function 


















