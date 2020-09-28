function nightfury(varargin)
% NIGHTFURY Sets the style of the current figure to the nightfury settings

aspect = 1./sqrt(2);
width = 29.7;
height = width * aspect;
discrete = [172, 223, 245;
    42, 174, 215;
    14, 120, 189;
    22, 59, 125;
    10, 9, 51;
    80, 14, 51;
    154, 54, 90;
    180, 62, 83;
    234, 101, 81;
    241, 133, 93;
    255, 171, 118;
    255, 199, 129;
    255, 232, 168] ./ 255;

map = lithiumIon.magma();

background = 0.2 .* [1 1 1];   

if strcmp(varargin, 'subplot')
                
    set(gca,'FontSize',14, ...
    'color',background, ...
    'ColorOrder', discrete)
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    set(gcf,'units','centimeter',...
        'position',[1.5,1.5,width*1.6,height*1.3], ...
        'color',background)
    set(gca, 'FontName', 'Helvetica')

else
    set(gca,'FontSize',28, ...
        'color',background, ...
        'ColorOrder', discrete)
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    set(gcf,'units','centimeter',...
        'position',[5,5,width,height], ...
        'color',background)
    set(gca, 'FontName', 'Helvetica')
end                  

end