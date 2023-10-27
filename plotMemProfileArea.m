%% Plot membrane profiles colored according to lambda
%
%   Julian Hassinger
%   Biophysics Graduate Group
%   George Oster Lab
%   University of California, Berkeleyy
%
%   Copyright 2016
%
%   Last Edited: 1/28/2016
%
%%

% Input(s):
%   Sol - solution array
%   t - area mesh points
%   R0 - nondimensionalization length
%   coatArea - area start and end points for regions coated with isotropic curvature generators, entered as array, e.g. [0 1; 2 3]
%   actArea - area start and end points for regions with applied forces, entered as array, e.g. [0 1; 2 3]
%   barArea - area start and end points for regions coated with anisotropic curvature generators, entered as array, e.g. [0 1; 2 3]
%   xLim - min and max limits for horizontal axis of the plot, can be left empty
%   yLim - min and max limits for vertical axis of the plot, can be left empty
%   plotTitle - title for the plots, can be left empty
%   dash - 1 for dashed lines
%   varargin
%       xLabelOn - turns the x label on (1) and off (0)
%       yLabelOn - turns the y label on (1) and off (0)
%       xTickLabelOn - turns the x ticks label on (1) and off (0)
%       yTickLabelOn - turns the y ticks label on (1) and off (0)


function plotMemProfileArea(Sol, t, R0, coatArea, actArea, barArea, xLim, yLim, plotTitle, dash, varargin)

%figure        % create a new figure
    
% set visual properties of the plot
fontsize = 25;  %54;  %70;
lineWidth = 6; %6;  %12;
axesWidth = 5;

% initializes labels and tick labels to on
xLabelOn = 1;
yLabelOn = 1;
xTickLabelOn = 1;
yTickLabelOn = 1;

% sets dash to 0 if empty
if isempty(dash)
    dash = 0;
end    

% handle arbitrary number of inputs
if (~isempty(varargin))
    for ii = 1:length(varargin)
        switch ii

            % first input sets xLabelOn
            case 1
            if (~isempty(varargin{1}))
                xLabelOn = varargin{1};
            end
            
            % second input sets yLabelOn
            case 2
            if (~isempty(varargin{2}))
                yLabelOn = varargin{2};
            end
            
            % third input sets xTickLabelOn
            case 3
            if (~isempty(varargin{3}))
                xTickLabelOn = varargin{3};
            end
            
            % fourth input sets yTickLabelOn
            case 4
            if (~isempty(varargin{4}))
                yTickLabelOn = varargin{4};
            end
            
        end
    end
end

% set the colors of the membrane, coat, forces, and anisotropic curvature

%memColor = 'black'; 
%memColor = [0.9967    0.7816    0.2007];
memColor = [0.9139    0.7258    0.3063];
coatColor = 'blue';
%coatColor = [0    0.4470    0.7410];
actColor = 'red';
%actColor = [0.8500    0.3250    0.0980];
%barColor = 'green';
barColor = [0    0.6    0];

% plot the membrane shape
if dash == 1
    
    plot(Sol(1,:)*R0, Sol(2,:)*R0, ':', -Sol(1,:)*R0, Sol(2,:)*R0, ':', 'Color', memColor, 'LineWidth', lineWidth);
    
else
    
    plot(Sol(1,:)*R0, Sol(2,:)*R0, -Sol(1,:)*R0, Sol(2,:)*R0,'Color', memColor, 'LineWidth', lineWidth);

end

hold on

% plots the isotropic-curvature-generating coat
if ~isempty(coatArea)
    % loops over multiple regions
    for ii = 1:size(coatArea,1)

        aMin = coatArea(ii,1);
        aMax = coatArea(ii,2);
        
        if dash == 1
            
            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ':', ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ':', 'Color', coatColor, 'LineWidth', lineWidth);
            
        else

            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, 'Color', coatColor, 'LineWidth', lineWidth);
            
        end

    end
end

% plots the regions of applied force
if ~isempty(actArea)
    % loops over multiple regions
    for ii = 1:size(actArea,1)

        aMin = actArea(ii,1);
        aMax = actArea(ii,2);
        
        if dash == 1
            
            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, '.', ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, '.', 'Color', actColor, 'LineWidth', lineWidth);
            
        else

            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, '.', ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, '.', 'Color', actColor, 'LineWidth', lineWidth);
        
        end

    end
end

% plots the anisotropic-curvature-generating coat
if ~isempty(barArea)
    % loops over multiple regions
    for ii = 1:size(barArea,1)

        aMin = barArea(ii,1);
        aMax = barArea(ii,2);
        
        if dash == 1
            
            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ':', ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ':', 'Color', barColor, 'LineWidth', lineWidth);
            
        else

            plot(Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, ...
                -Sol(1,t>=aMin & t<=aMax)*R0, Sol(2,t>=aMin & t<=aMax)*R0, 'Color', barColor, 'LineWidth', lineWidth);
        
        end

    end
end

% adds x label
if xLabelOn == 1
    xlabel('R (nm)', 'FontSize',fontsize, 'FontName', 'Helvetica');
end

% adds y label
if yLabelOn == 1
    ylabel({'Z (nm)'}, 'FontSize',fontsize, 'FontName', 'Helvetica');
end

% sets plot aesthetic properties
set(gca,'FontSize',fontsize, 'FontName', 'Helvetica', 'XMinorTick', 'on', 'YMinorTick', 'on', 'Linewidth', axesWidth);

% turns off x tick label
if xTickLabelOn == 0
    set(gca, 'XTickLabel', []);
end

% turns off y tick label
if yTickLabelOn == 0
    set(gca, 'YTickLabel', []);
end

% sets plots title
if ~isempty(plotTitle)
    title(plotTitle, 'FontSize', fontsize+4, 'FontName', 'Helvetica');
end

% sets x limits of plot
if ~isempty(xLim)
    xlim(xLim);
end

% sets y limits of plot
if ~isempty(yLim)
    ylim(yLim);
end

% sets the aspect ratio of the plot to equalize lengths
daspect([1 1 1]);


hold off

end