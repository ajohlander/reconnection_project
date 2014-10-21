function out = c_4_v_timing_mva(x1,x2,x3,x4,R,column)
%C_4_V_TIMING_MVA Performs timing and minimum variance analysis on four s/c
%field data with a graphical user interface.
%   C_4_V_TIMING_MVA(b1,b2,b3,b4,R,column) interactive discontinuity
%   analysis analyzer on magnetic field b1,...b4 with position R using
%   column number 'column'. R has the form R.R1,...R.R4.
%   C_4_V_TIMING_MVA('B?',R) uses B1, B2, B3 and B4 from workspace
%   C_4_V_TIMING_MVA('B?',R, column) also uses column. 2=x, 3=y, 4=z.
%
%   The user is asked to define an interval in which the discontinuity
%   analysis should be made. Two new figures open when this is done.
%   In the first window, the normal vectors obtained from MVA and timing is
%   shown in 3D space. Normal vectors that do not fulfill l2/l3>5 is shown
%   as dotted arrows. The dotted arrows are not used in determining the
%   maximum angle.
%   In the second figure, the magnetic field is plotted in the
%   LMN-system, where L=maximum, M=intermiediate and N = minimum.
%
%   A prompt asks the user to save the variables: velocity of discontinuity V,
%   uncertainty in velocity dV, all 5 normal vectors in n, where
%   n.nTiming is the normal vector from timing and n.n1,...n4 is normal
%   vectors obtained from MVA, all vectors from MVA v, and eigenvalues from
%   MVA. The variables do not appear until the main figure is closed.



% Getting parameters
if(nargin<=3 && ischar(x1)), % either action as parameter or string variable
    if strfind(x1,'?')
        var_str = x1;
        
        evalin('base',['if ~exist(''' irf_ssub(var_str,1) '''), c_load(''' var_str ''');end' ]);
        c_eval('b?=evalin(''base'',irf_ssub(var_str,?));');
    end
    
    if nargin == 2
        column = 2;
    elseif nargin == 3
        R = x2;
        column = x3;
    end;
    
elseif   (nargin ==4) || (nargin == 5)
    if nargin==4
        irf.log('notice','Using second column');column=2;
    end
    b1 = x1;
    b2 = x2;
    b3 = x3;
    b4 = x4;
    
end

tint = [b1(1,1),b1(end,1)];


%---------------GUI-------------------
fGUI = irf_plot(4,'newfigure');
set(gcf,'PaperUnits','centimeters')
xSize = 15; ySize = 15;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*50 ySize*50])
set(gcf,'paperpositionmode','auto') % to get the same printing as on screen
clear xLeft xSize sLeft ySize yTop

h = zeros(1,4);
h(1) = irf_panel('C1');
h(2) = irf_panel('C2');
h(3) = irf_panel('C3');
h(4) = irf_panel('C4');
hold(h(1))
hold(h(2))
hold(h(3))
hold(h(4))

%Bz
pb = zeros(1,4);
pb(1) = plot(h(1),b1(:,1),b1(:,column));
pb(2) = plot(h(2),b2(:,1),b2(:,column));
pb(3) = plot(h(3),b3(:,1),b3(:,column));
pb(4) = plot(h(4),b4(:,1),b4(:,column));

%line at 0
fLine = plot(h(1),tint,[0,0]);
set(fLine,'Color',[.5 .5 .5])
fLine = plot(h(2),tint,[0,0]);
set(fLine,'Color',[.5 .5 .5])
fLine = plot(h(3),tint,[0,0]);
set(fLine,'Color',[.5 .5 .5])
fLine = plot(h(4),tint,[0,0]);
set(fLine,'Color',[.5 .5 .5])

%Meta
xlim(h(1),tint)
xlim(h(2),tint)
xlim(h(3),tint)
xlim(h(4),tint)


ystr = '';

switch column
    case 2
        ystr = 'B_{x} [nT]';
    case 3
        ystr = 'B_{y} [nT]';
    case 4
        ystr = 'B_{z} [nT]';
end

ylabel(h(1),ystr,'FontSize',16)
ylabel(h(2),ystr,'FontSize',16)
ylabel(h(3),ystr,'FontSize',16)
ylabel(h(4),ystr,'FontSize',16)

irf_timeaxis(h)

% Loop for the GUI
while true
    % Makes the GUI-window the current figure
    nFig = get(fGUI,'Parent');
    figure(nFig{1})
    
    clear x % For good measure
    try
        [x,~] = ginput(2);
    catch % If user closes window, the program closes
        disp('Figure closed. Any saved variables can be found in the Workspace');
        return
    end
    
    % Remove any existing colored area
    if (exist('a','var')==1)
        delete(a)
    end
    % The s/c that was clicked is the reference s/c
    clickedAx = gca;
    nCluster = find(h == clickedAx);
    
    b = 0;
    switch nCluster
        case 1
            b = b1;
        case 2
            b = b2;
        case 3
            b = b3;
        case 4
            b = b4;
    end
    
    % Finds the time interval from user data
    clickInd1 = find_closest_index(x(1),b(:,1));
    clickInd2 = find_closest_index(x(2),b(:,1));
    
    M = [clickInd1, clickInd2];
    
    if(M(1)>M(2))
        M = fliplr(M);
    end
    
    absB = sqrt(sum(b(:,2:4)'.^2));
    maxB = max(absB);
    
    %Draws a colored area that marks the time interval chosen
    a = zeros(1,4);
    for i = 1:4
        a(i) = area(h(i),[b(M(1,1),1) b(M(1,2),1)], [maxB maxB],-maxB, 'FaceColor', [0.5,1,0.5]);
        set(a(i),'EdgeColor','none')
        uistack(pb(i),'top')
    end
    
    
    % Calls for the timing function
    [V,dV,n,v,l] = disc_timing(b1,b2,b3,b4,R,M,nCluster);
    
    % Asks the user to save the variables to the Workspace.
    % The function does not return any variables due to the loop structure.
    checkLabels = {'Save velocity to variable named:' ...
        'Save dV to variable named:' 'Save normal vectors to variable named:'...
        'Save MVA matrices to variable named:' 'Save MVA eigenvalues to variable named:'};
    varNames = {'V','dV','n','v','l'};
    items = {V,dV,n,v,l};
    export2wsdlg(checkLabels,varNames,items,...
        'Save outputs to Workspace');
    
    out = V;
end

end

