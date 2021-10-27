%Dragline Draw - Developed by The.Amazing (copyright 2021)
%By using this software you agree to not distribute it commercially, or
%pass it off as your own.
%Enjoy- this is version 1.0 (VERY UNREFINED)

clear all
[points,nvecs,nvecpoints] = stlread('Test_Contour.stl',45);
ZOFFSET = 0; %vertical offset (mm)
flipped = 0; %flip geometry horizontally? (0=No,1=Yes)

points(:,3) = points(:,3)+ZOFFSET;
nvecpoints(:,3) = nvecpoints(:,3)+ZOFFSET;
if flipped == 1
    points(:,1) = -points(:,1);
    nvecpoints(:,1) = -nvecpoints(:,1);
    nvecs(:,1) = -nvecs(:,1);
end

stepsize_contour = 0.13; %sampling distance between points for contoured infill, mm
stepsize_contour_border = 0.05; %sampling distance between points for contoured borders, mm

layerheight = 0.3; %non-contoured layer height, mm (for planar layers)
linewidth = 0.4; %nozzlewidth, mm
close_linespacing = 0.4; %spacing between close contoured lines, mm (not upper contours)

support_interface_offset = 0.3; %Gap between support and upper contoured layers (mm)

sampling_dist = 0.01; %final point spacing for contoured lines, mm (interpolated at end of generation)

upper_layers_flowfactor = 3.3; %flow rate multiplier for upper contoured layers
upper_layers_borderfactor = 4; %flow rate multiplier for upper contoured layer borders

flowfactor = 1.3; %flow rate multiplier for all other layers

stretch_down = 0; %1: stretch paths toward the edge of the XY-coordinate-region of the part, 0: stretch paths toward the bottom

clip_paths_every_iteration = 1; %1: trim paths every iteration (slower, but less errors), 0: trim paths after all iterations (faster)

support_temp = 230; %support material extruder temperature (Fahrenheit)
mesh_temp = 200; %upper layer material extruder temperature (Fahrenheit)

topcontour_linespacing = 1.2; %upper layers spacing between paths (mm)
num = 25; %number of samples for running average smoothing (more = more smooth, less accurate)
filamentD = 1.75; %filament diameter (mm)

middleX = mean(points(:,1));
middleY = mean(points(:,2));

stretchFX = max(points(:,1))+4;
stretchFY = max(points(:,2))+4;

stretchBX = min(points(:,1))-4;
stretchBY = min(points(:,2))-4;

lims = [min(points(:,1)),max(points(:,1));min(points(:,2)),max(points(:,2))];

clearZ = max(points(:,3))+4;

infillspacing = 3; %linespacing for infill/support material mm
skinlayer = 0; %number of layers of outer skin for the support material/planar layers
wall_lines = 1; %number of wall lines for the support material/planar layers
wallsmoothnum = 27; %number of samples for running average smoothing of walls for the support/planar layers  (more = more smooth, less accurate)

flatbottom = 1; %1: part sits flat on build plate, 0: part is upper layers only (use for mesh lens)
bordersmoothnum = 40; %number of smaples for running average smoothing for contoured borders (more = more smooth, less accurate)
contourborderlines = 3; %number of contoured border lines

contourlayerheight = 0.3; %contoured layers layer height, mm
contourthickness = 4; %total contoured thickness, mm
num_contourlayers = contourthickness/contourlayerheight; %number of contoured layers
num_topcontour = 2; %number of upper contoured layers (not support)
num_topborder = 2; %number of border lines in the upper contoured layers

contourlayers = cell(num_contourlayers,1);
layers = cell(num_contourlayers,1);
%tic
tic

for i = 1:num_contourlayers+1
    clf
    Lpoints = nvecpoints-(i-1)*contourlayerheight*nvecs;
    contourlayers{i} = [Lpoints nvecs];
end


pointspacing = stepsize_contour_border;
printshape = contourlayers{1};
[xq,yq] = meshgrid(min(printshape(:,1)):pointspacing:max(printshape(:,1)),min(printshape(:,2)):pointspacing:max(printshape(:,2)));
zq = griddata(printshape(:,1),printshape(:,2),printshape(:,3),xq,yq);

shapeX = reshape(xq,[],1);
shapeY = reshape(yq,[],1);
shapeZ = reshape(zq,[],1);
shapeX = shapeX(~isnan(shapeZ));
shapeY = shapeY(~isnan(shapeZ));
shapeZ = shapeZ(~isnan(shapeZ));

XYBI = boundary(shapeX,shapeY,1.0);

Xshapebound = shapeX(XYBI);
Yshapebound = shapeY(XYBI);

disp('contours complete')
%toc
%%
Rstop_cond = 0;
Lstop_cond = 0;
stepsize = stepsize_contour;
for i = 1:num_contourlayers
    if i > num_contourlayers-num_topcontour
        linespacing = close_linespacing;
    else
        linespacing = topcontour_linespacing;
    end
        
    contour = contourlayers{i};
    if mod(i,2)
        sortdimension = 2;
        perpdimension = 1;
        stretchF = stretchFY;
        stretchB = stretchBY;
    else
        sortdimension = 1;
        perpdimension = 2;
        stretchF = stretchFX;
        stretchB = stretchBX;
    end
    initnumline = min(contour(:,sortdimension)):stepsize:max(contour(:,sortdimension));
    otherdim = ones(size(initnumline)).*(min(contour(:,perpdimension))+max(contour(:,perpdimension)))./2;
    if mod(i,2)
        [initZ,initnX,initnY,initnZ] = STLSURF(otherdim',initnumline',contour(:,1:3),contour(:,4:6));
        initialpath = [otherdim',initnumline',initZ,initnX,initnY,initnZ];
    else
        [initZ,initnX,initnY,initnZ] = STLSURF(initnumline',otherdim',contour(:,1:3),contour(:,4:6));
        initialpath = [initnumline',otherdim',initZ,initnX,initnY,initnZ];
    end
    
    initialpath = sortrows(initialpath,sortdimension);
    
    initialpath = runningaverage(initialpath,num);  
    if clip_paths_every_iteration == 1
        insideshape = inpolygon(initialpath(:,1),initialpath(:,2),Xshapebound,Yshapebound);
        initialpath = initialpath(insideshape,:);
    end
    
    figure(1)
    clf
    plot3(contour(:,1),contour(:,2),contour(:,3),'c.')
    hold on
    axis equal
    zlim([0,20]);
    plot3(initialpath(:,1),initialpath(:,2),initialpath(:,3),'b-');
    pause(0.1)  
    disp('found initial path')
    
    sidevec = cross(initialpath(2:end,1:3)-initialpath(1:end-1,1:3),initialpath(1:end-1,4:6),2);
    sidepointsR=initialpath(1:end-1,1:3)+linespacing.*sidevec./sqrt(sidevec(:,1).^2+sidevec(:,2).^2+sidevec(:,3).^2);
    sidepointsL=initialpath(1:end-1,1:3)-linespacing.*sidevec./sqrt(sidevec(:,1).^2+sidevec(:,2).^2+sidevec(:,3).^2);
    
    [RZ,RnX,RnY,RnZ] = STLSURF(sidepointsR(:,1),sidepointsR(:,2),contour(:,1:3),contour(:,4:6));
    [LZ,LnX,LnY,LnZ] = STLSURF(sidepointsL(:,1),sidepointsL(:,2),contour(:,1:3),contour(:,4:6));
    newpathR = [sidepointsR(:,1),sidepointsR(:,2),RZ,RnX,RnY,RnZ];
    newpathL = [sidepointsL(:,1),sidepointsL(:,2),LZ,LnX,LnY,LnZ];
    
    newpathR = sortrows(newpathR,sortdimension);
    newpathL = sortrows(newpathL,sortdimension);
    newpathL = runningaverage(newpathL,num);
    newpathR = runningaverage(newpathR,num);
    
    if stretch_down == 1
        if Rstop_cond == 0
            if newpathR(end,3)>0
                disp('stretching right path forward')
                linepoints = newpathR(end-1:end,1:3);
                linediff = linepoints(2,:)-linepoints(1,:);
                offset = newpathR(end,3)-0;
                scale = offset/abs(linediff(3));
                endpoint = scale*linediff+linepoints(1,:);
                [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                newpathR = [newpathR;newline];
                %pause
            end
            if newpathR(1,3)>0
                disp('stretching right path backward')
                linepoints = newpathR(1:2,1:3);
                linediff = linepoints(1,:)-linepoints(2,:);
                offset = newpathR(1,3)-0;
                scale = offset/abs(linediff(3));
                endpoint = scale*linediff+linepoints(2,:);
                [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                newpathR = [newline;newpathR];
                %pause
            end
        end

        if Lstop_cond == 0
            if newpathL(end,3)>0
                disp('stretching left path forward')
                linepoints = newpathL(end-1:end,1:3);
                linediff = linepoints(2,:)-linepoints(1,:);
                offset = newpathL(end,3)-0;
                scale = offset/abs(linediff(3));
                endpoint = scale*linediff+linepoints(1,:);
                [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                newpathL = [newpathL;newline];
                %pause
            end
            if newpathL(1,3)>0
                disp('stretching left path backward')
                linepoints = newpathL(1:2,1:3);
                linediff = linepoints(1,:)-linepoints(2,:);
                offset = newpathL(1,3)-0;
                scale = offset/abs(linediff(3));
                endpoint = scale*linediff+linepoints(2,:);
                [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                newpathL = [newline;newpathL];
                %pause
            end
        end
    else
        if Rstop_cond == 0
            if newpathR(end,sortdimension)<stretchF
                disp('stretching right path forward')
                linepoints = newpathR(end-1:end,1:3);
                linediff = linepoints(2,:)-linepoints(1,:);
                offset = stretchF-newpathR(end,sortdimension);
                scale = offset/abs(linediff(sortdimension));
                endpoint = scale*linediff+linepoints(1,:);
                [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                newpathR = [newpathR;newline];
                %pause
            end
            if newpathR(1,sortdimension)>stretchB
                disp('stretching right path backward')
                linepoints = newpathR(1:2,1:3);
                linediff = linepoints(1,:)-linepoints(2,:);
                offset = newpathR(1,sortdimension)-stretchB;
                scale = offset/abs(linediff(sortdimension));
                endpoint = scale*linediff+linepoints(2,:);
                [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                newpathR = [newline;newpathR];
                %pause
            end
        end

        if Lstop_cond == 0
            if newpathL(end,sortdimension)<stretchF
                disp('stretching left path forward')
                linepoints = newpathL(end-1:end,1:3);
                linediff = linepoints(2,:)-linepoints(1,:);
                offset = stretchF-newpathL(end,sortdimension);
                scale = offset/abs(linediff(sortdimension));
                endpoint = scale*linediff+linepoints(1,:);
                [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                newpathL = [newpathL;newline];
                %pause
            end
            if newpathL(1,sortdimension)>stretchB
                disp('stretching left path backward')
                linepoints = newpathL(1:2,1:3);
                linediff = linepoints(1,:)-linepoints(2,:);
                offset = newpathL(1,sortdimension)-stretchB;
                scale = offset/abs(linediff(sortdimension));
                endpoint = scale*linediff+linepoints(2,:);
                [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                newpathL = [newline;newpathL];
                %pause
            end
        end
    end
    pathD = [sqrt((newpathR(2:end,1)-newpathR(1:end-1,1)).^2+(newpathR(2:end,2)-newpathR(1:end-1,2)).^2+(newpathR(2:end,3)-newpathR(1:end-1,3)).^2);0];
    while sum(pathD>stepsize.*2)>0
        skipI = find(pathD>stepsize.*2);
        for sI = 1:length(skipI)
            newrow = (newpathR(skipI(sI),:)+newpathR(skipI(sI)+1,:))./2;
            newpathR=[newpathR(1:skipI(sI),:);newrow;newpathR(skipI(sI)+1:end,:)];
        end
        pathD = [sqrt((newpathR(2:end,1)-newpathR(1:end-1,1)).^2+(newpathR(2:end,2)-newpathR(1:end-1,2)).^2+(newpathR(2:end,3)-newpathR(1:end-1,3)).^2);0];
    end

    pathD = [sqrt((newpathL(2:end,1)-newpathL(1:end-1,1)).^2+(newpathL(2:end,2)-newpathL(1:end-1,2)).^2+(newpathL(2:end,3)-newpathL(1:end-1,3)).^2);0];
    while sum(pathD>stepsize.*2)>0
        skipI = find(pathD>stepsize.*2);
        for sI = 1:length(skipI)
            newrow = (newpathL(skipI(sI),:)+newpathL(skipI(sI)+1,:))./2;
            newpathL=[newpathL(1:skipI(sI),:);newrow;newpathL(skipI(sI)+1:end,:)];
        end
        pathD = [sqrt((newpathL(2:end,1)-newpathL(1:end-1,1)).^2+(newpathL(2:end,2)-newpathL(1:end-1,2)).^2+(newpathL(2:end,3)-newpathL(1:end-1,3)).^2);0];
    end
    if clip_paths_every_iteration == 1
        insideshapeL = inpolygon(newpathL(:,1),newpathL(:,2),Xshapebound,Yshapebound);
        newpathL = newpathL(insideshapeL,:);
        insideshapeR = inpolygon(newpathR(:,1),newpathR(:,2),Xshapebound,Yshapebound);
        newpathR = newpathR(insideshapeR,:);
    end
    disp('found left and right paths')
    
    plot3(newpathL(:,1),newpathL(:,2),newpathL(:,3),'r-');
    plot3(newpathR(:,1),newpathR(:,2),newpathR(:,3),'m-');
    paths = {newpathL,initialpath,newpathR};
    pause(0.1)
    iteration = 1;
    stopcondition = 0;
    Rstop_cond = 0;
    Lstop_cond = 0;
    while (stopcondition ==0)
        
        disp(max([newpathL(:,3);newpathR(:,3)]));
        sidevecR = cross(newpathR(2:end,1:3)-newpathR(1:end-1,1:3),newpathR(1:end-1,4:6),2);
        sidevecL = cross(newpathL(2:end,1:3)-newpathL(1:end-1,1:3),newpathL(1:end-1,4:6),2);
        sidepointsR=newpathR(1:end-1,1:3)+linespacing.*sidevecR./sqrt(sidevecR(:,1).^2+sidevecR(:,2).^2+sidevecR(:,3).^2);
        sidepointsL=newpathL(1:end-1,1:3)-linespacing.*sidevecL./sqrt(sidevecL(:,1).^2+sidevecL(:,2).^2+sidevecL(:,3).^2);
        [RZ,RnX,RnY,RnZ] = STLSURF(sidepointsR(:,1),sidepointsR(:,2),contour(:,1:3),contour(:,4:6));
        [LZ,LnX,LnY,LnZ] = STLSURF(sidepointsL(:,1),sidepointsL(:,2),contour(:,1:3),contour(:,4:6));
        newpathR = [sidepointsR(:,1),sidepointsR(:,2),RZ,RnX,RnY,RnZ];
        newpathL = [sidepointsL(:,1),sidepointsL(:,2),LZ,LnX,LnY,LnZ];

        newpathR = sortrows(newpathR,sortdimension);
        newpathL = sortrows(newpathL,sortdimension);
        
        newpathL = runningaverage(newpathL,num);
        newpathR = runningaverage(newpathR,num);
        
        if stretch_down == 1
            if Rstop_cond == 0
                if newpathR(end,3)>0
                    disp('stretching right path forward')
                    linepoints = newpathR(end-1:end,1:3);
                    linediff = linepoints(2,:)-linepoints(1,:);
                    offset = newpathR(end,3)-0;
                    scale = offset/abs(linediff(3));
                    endpoint = scale*linediff+linepoints(1,:);
                    [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                    newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                    newpathR = [newpathR;newline];
                    %pause
                end
                if newpathR(1,3)>0
                    disp('stretching right path backward')
                    linepoints = newpathR(1:2,1:3);
                    linediff = linepoints(1,:)-linepoints(2,:);
                    offset = newpathR(1,3)-0;
                    scale = offset/abs(linediff(3));
                    endpoint = scale*linediff+linepoints(2,:);
                    [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                    newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                    newpathR = [newline;newpathR];
                    %pause
                end
            end
            
            if Lstop_cond == 0
                if newpathL(end,3)>0
                    disp('stretching left path forward')
                    linepoints = newpathL(end-1:end,1:3);
                    linediff = linepoints(2,:)-linepoints(1,:);
                    offset = newpathL(end,3)-0;
                    scale = offset/abs(linediff(3));
                    endpoint = scale*linediff+linepoints(1,:);
                    [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                    newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                    newpathL = [newpathL;newline];
                    %pause
                end
                if newpathL(1,3)>0
                    disp('stretching left path backward')
                    linepoints = newpathL(1:2,1:3);
                    linediff = linepoints(1,:)-linepoints(2,:);
                    offset = newpathL(1,3)-0;
                    scale = offset/abs(linediff(3));
                    endpoint = scale*linediff+linepoints(2,:);
                    [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                    newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                    newpathL = [newline;newpathL];
                    %pause
                end
            end
        else
            if Rstop_cond == 0
                if newpathR(end,sortdimension)<stretchF
                    disp('stretching right path forward')
                    linepoints = newpathR(end-1:end,1:3);
                    linediff = linepoints(2,:)-linepoints(1,:);
                    offset = stretchF-newpathR(end,sortdimension);
                    scale = offset/abs(linediff(sortdimension));
                    endpoint = scale*linediff+linepoints(1,:);
                    [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                    newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                    newpathR = [newpathR;newline];
                    %pause
                end
                if newpathR(1,sortdimension)>stretchB
                    disp('stretching right path backward')
                    linepoints = newpathR(1:2,1:3);
                    linediff = linepoints(1,:)-linepoints(2,:);
                    offset = newpathR(1,sortdimension)-stretchB;
                    scale = offset/abs(linediff(sortdimension));
                    endpoint = scale*linediff+linepoints(2,:);
                    [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                    newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                    newpathR = [newline;newpathR];
                    %pause
                end
            end
            
            if Lstop_cond == 0
                if newpathL(end,sortdimension)<stretchF
                    disp('stretching left path forward')
                    linepoints = newpathL(end-1:end,1:3);
                    linediff = linepoints(2,:)-linepoints(1,:);
                    offset = stretchF-newpathL(end,sortdimension);
                    scale = offset/abs(linediff(sortdimension));
                    endpoint = scale*linediff+linepoints(1,:);
                    [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                    newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                    newpathL = [newpathL;newline];
                    %pause
                end
                if newpathL(1,sortdimension)>stretchB
                    disp('stretching left path backward')
                    linepoints = newpathL(1:2,1:3);
                    linediff = linepoints(1,:)-linepoints(2,:);
                    offset = newpathL(1,sortdimension)-stretchB;
                    scale = offset/abs(linediff(sortdimension));
                    endpoint = scale*linediff+linepoints(2,:);
                    [addZ,addnX,addnY,addnZ] = STLSURF(endpoint(:,1),endpoint(:,2),contour(:,1:3),contour(:,4:6));
                    newline = [endpoint(:,1),endpoint(:,2),addZ,addnX,addnY,addnZ];
                    newpathL = [newline;newpathL];
                    %pause
                end
            end
        end
        
        pathD = [sqrt((newpathR(2:end,1)-newpathR(1:end-1,1)).^2+(newpathR(2:end,2)-newpathR(1:end-1,2)).^2+(newpathR(2:end,3)-newpathR(1:end-1,3)).^2);0];
        while sum(pathD>stepsize.*2)>0
            skipI = find(pathD>stepsize.*2);
            for sI = 1:length(skipI)
                newrow = (newpathR(skipI(sI),:)+newpathR(skipI(sI)+1,:))./2;
                newpathR=[newpathR(1:skipI(sI),:);newrow;newpathR(skipI(sI)+1:end,:)];
            end
            pathD = [sqrt((newpathR(2:end,1)-newpathR(1:end-1,1)).^2+(newpathR(2:end,2)-newpathR(1:end-1,2)).^2+(newpathR(2:end,3)-newpathR(1:end-1,3)).^2);0];
        end

        pathD = [sqrt((newpathL(2:end,1)-newpathL(1:end-1,1)).^2+(newpathL(2:end,2)-newpathL(1:end-1,2)).^2+(newpathL(2:end,3)-newpathL(1:end-1,3)).^2);0];
        while sum(pathD>stepsize.*2)>0
            skipI = find(pathD>stepsize.*2);
            for sI = 1:length(skipI)
                newrow = (newpathL(skipI(sI),:)+newpathL(skipI(sI)+1,:))./2;
                newpathL=[newpathL(1:skipI(sI),:);newrow;newpathL(skipI(sI)+1:end,:)];
            end
            pathD = [sqrt((newpathL(2:end,1)-newpathL(1:end-1,1)).^2+(newpathL(2:end,2)-newpathL(1:end-1,2)).^2+(newpathL(2:end,3)-newpathL(1:end-1,3)).^2);0];
        end
        if clip_paths_every_iteration == 1
            insideshapeL = inpolygon(newpathL(:,1),newpathL(:,2),Xshapebound,Yshapebound);
            newpathL = newpathL(insideshapeL,:);
            insideshapeR = inpolygon(newpathR(:,1),newpathR(:,2),Xshapebound,Yshapebound);
            newpathR = newpathR(insideshapeR,:);
        end
        newpathL = runningaverage(newpathL,5);
        newpathR = runningaverage(newpathR,5);
        plot3(newpathL(:,1),newpathL(:,2),newpathL(:,3),'r-');
        plot3(newpathR(:,1),newpathR(:,2),newpathR(:,3),'m-');
        paths = {newpathL,paths{1:end},newpathR};
        pause(0.1)
        disp('completed iteration')
        disp(iteration)
        iteration=iteration+1;
        if mean(newpathR(:,perpdimension))>mean(newpathL(:,perpdimension))
            upperpath = newpathR;
            upper = 'R';
            lowerpath = newpathL;
        else
            upperpath = newpathL;
            upper = 'L';
            lowerpath = newpathR;
        end
        if abs(min(upperpath(:,perpdimension))-lims(perpdimension,2))<close_linespacing*10
            if upper == 'R'
                Rstop_cond = 1;
            else
                Lstop_cond = 1;
            end
        end
         
        if abs(max(lowerpath(:,perpdimension))-lims(perpdimension,1))<close_linespacing*10
            if upper == 'R'
                Lstop_cond = 1;
            else
                Rstop_cond = 1;
            end
        end
        if Rstop_cond+Lstop_cond==2
            stopcondition = 1;
        end
%         if (max([newpathL(:,3);newpathR(:,3)])<layerheight/2)
%             stopcondition = 1;
%         end
    end
    layers{i} = paths;
%     pause
end
disp('initial contoured paths complete')
%%
if flatbottom == 1
    for i = 1:length(layers)
        paths = layers{i};
        for j = 1:length(paths)
            path = paths{j};
            path(path(:,3)<0,:) = [];
            paths{j} = path;
        end
        layers{i} = paths;
    end
end
%% 
borderlayers = cell(num_contourlayers,1);
linespacing = close_linespacing;
stepsize = stepsize_contour_border;
for Ci = 1:length(borderlayers)
    
    borderpaths = cell(1,contourborderlines+1);
    contour = contourlayers{Ci};
    layerforboundary = layers{Ci};
    if flatbottom == 0
        bi = boundary(contour(:,1),contour(:,2),1);
        initboundline = contour(bi,1:2);
    else
        bpoints = [];
        I = 1;
        for i = 1:length(layerforboundary)
            pathforboundary = layerforboundary{i};
            bpoints = [bpoints;pathforboundary(:,1:2)];
        end
        bi = boundary(bpoints(:,1),bpoints(:,2),0.1);
        initboundline = bpoints(bi,1:2);
    end
    [BZ,BnX,BnY,BnZ] = STLSURF(initboundline(:,1),initboundline(:,2),contour(:,1:3),contour(:,4:6));
    initboundpath = [initboundline(:,1),initboundline(:,2),BZ,BnX,BnY,BnZ];
    initboundpath(end+1,:) = initboundpath(1,:);
    pathD = [sqrt((initboundpath(2:end,1)-initboundpath(1:end-1,1)).^2+(initboundpath(2:end,2)-initboundpath(1:end-1,2)).^2+(initboundpath(2:end,3)-initboundpath(1:end-1,3)).^2);0];
    while max(pathD)>stepsize*2
        disp(length(initboundpath));
        [mx,mxi]=max(pathD);
        newpoint = 0.5.*(initboundpath(mxi+1,1:6)+initboundpath(mxi,1:6));
        initboundpath=[initboundpath(1:mxi,:);newpoint;initboundpath(mxi+1:end,:)];
        pathD = [sqrt((initboundpath(2:end,1)-initboundpath(1:end-1,1)).^2+(initboundpath(2:end,2)-initboundpath(1:end-1,2)).^2+(initboundpath(2:end,3)-initboundpath(1:end-1,3)).^2);0];
%         plot3(initboundpath(:,1),initboundpath(:,2),initboundpath(:,3),'r.','LineWidth',1);
        
    end
    initboundpath = runningaverage(initboundpath,bordersmoothnum);
    initboundpath(end+1,:) = initboundpath(1,:);
    while max(pathD)>stepsize*2
        disp(length(initboundpath));
        [mx,mxi]=max(pathD);
        newpoint = 0.5.*(initboundpath(mxi+1,1:6)+initboundpath(mxi,1:6));
        initboundpath=[initboundpath(1:mxi,:);newpoint;initboundpath(mxi+1:end,:)];
        pathD = [sqrt((initboundpath(2:end,1)-initboundpath(1:end-1,1)).^2+(initboundpath(2:end,2)-initboundpath(1:end-1,2)).^2+(initboundpath(2:end,3)-initboundpath(1:end-1,3)).^2);0];
%         plot3(initboundpath(:,1),initboundpath(:,2),initboundpath(:,3),'r.','LineWidth',1);
        
    end
    borderpaths{1,1} = initboundpath;
    newborderpath = initboundpath;
    figure(2)
    clf
%     plot3(contour(:,1),contour(:,2),contour(:,3),'c.')
    plot3(initboundpath(:,1),initboundpath(:,2),initboundpath(:,3),'b-','LineWidth',1);

    hold on
    axis equal
    
    disp('found initial border')
    for Pi = 2:contourborderlines+1
        sidevecs = cross(newborderpath(2:end,1:3)-newborderpath(1:end-1,1:3),newborderpath(1:end-1,4:6),2);
        sidepoints=newborderpath(1:end-1,1:3)-linespacing.*sidevecs./sqrt(sidevecs(:,1).^2+sidevecs(:,2).^2+sidevecs(:,3).^2);
        [BZ,BnX,BnY,BnZ] = STLSURF(sidepoints(:,1),sidepoints(:,2),contour(:,1:3),contour(:,4:6));
        newborderpath = [sidepoints(:,1),sidepoints(:,2),BZ,BnX,BnY,BnZ];
        
        newborderpath(end+1,:) = newborderpath(1,:);
        pathD = [sqrt((newborderpath(2:end,1)-newborderpath(1:end-1,1)).^2+(newborderpath(2:end,2)-newborderpath(1:end-1,2)).^2+(newborderpath(2:end,3)-newborderpath(1:end-1,3)).^2);0];
        while max(pathD)>stepsize*2
            disp(length(newborderpath));
            [mx,mxi]=max(pathD);
            newpoint = 0.5.*(newborderpath(mxi+1,1:6)+newborderpath(mxi,1:6));
            newborderpath=[newborderpath(1:mxi,:);newpoint;newborderpath(mxi+1:end,:)];
            pathD = [sqrt((newborderpath(2:end,1)-newborderpath(1:end-1,1)).^2+(newborderpath(2:end,2)-newborderpath(1:end-1,2)).^2+(newborderpath(2:end,3)-newborderpath(1:end-1,3)).^2);0];
        end
        newborderpath = runningaverage(newborderpath,bordersmoothnum);
        newborderpath(end+1,:) = newborderpath(1,:);
        while max(pathD)>stepsize*2
            disp(length(newborderpath));
            [mx,mxi]=max(pathD);
            newpoint = 0.5.*(newborderpath(mxi+1,1:6)+newborderpath(mxi,1:6));
            newborderpath=[newborderpath(1:mxi,:);newpoint;newborderpath(mxi+1:end,:)];
            pathD = [sqrt((newborderpath(2:end,1)-newborderpath(1:end-1,1)).^2+(newborderpath(2:end,2)-newborderpath(1:end-1,2)).^2+(newborderpath(2:end,3)-newborderpath(1:end-1,3)).^2);0];
        end
        plot3(newborderpath(:,1),newborderpath(:,2),newborderpath(:,3),'b-','LineWidth',1);
        
        borderpaths{Pi,1} = newborderpath;
        disp(Pi);
        
    end
    borderpaths(:,2:end) = [];
    borderlayers{Ci} = borderpaths;
    
end
%%

stepsize = 0.1;
for i = 1:length(layers)
    paths = layers{i};
    borderpaths = borderlayers{i};
    innerborder = borderpaths{end-1};
    cutboundary = innerborder(:,1:2);
    for j = 1:length(paths)
        path = paths{j};
        inbound = inpolygon(path(:,1),path(:,2),cutboundary(:,1),cutboundary(:,2));
        path(~inbound,:) = [];
        if mod(j,2)
            path = path(end:-1:1,:);
        end
        paths{j} = path;
    end
    pathsizes = cellsizes(paths);
    paths(pathsizes==0) = [];
    firstpath = paths{1};
    startpoint = firstpath(1,1:3);
    for j = length(borderpaths):-1:1
        path = borderpaths{j};
        dists = sqrt((path(:,1)-startpoint(1,1)).^2+(path(:,2)-startpoint(1,2)).^2+(path(:,3)-startpoint(1,3)).^2);
        [~,minI]=min(dists);
        path = [path(minI:end,:);path(1:minI+1,:)];
        startpoint = path(1,1:3);
        borderpaths{j} = path;
    end
    paths = paths(end:-1:1);
    paths(end+1:end+length(borderpaths)-1) = borderpaths(end-1:-1:1);
    paths = paths(end:-1:1);
    layers{i} = paths;
end
%% 
pointspacing = stepsize;
TwoDPrintMesh = contourlayers{end};
[xq,yq] = meshgrid(min(TwoDPrintMesh(:,1)):pointspacing:max(TwoDPrintMesh(:,1)),min(TwoDPrintMesh(:,2)):pointspacing:max(TwoDPrintMesh(:,2)));
zq = griddata(TwoDPrintMesh(:,1),TwoDPrintMesh(:,2),TwoDPrintMesh(:,3),xq,yq);
centerx = (min(TwoDPrintMesh(:,1))+max(TwoDPrintMesh(:,1)))./2;
rangex = max(TwoDPrintMesh(:,1))-min(TwoDPrintMesh(:,1));
centery = (min(TwoDPrintMesh(:,2))+max(TwoDPrintMesh(:,2)))./2;
rangey = max(TwoDPrintMesh(:,2))-min(TwoDPrintMesh(:,2));
XlineFW = centerx-0.55*rangex:stepsize:centerx+0.55*rangex;
XlineBW = centerx+0.55*rangex:-stepsize:centerx-0.55*rangex;
YlineFW = centery-0.55*rangey:stepsize:centery+0.55*rangey;
YlineBW = centery+0.55*rangey:-stepsize:centery-0.55*rangey;

CY = centery-0.55*rangey;
iY = 1;
infillX = [];
while CY < centery+0.55*rangey
    if mod(iY,2)
        newlines = [XlineFW',XlineFW'.*0+CY];
    else
        newlines = [XlineBW',XlineBW'.*0+CY];
    end
    infillX = [infillX;newlines];
    CY = CY + infillspacing;
    iY = iY+1;
end

CX = centerx-0.55*rangex;
iX = 1;
infillY = [];
while CX < centerx+0.55*rangex
    if mod(iX,2)
        newlines = [YlineFW'.*0+CX,YlineFW'];
    else
        newlines = [YlineBW'.*0+CX,YlineBW'];
    end
    infillY = [infillY;newlines];
    CX = CX + infillspacing;
    iX = iX+1;
end

CY = centery-0.55*rangey;
iY = 1;
solidX = [];
while CY < centery+0.55*rangey
    if mod(iY,2)
        newlines = [XlineFW',XlineFW'.*0+CY];
    else
        newlines = [XlineBW',XlineBW'.*0+CY];
    end
    solidX = [solidX;newlines];
    CY = CY + linespacing;
    iY = iY+1;
end

CX = centerx-0.55*rangex;
iX = 1;
solidY = [];
while CX < centerx+0.55*rangex
    if mod(iX,2)
        newlines = [YlineFW'.*0+CX,YlineFW'];
    else
        newlines = [YlineBW'.*0+CX,YlineBW'];
    end
    solidY = [solidY;newlines];
    CX = CX + linespacing;
    iX = iX+1;
end

% figure(6)
% clf
% plot(infillX(:,1),infillX(:,2),'b.');
% axis equal
% figure(7)
% clf
% plot(infillY(:,1),infillY(:,2),'r.');
% axis equal

pointsX = reshape(xq,[],1);
oX = min(pointsX);
pointsY = reshape(yq,[],1);
oY = min(pointsY);
pointsZ = reshape(zq,[],1);

currentheight = floor(max(pointsZ)/layerheight)*layerheight;
reglayers = {};
I = 1;
shapelayers = {};
figure(5)
clf
while currentheight>0
    hold off
    CS = rand(1,3);
    shapelayer = 0.*zq;
    shapelayer(zq>currentheight) = 1;
    %&(zq<currentheight+2.5*4*wall_lines*linespacing)
    shapelayer(isnan(shapelayer))=0;
    SE = strel('disk',5,8);
    wallSE = strel('disk',round(linespacing/pointspacing),8);
    pathSE = strel('disk',2,8);
    shapelayer = imdilate(shapelayer,SE);
    shapelayer = imerode(shapelayer,SE);
    shapelayers{I,1} = shapelayer;
    CC = bwconncomp(shapelayer);
    numshapes = length(CC.PixelIdxList);
%     disp(numshapes)
    shapes = cell(numshapes,1);
    walls = {};
    innerwalls = {};
    for si = 1:numshapes
        kind = CC.PixelIdxList{si};
        [row,col] = ind2sub(size(shapelayer),kind);
        shapemat = shapelayer.*0;
        shapemat(kind) = 1;
        for Wi = 1:wall_lines
            [row,col] = ind2sub(size(shapemat),find(imdilate(edge(shapemat,'nothinning'),pathSE)));
%             imshow(edge(shapemat,'nothinning'));
            ycoord = row.*pointspacing;
            xcoord = col.*pointspacing;
            bI = boundary(xcoord,ycoord,0.9);
            wall = [xcoord(bI),ycoord(bI),0.*ycoord(bI)+currentheight];
            if size(wall,1)>0
                wall(end+1,:) = wall(1,:);
                pathD = [sqrt((wall(2:end,1)-wall(1:end-1,1)).^2+(wall(2:end,2)-wall(1:end-1,2)).^2+(wall(2:end,3)-wall(1:end-1,3)).^2);0];
                while max(pathD)>stepsize*2
                    [mx,mxi]=max(pathD);
                    newpoint = 0.5.*(wall(mxi+1,1:3)+wall(mxi,1:3));
                    wall=[wall(1:mxi,:);newpoint;wall(mxi+1:end,:)];
                    pathD = [sqrt((wall(2:end,1)-wall(1:end-1,1)).^2+(wall(2:end,2)-wall(1:end-1,2)).^2+(wall(2:end,3)-wall(1:end-1,3)).^2);0];
                end
                wall = runningaverage(wall,5);
                wall(end+1,:) = wall(1,:);
                while max(pathD)>stepsize*2
                    [mx,mxi]=max(pathD);
                    newpoint = 0.5.*(wall(mxi+1,1:3)+wall(mxi,1:3));
                    wall=[wall(1:mxi,:);newpoint;wall(mxi+1:end,:)];
                    pathD = [sqrt((wall(2:end,1)-wall(1:end-1,1)).^2+(wall(2:end,2)-wall(1:end-1,2)).^2+(wall(2:end,3)-wall(1:end-1,3)).^2);0];
                end
                wall(:,1) = wall(:,1)+oX;
                wall(:,2) = wall(:,2)+oY;
                if Wi<=wall_lines
                    walls = [walls;{wall}];
                end
                shapemat = imerode(shapemat,wallSE);
                shapemat = imerode(shapemat,wallSE);
                shapemat = imdilate(shapemat,wallSE);
            end
        end
        innerwalls{si,1} = wall;
    end
    disp('walls done')
%     disp(length(innerwalls));
    totalshape = shapelayers{I};
    for Wi = 1:wall_lines
        totalshape = imerode(totalshape,wallSE);
        totalshape = imerode(totalshape,wallSE);
        totalshape = imdilate(totalshape,wallSE);
    end
    if ((currentheight>=skinlayer*layerheight)&&(currentheight<=max(pointsZ)-skinlayer*layerheight))
        covered = shapelayers{I-skinlayer};
        for Wi = 1:wall_lines
            covered = imerode(covered,wallSE);
            covered = imerode(covered,wallSE);
            covered = imdilate(covered,wallSE);
        end
        infilled = covered.*totalshape;
        exposed=totalshape-infilled;
        CCe = bwconncomp(exposed);
        CCi = bwconncomp(infilled);
        numshapesEx = length(CCe.PixelIdxList);
    else
        exposed=totalshape;
        infilled=0.*exposed;
    end
    solid = exposed.*0;
    infill = infilled.*0;
    solidPixels = linespacing/pointspacing;
    infillPixels = infillspacing/pointspacing;
    ZIGs = {};
    ZAGs = {};
    if mod(I,2)
        infill(1:infillPixels:end,:) = 1;
        solid(:,1:solidPixels:end) = 1;
        ZIGim = infill&infilled;
        ZAGim = solid&exposed;
        [ZIGx,ZIGy] = find(ZIGim');
        [ZAGy,ZAGx] = find(ZAGim);
        Xsteps = find([0;ZAGx(2:end)-ZAGx(1:end-1)]);
        for Fi = 1:2:length(Xsteps)-1
            revline = ZAGy(Xsteps(Fi):Xsteps(Fi+1)-1);
            revline = revline(end:-1:1,:);
            ZAGy(Xsteps(Fi):Xsteps(Fi+1)-1) = revline;
        end
        Ysteps = find([0;ZIGy(2:end)-ZIGy(1:end-1)]);
        for Fi = 1:2:length(Ysteps)-1
            revline = ZIGx(Ysteps(Fi):Ysteps(Fi+1)-1);
            revline = revline(end:-1:1,:);
            ZIGx(Ysteps(Fi):Ysteps(Fi+1)-1) = revline;
        end
    else
        infill(:,1:infillPixels:end) = 1;
        solid(1:solidPixels:end,:) = 1;
        ZIGim = infill&infilled;
        ZAGim = solid&exposed;
        [ZIGy,ZIGx] = find(ZIGim);
        [ZAGx,ZAGy] = find(ZAGim');
        Xsteps = find([0;ZIGx(2:end)-ZIGx(1:end-1)]);
        for Fi = 1:2:length(Xsteps)-1
            revline = ZIGy(Xsteps(Fi):Xsteps(Fi+1)-1);
            revline = revline(end:-1:1,:);
            ZIGy(Xsteps(Fi):Xsteps(Fi+1)-1) = revline;
        end
        Ysteps = find([0;ZAGy(2:end)-ZAGy(1:end-1)]);
        for Fi = 1:2:length(Ysteps)-1
            revline = ZAGx(Ysteps(Fi):Ysteps(Fi+1)-1);
            revline = revline(end:-1:1,:);
            ZAGx(Ysteps(Fi):Ysteps(Fi+1)-1) = revline;
        end
    end
    ZIG = [ZIGx.*pointspacing+oX ZIGy.*pointspacing+oY ZIGy.*0+currentheight];
    ZAG = [ZAGx.*pointspacing+oX ZAGy.*pointspacing+oY ZAGy.*0+currentheight];
    ZIGleft = ZIG;
    ZAGleft = ZAG;
    ZIGnew = ZIG;
    ZAGnew = ZAG;
    
    if size(wall,1)>0
        refpoint = wall(end,:);
    else
        refpoint = ZAG(1,:);
    end
%     for ZAi = 1:size(ZAG,1)
%         dists = sqrt((ZAGleft(:,1)-refpoint(:,2)).^2+(ZAGleft(:,1)-refpoint(:,2)).^2);
%         [~,zamI]=min(dists);
%         ZAGnew(ZAi,:) = ZAGleft(zamI,:);
%         refpoint = ZAGleft(zamI,:);
%         ZAGleft(zamI,:) = [];
%     end
%     for Zii = 1:size(ZIG,1)
%         dists = sqrt((ZIGleft(:,1)-refpoint(:,2)).^2+(ZIGleft(:,1)-refpoint(:,2)).^2);
%         [~,zimI]=min(dists);
%         ZIGnew(Zii,:) = ZIGleft(zimI,:);
%         refpoint = ZIGleft(zimI,:);
%         ZIGleft(zimI,:) = [];
%     end
    
    reglayers{I,1}=[walls;{ZAGnew};{ZIGnew}];
%         plot3(ZIG(:,1),ZIG(:,2),ZIG(:,3),'r-')
    disp(currentheight);
    currentheight=currentheight-layerheight;
    I = I+1;
    
end
reglayers = reglayers(end:-1:1,:);
%% 
figure(5)
clf
CS = rand(1,3);
for i = 1:length(reglayers)
    hold off
    axis equal
    
%     mesh(cX{i},cY{i},cZ{i});
    regpaths = reglayers{i};
    pathsizes = cellsizes(regpaths);
    regpaths(pathsizes==0) = [];
    reglayers{i} = regpaths;
    for j = 1:length(regpaths)
        regpath = regpaths{j};
        plot3(regpath(1,1),regpath(1,2),regpath(1,3),'r*')
        axis equal
        hold on
        grid on
        plot3(regpath(1:end/4,1),regpath(1:end/4,2),regpath(1:end/4,3),'.','Color',CS);
        pause(0.001)
        plot3(regpath(end/4+1:end/2,1),regpath(end/4+1:end/2,2),regpath(end/4+1:end/2,3),'.','Color',CS);
        pause(0.001)
        plot3(regpath(end/2+1:3*end/4,1),regpath(end/2+1:3*end/4,2),regpath(end/2+1:3*end/4,3),'.','Color',CS);
        pause(0.001)
        plot3(regpath(3*end/4+1:end,1),regpath(3*end/4+1:end,2),regpath(3*end/4+1:end,3),'.','Color',CS);
        disp(length(regpath));

%         if j>4
%             for k = 1:1:size(regpath,1)
%                 plot3(regpath(k,1),regpath(k,2),regpath(k,3),'.','Color',CS);axis equal
%                 axis([-40 100 -20 120 0 20]);
%                 hold on
%                 grid on
%                 pause(0.001);
%                 if CS(3)<1
%                     CS = CS+[0 0 0.01];
%                 else
%                     if CS(2)<1
%                         CS = CS+[0 0.01 0];
%                     else
%                         if CS(1)<1
%                             CS = CS+[0.01 0 0];
%                         end
%                     end
%                 end
%                 if CS(1)>1
%                     CS = [0 CS(2) CS(3)];
%                 end
%                 if CS(2)>1
%                     CS = [CS(1) 0 CS(3)];
%                 end
%                 if CS(3)>1
%                     CS = [CS(1) CS(2) 0];
%                 end
%             end
        
    end
    
end

%%
figure(4);
clf
for i = length(layers):-1:1
    
%     mesh(cX{i},cY{i},cZ{i});
    paths = layers{i};
    CS = rand(1,3);
    pathsizes = cellsizes(paths);
    paths(pathsizes==0) = [];
    layers{i} = paths;
    for j = 1:length(paths)
        path = paths{j};
        try
            path = splineinterpolate(path,sampling_dist);
        catch
            disp('unable to spline')
        end
        hold on
        axis equal
        grid on
        plot3(path(1,1),path(1,2),path(1,3),'r*')
        plot3(path(:,1),path(:,2),path(:,3),'Color',CS);
        pause(0.001)
%         plot3(path(end/4+1:end/2,1),path(end/4+1:end/2,2),path(end/4+1:end/2,3),'Color',CS);
%         pause(0.001)
%         plot3(path(end/2+1:3*end/4,1),path(end/2+1:3*end/4,2),path(end/2+1:3*end/4,3),'Color',CS);
%         pause(0.001)
%         plot3(path(3*end/4+1:end,1),path(3*end/4+1:end,2),path(3*end/4+1:end,3),'Color',CS);
%         disp(length(path));
%         pause(0.001)
    end
end

%% 
%GCODE GENERATION
output_filename = 'wacky.gcode';
outID = fopen(output_filename,'w');
printbeddims = [210 210];

print_feedrate = 3600;
travel_feedrate = 4500;
Z_feedrate = 300;

fprintf(outID,";Layer height: 0.2");
fprintf(outID,'\n');
fprintf(outID,";Generated with Dragline Draw");
fprintf(outID,'\n');
fprintf(outID,"M190 S60");
fprintf(outID,'\n');
fprintf(outID,strcat("M104 S",string(support_temp)));
fprintf(outID,'\n');
fprintf(outID,strcat("M109 S",string(support_temp)));
fprintf(outID,'\n');
fprintf(outID,"M82 ;absolute extrusion mode");
fprintf(outID,'\n');
fprintf(outID,"G21 ;metric values");
fprintf(outID,'\n');
fprintf(outID,"G90 ;absolute positioning");
fprintf(outID,'\n');
fprintf(outID,"M82 ;set extruder to absolute mode");
fprintf(outID,'\n');
fprintf(outID,"M107 ;start with the fan off");
fprintf(outID,'\n');
fprintf(outID,"G28 X0 Y0 ;move X/Y to min endstops");
fprintf(outID,'\n');
fprintf(outID,"G28 Z0 ;move Z to min endstops");
fprintf(outID,'\n');
fprintf(outID,"G1 Z15.0 F9000 ;move the platform down 15mm");
fprintf(outID,'\n');
fprintf(outID,"G92 E0 ;zero the extruded length");
fprintf(outID,'\n');
fprintf(outID,"G1 F200 E3 ;extrude 3mm of feed stock");
fprintf(outID,'\n');
fprintf(outID,"G92 E0 ;zero the extruded length again");
fprintf(outID,'\n');
fprintf(outID,"G1 F9000");
fprintf(outID,'\n');
fprintf(outID,"G5");
fprintf(outID,'\n');
fprintf(outID,"G92 E0");
fprintf(outID,'\n');
fprintf(outID,"G1 F1500 E-6.5");
fprintf(outID,'\n');

epoint = 0;

trav = 0;

fprintf(outID,strcat("G1 F",string(print_feedrate)));
fprintf(outID,'\n');

factor = flowfactor;

for i = 1:length(reglayers)
    paths = reglayers{i};
    for j = 1:length(paths)
        path = paths{j};
        for k = 1:size(path,1)
            xpoint = path(k,1)-middleX+printbeddims(1)./2;
            ypoint = path(k,2)-middleY+printbeddims(2)./2;
            zpoint = path(k,3);
            if k == 1
                feedrate = travel_feedrate;
            else
                diststep = sqrt((xpoint-oldxpoint)^2+(ypoint-oldypoint)^2+(zpoint-oldzpoint)^2);
                if diststep<stepsize*2
                    Vstep = diststep*layerheight*linespacing;
                    estep = 3*Vstep/(pi*filamentD^2);
                    epoint = epoint+estep.*factor;
                    feedrate = print_feedrate;
                else
                    feedrate = travel_feedrate;
                    trav = 1;
                    
                    fprintf(outID,strcat("G1 F",string(travel_feedrate)));
                    fprintf(outID,'\n');
                end
            end
            newline = strcat("G1 X",string(xpoint)," Y",string(ypoint)," Z",string(zpoint)," E",string(epoint));
            if trav == 1
                fprintf(outID,strcat("G1 F",string(print_feedrate)));
                fprintf(outID,'\n');
                trav = 0;
            end

            fprintf(outID,newline);
            fprintf(outID,'\n');
            oldxpoint = xpoint;
            oldypoint = ypoint;
            oldzpoint = zpoint;
        end
    end
    fprintf(outID,strcat(";LAYER:",string(i)));
    fprintf(outID,'\n');
    disp(i);
end


fprintf(outID,strcat("G1 Z",string(clearZ)," F",string(Z_feedrate)));
fprintf(outID,'\n');

layernum = 1;


fprintf(outID,strcat("G1 F",string(print_feedrate)));
fprintf(outID,'\n');
printborder = 1;
for i = length(layers):-1:1
    if i == length(layers)-num_topcontour
        fprintf(outID,strcat("M104 S",string(mesh_temp)));
        fprintf(outID,'\n');
        fprintf(outID,strcat("M109 S",string(mesh_temp)));
        fprintf(outID,'\n');
        fprintf(outID,strcat("G4 P60000"));
        fprintf(outID,'\n');
        factor = upper_layers_borderfactor;
    end
    if i == num_topcontour-num_topborder
        printborder = 0;
    end
    paths = layers{i};
    if printborder == 1
        startpathindex = 1;
    else
        startpathindex = contourborderlines;
    end
    for j = startpathindex:length(paths)
        if j>contourborderlines
            if i <= length(layers)-num_topcontour
                factor = upper_layers_flowfactor;
                alt_offset = support_interface_offset;
            else
                factor = flowfactor;
                alt_offset = 0;
            end
        else
            if i <= length(layers)-num_topcontour
                factor = upper_layers_borderfactor;
                alt_offset = support_interface_offset;
            else
                factor = flowfactor;
                alt_offset = 0;
            end
        end
        path = paths{j};
        for k = 1:size(path,1)
            if (k == 1)
                fprintf(outID,strcat("G1 X",string(path(k,1)-middleX+printbeddims(1)./2)," Y",string(path(k,2)-middleY+printbeddims(2)./2)," F",string(travel_feedrate)));
                fprintf(outID,'\n');
            end
            xpoint = path(k,1)-middleX+printbeddims(1)./2;
            ypoint = path(k,2)-middleY+printbeddims(2)./2;
            zpoint = path(k,3);
            if k == 1
                feedrate = travel_feedrate;
            else
                diststep = sqrt((xpoint-oldxpoint)^2+(ypoint-oldypoint)^2+(zpoint-oldzpoint)^2);
                if diststep<stepsize*100
                    Vstep = diststep*layerheight*linespacing;
                    estep = 3*Vstep/(pi*filamentD^2);
                    epoint = epoint+estep.*factor;
                    feedrate = print_feedrate;
                else
                    feedrate = travel_feedrate;
                end
            end
            newline = strcat("G1 X",string(xpoint)," Y",string(ypoint)," Z",string(zpoint+alt_offset)," E",string(epoint));
            fprintf(outID,newline);
            fprintf(outID,'\n');
            oldxpoint = xpoint;
            oldypoint = ypoint;
            oldzpoint = zpoint;
        end
    end
    fprintf(outID,"G1 E-1 F300");
    fprintf(outID,'\n');
    fprintf(outID,strcat("G1 Z",string(clearZ)," F",string(Z_feedrate)));
    fprintf(outID,'\n');
    fprintf(outID,strcat(";LAYER:",string(layernum+length(reglayers))));
    fprintf(outID,'\n');
    disp(layernum);
    layernum = layernum+1;
end
fprintf(outID,"G1 X0 Y0");
fprintf(outID,'\n');
fprintf(outID,"M104 S0");
fprintf(outID,'\n');
fprintf(outID,"M140 S0");
fprintf(outID,'\n');
fprintf(outID,"M84 ;steppers off");
fprintf(outID,'\n');
fprintf(outID,"M107");
fprintf(outID,'\n');

fclose(outID);

function [outpath] = runningaverage(inpath,num)
    pathsum = 0.*inpath(1:end-num+1,:);
%     disp(length(inpath))
%     disp(length(pathsum))
    for i = 1:num
        pathsum = pathsum+inpath(i:end-num+i,:);
    end
    pathavg = pathsum./num;
    outpath = pathavg;
%     endnum = ceil((num-1)./2);
%     disp(endnum)
%     startnum = floor((num-1)./2);
%     disp(startnum)
%     outpath = [inpath(1:startnum,:);pathavg;inpath(end-endnum,:)];
%     if length(outpath)~=length(inpath)
%         disp(length(inpath));
%         disp(length(outpath));
%     end
end

function [Z,nX,nY,nZ] = STLSURF(X,Y,centroids,normals)
    nX = zeros(size(X));
    nY = zeros(size(X));
    nZ = zeros(size(X));
    cX = zeros(size(X));
    cY = zeros(size(X));
    cZ = zeros(size(X));
    Z = zeros(size(X));
    xcentroids = centroids(:,1);
    ycentroids = centroids(:,2);
    zcentroids = centroids(:,3);
    if length(X)~=length(Y)
        error("first two arguments must have same length");
    end
    for i = 1:length(X)
        [~,dI]=min(sqrt((X(i)-xcentroids).^2+(Y(i)-ycentroids).^2));
        nX(i,1) = normals(dI,1);
        nY(i,1) = normals(dI,2);
        nZ(i,1) = normals(dI,3);
        cX(i,1) = xcentroids(dI);
        cY(i,1) = ycentroids(dI);
        cZ(i,1) = zcentroids(dI);
        D=nX(i,1)*cX(i,1)+nY(i,1)*cY(i,1)+nZ(i,1)*cZ(i,1);
        Z(i,1) = (D-nX(i,1)*X(i)-nY(i,1)*Y(i))/nZ(i,1);
    end
end

function sizes = cellsizes(CELL)
    sizes = zeros(size(CELL));
    for i = 1:length(CELL)
        sizes(i) = size(CELL{i},1);
    end
end

function pathout = splineinterpolate(pathin,dist)
    spline_obj = cscvn(pathin');
    xcoeffs = spline_obj.coefs(1:6:end,:);
    ycoeffs = spline_obj.coefs(2:6:end,:);
    zcoeffs = spline_obj.coefs(3:6:end,:);
    nxcoeffs = spline_obj.coefs(4:6:end,:);
    nycoeffs = spline_obj.coefs(5:6:end,:);
    nzcoeffs = spline_obj.coefs(6:6:end,:);
    tbreaks = spline_obj.breaks;
    originaldists = [0;sqrt((pathin(2:end,1)-pathin(1:end-1,1)).^2+(pathin(2:end,2)-pathin(1:end-1,2)).^2+(pathin(2:end,3)-pathin(1:end-1,3)).^2)];
    totaldist = sum(originaldists);
    numcurves = round(totaldist./dist);
    tstep = (max(tbreaks)-min(tbreaks))/numcurves;
    t = min(tbreaks):tstep:max(tbreaks);
    pathout = zeros(length(t),6);
    for ti = 1:length(t)
        tdiff = tbreaks-t(ti);
        curvei = max(find(tdiff<=0));
        if curvei > length(xcoeffs)
            curvei = length(xcoeffs);
        end
        T = t(ti)-tbreaks(curvei);
        Tpoly = [T.^3,T.^2,T,1];
        pathout(ti,1) = sum(xcoeffs(curvei,:).*Tpoly);
        pathout(ti,2) = sum(ycoeffs(curvei,:).*Tpoly);
        pathout(ti,3) = sum(zcoeffs(curvei,:).*Tpoly);
        pathout(ti,4) = sum(nxcoeffs(curvei,:).*Tpoly);
        pathout(ti,5) = sum(nycoeffs(curvei,:).*Tpoly);
        pathout(ti,6) = sum(nzcoeffs(curvei,:).*Tpoly);
        disp(curvei)
    end
end