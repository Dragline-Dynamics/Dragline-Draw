# Dragline-Draw
A multipurpose contoured 3D-printing slicer, developed for the creation of mesh lenses.

VERSION 1
A MATLAB script, several bugs to be worked out.

In the first section of the code, put in settings:
```MATLAB
ZOFFSET = 0; %vertical offset (mm)
flipped = 0; %flip geometry horizontally? (0=No,1=Yes)

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
```
