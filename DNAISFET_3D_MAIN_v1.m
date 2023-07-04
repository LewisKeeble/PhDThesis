function out = model

%% FLAGS
 
% Set flag to choose simulation type. 0 = no surface reactions, no EDL. 1 =
% surface reactions, no EDL. 2 = surface reactions, EDL
Simulation_Flag = 1;

%Flag to indicate whether the electrodes are 'activated' and DNA is attracted
%to electrode trapping regions
%0 = no activation, 1 = activation
ElectrodeActivation_Flag = 0;

%Flag to determine if there are square castellations added to the line
%electrodes.
%Castellation_flag = 0 if no castellations, = 1 for castellations
Castellation_Flag = 0;

%Flag that selects how trapped DNA is determined. 0 = all particles
%trapped, 1 = percentage with smallest z coordinate trapped, 2 = trapping region defined, 3 csv input.
Trapping_Flag = 2;

%% USER INPUT VARIABLES

% REQUIRED FOR ALL SIMULATION TYPES

% Enter number of DNA molecules in solution
DNACluster_N = 1000

% Starting time to each cluster is Normal distribution...
mu_start_time = 13 % in minutes
sigma_start_time = 3 % in minutes

% Create starting time dist (and sample from it)
pd = makedist('Normal', 'mu', mu_start_time, 'sigma', sigma_start_time)
sigmoid_start = 60 * sort(random(pd, 1, DNACluster_N))

% Adjust sigmoid_start so all times are relative to first time
sigmoid_start_relative = sigmoid_start - sigmoid_start(1)

ISFET_Nx = 64
ISFET_Ny = 64
ISFET_Width = 9.3e-6
ISFET_Length = 10.3e-6
Electrode_Width = 0
ISFETElectrode_Separation = 1e-6
ISFET_StartSeparation = 2e-6
ISFET_EndSeparation = 2e-6
ISFETISFET_Separation = 2e-6

Electrode_N = ISFET_Nx + 1

ElectrodeElectrode_Separation = 2*ISFETElectrode_Separation + ISFET_Width

Chip_XSize = ISFET_Nx*ElectrodeElectrode_Separation + Electrode_N*Electrode_Width

Chip_YSize = ISFET_StartSeparation + ISFET_Ny*ISFET_Length + (ISFET_Ny-1)*(ISFETISFET_Separation) + ISFET_EndSeparation

DNACluster_StartRadius = 2.5e-6;

% Add distance between clusters and walls to ease meshing
ClusterWall_Tolerance = 2*DNACluster_StartRadius;

SolutionHeight = 3e-3

ChipWall_Separation_Xpos = (1.75e-3 + 2*0.262e-3 - Chip_XSize)/2;
ChipWall_Separation_Xneg = ChipWall_Separation_Xpos;
ChipWall_Separation_Ypos = 0;
ChipWall_Separation_Yneg = 1.5e-3 - Chip_YSize;

Solution_XSize = Chip_XSize + ChipWall_Separation_Xpos + ChipWall_Separation_Xneg;
Solution_YSize = Chip_YSize + ChipWall_Separation_Ypos + ChipWall_Separation_Yneg;

DNACluster_Velocity = 1.73e-6;

Reaction_TimeLength = 2400;

T_Reaction = 336.15; % Kelvin

pH = 7; % Initial pH of solution

D_H = 1e-8; % Diffusivity of proton in water

D_HDNA = 0.25e-8; % Diffusivity of proton in DNA solution

Burst = 4e12; % Constant of proportionality between cluster surface area and proton production

H_Max = 2e-6; % Maximum mesh element size

H_Min = 2e-7; % Minimum mesh element size

Tolerance = 1E-3;

RIParameter = 2*DNACluster_Velocity; % Reinitialisation Parameter

% REQUIRED FOR DNA TRAPPING (ElectrodeActivation_Flag = 1)

% Determine starting position of DNA clusters with added trapping.
% Variables required based on flag setting for Trapping_Flag
TrappingRegion_Height = 1e-3; % Trapping_Flag = 2
TrappingRegion_ChipSeparation = 0.5e-4; %separation between the chip and the trapping region in the x and y directions

Trapping_Percentage = 50; %Trapping_Flag = 1

Trapping_File = 'TrappingPositions.txt'; %Trapping_Flag = 3

Electrode_Thickness = 0;
Castellation_NPerPixel = 0;
Castellation_Length = 0;
Castellation_Width = 0;
Castellation_Separation = 0;

% REQUIRED FOR SURFACE REACTIONS (Simulation_Flag > 0)

% Equilibrium constants for the protonation of surface states
pK_SiO = 7.5;
pK_SiOH = -2.9;

% Density of surface sites
SurfaceSiteDensity = 5e18;

%% END OF USER INPUT

% Create random co-ordinates for DNA clusters
DNACluster_CentreX = ClusterWall_Tolerance + (Solution_XSize-2*ClusterWall_Tolerance)*rand(1,DNACluster_N);
DNACluster_CentreY = ClusterWall_Tolerance + (Solution_YSize-2*ClusterWall_Tolerance)*rand(1,DNACluster_N);
DNACluster_CentreZ = ClusterWall_Tolerance + (SolutionHeight-2*ClusterWall_Tolerance)*rand(1,DNACluster_N);

%If electrodes are not activated, the random coordinates are all kept
%and returned to the main program
if ElectrodeActivation_Flag == 1
    
    if Trapping_Flag < 3

        Trapped_Index = 1:DNACluster_N;

        %Percentile with smallest z trapped
        if Trapping_Flag == 1

            %Set z coord of lowest percentile to ISFET surface
            Trapped_Index = find(DNACluster_CentreZ < prctile( DNACluster_CentreZ, Trapping_Percentage));

        elseif Trapping_Flag == 2

            %Find x and y bounds of trapping region
            TrappingRegion_XStart = TrappingRegion_ChipSeparation;
            TrappingRegion_XEnd = Chip_XSize + TrappingRegion_ChipSeparation;

            TrappingRegion_YStart = TrappingRegion_ChipSeparation;
            TrappingRegion_YEnd = Chip_YSize + TrappingRegion_ChipSeparation;

            Trapped_Index = find( DNACluster_CentreX >= TrappingRegion_XStart & DNACluster_CentreX <= TrappingRegion_XEnd & DNACluster_CentreY >= TrappingRegion_YStart & DNACluster_CentreY <= TrappingRegion_YEnd & DNACluster_CentreZ <= TrappingRegion_Height);
        end

        %If all DNA trapped (trapping_flag = 0) 
        %place DNA randomly along the edges of the line electrodes, to which they would be attracted
        if Castellation_Flag == 1

            Castellation_PixelSelection = randi( Castellation_NPerPixel, [1, length(Trapped_Index)]) - Castellation;
            Castellation_RowSelection = randi( ISFET_Ny, [1, length(Trapped_Index)]);
            Castellation_ColumnSelection = randi ( Electrode_N, [1, length(Trapped_Index)]);
            Castellation_SideSelection = randi( 2, [1, length(Trapped_Index)]);

            Castellation_SideSelection( (Castellation_SideSelection == 1 & Castellation_ColumnSelection == 1)) = 2;
            Castellation_SideSelection( (Castellation_SideSelection == 2 & Castellation_ColumnSelection == ISFET_Nx)) = 1;
            Castellation_SideSelection( Castellation_SideSelection == 1) = -1;
            Castellation_SideSelection( Castellation_SideSelection == 2) = 1;

            %Set x and y coords of DNA cluster centre from castellation
            %trapping region
            DNACluster_CentreX( Trapped_Index) = (ElectrodeElectrode_Separation + Electrode_Width)*(Castellation_ColumnSelection - 1) + 0.5*Electrode_Width + Castellation_SideSelection*(0.5*Electrode_Width + Castellation_Length);
            DNACluster_CentreY( Trapped_Index) = ISFET_StartDisplacement + ISFET_Length*0.5 + (Castellation_RowSelection - 1)*(ISFET_Length + ISFETISFET_Separation) - (Castellation_NPerPixel*Castellation_Width + (Castellation_NPerPixel-1)*Castellation_Separation)/2 + Castellation_PixelSelection*Castellation_Width + (Castellation_PixelSelection - 1)*Castellation_Width;

            %Clusters start on top of electrode
            DNACluster_CentreZ( Trapped_Index) = Electrode_Thickness;

        elseif Castellation_Flag == 0

                %Select a random electrode for cluster origin. Then randomly select
                %one of the two edges of the electrode for cluster origin.
                Electrode_Origin = randi(Electrode_N, [1,length(Trapped_Index)]) - 1;
                Edge_Origin = randi(2, [1,length(Trapped_Index)]) - 1;

                %Compute x position. Chip layout goes in repeating units of
                %electrode_width, electrodeisfet_separation, isfet_width,
                %electrodeisfet_separation
                DNACluster_CentreX(Trapped_Index) = Electrode_Origin*(Electrode_Width + ElectrodeElectrode_Separation) + Edge_Origin*Electrode_Width;

                %Clusters can take any y position
                DNACluster_CentreY(Trapped_Index) = (ISFET_Ny*ISFET_Length + (ISFET_Ny - 1)*ISFETISFET_Separation + ISFET_StartSeparation + ISFET_EndSeparation)*rand(1,length(Trapped_Index));

                %Clusters start on top of electrode
                DNACluster_CentreZ(Trapped_Index) = Electrode_Thickness;
        end
        
    else
        
        %Read in csv file that has x, y, z position of cluster centre along
        %columns, and the different clusters down rows
        DNACluster_CentreXYZ = csvread(Trapping_File);
        
        DNACluster_CentreX = DNACluster_CentreXYZ(:,1);
        DNACluster_CentreY = DNACluster_CentreXYZ(:,2);
        DNACluster_CentreZ = DNACluster_CentreXYZ(:,3);
    end
end

% Remove DNA that will not amplify before being overlapped by another
% cluster

% Save index of overlapped DNA
CoordsDelete = [];

for i = 1:DNACluster_N
    for j = 1:DNACluster_N
        
        % Only interested in comparing clusters j initiating before
        % cluster i
        if j < i
            
            % Find radius of cluster j at initiation time of cluster i
            radius_at_initiation = (sigmoid_start(i) - sigmoid_start(j))*DNACluster_Velocity + DNACluster_StartRadius;
            
            % Determine distance from cluster i to cluster j
            cluster_separation = sqrt((DNACluster_CentreX(i) - DNACluster_CentreX(j))^2 + (DNACluster_CentreY(i) - DNACluster_CentreY(j))^2 + (DNACluster_CentreX(i) - DNACluster_CentreZ(j))^2);
            
            % If cluster i is within cluster j at initiation, remove
            % cluster i from the simulation
            if cluster_separation < radius_at_initiation
                
               CoordsDelete(end+1) = i;
            end
        end
    end
end

% Delete overlapped clusters
CoordsDelete = unique(CoordsDelete);
DNACluster_CentreX(CoordsDelete) = [];
DNACluster_CentreY(CoordsDelete) = [];
DNACluster_CentreZ(CoordsDelete) = [];

% Save new number of 'active' DNA clusters
DNACluster_N = DNACluster_N - length(CoordsDelete);

% Open new COMSOL model

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\lk3618\OneDrive - Imperial College London\Sims\COMSOL\DNAISFET')

model.component.create('comp1', true)

% Enter parameter values
model.param.set('pK_SiO', sprintf('%d', pK_SiO))
model.param.descr('pK_SiO', 'Equilibrium constant of SiO protonation.')
model.param.set('pK_SiOH', sprintf('%d', pK_SiOH))
model.param.descr('pK_SiOH', 'Equilibrium constant of SiOH protonation')
model.param.set('ISFET_Nx', sprintf('%d', ISFET_Nx))
model.param.descr('ISFET_Nx', 'Number of ISFETs in X direction')
model.param.set('ISFET_Ny', sprintf('%d', ISFET_Ny))
model.param.descr('ISFET_Ny', 'Number of ISFETs in Y direction')
model.param.set('ISFET_Width', sprintf('%d', ISFET_Width))
model.param.descr('ISFET_Width', 'Width of ISFET')
model.param.set('ISFET_Length', sprintf('%d', ISFET_Length))
model.param.descr('ISFET_Length', 'Length of ISFET')
model.param.set('ISFETISFET_Separation', sprintf('%d', ISFETISFET_Separation))
model.param.descr('ISFETISFET_Separation', 'Separation in the Y direction between adjacent ISFETs')
model.param.set('ISFETElectrode_Separation', sprintf('%d', ISFETElectrode_Separation))
model.param.descr('ISFETElectrode_Separation', 'Separation in the X direction between electrodes and neighbouring ISFETs. If Electrode_Width = 0, then 2*ISFETElectrode_Separation will equal the separation between neighbouring ISFETs in the X direction.')
model.param.set('Electrode_Width', sprintf('%d', Electrode_Width))
model.param.descr('Electrode_Width', 'Width of electrodes')
model.param.set('Electrode_N', sprintf('%d', Electrode_N))
model.param.descr('Electrode_N', 'Number of electrodes. Only in X direction as line electrodes run along Y axis.')
model.param.set('ISFET_StartSeparation', sprintf('%d', ISFET_StartSeparation))
model.param.descr('ISFET_StartSeparation', 'Separation between first ISFET in Y direction and the reaction chamber wall in the plane Y = 0.')
model.param.set('ISFET_EndSeparation', sprintf('%d', ISFET_EndSeparation))
model.param.descr('ISFET_EndSeparation', 'Separation between last ISFET in Y direction and the reaction chamber wall in the plane Y = Chip_YSize.')
model.param.set('Chip_XSize', sprintf('%d', Chip_XSize))
model.param.descr('Chip_XSize', 'Total size of the chip in the X direction.')
model.param.set('ElectrodeElectrode_Separation', sprintf('%d', ElectrodeElectrode_Separation))
model.param.descr('ElectrodeElectrode_Separation', 'Total separation between adjacent electrodes in the X direction.')
model.param.set('Solution_XSize', sprintf('%d', Solution_XSize))
model.param.descr('Solution_XSize', 'Total size of the reaction chamber in the X direction.')
model.param.set('Solution_YSize', sprintf('%d', Solution_YSize))
model.param.descr('Solution_YSize', 'Total size of the reaction chamber in the Y direction.')
model.param.set('Chip_YSize', sprintf('%d', Chip_YSize))
model.param.descr('Chip_YSize', 'Total size of the chip in the Y direction.')
model.param.set('SolutionHeight', sprintf('%d', SolutionHeight))
model.param.descr('SolutionHeight', 'Size of the reaction chamber in the z direction.')
model.param.set('DNACluster_StartRadius', sprintf('%d', DNACluster_StartRadius))
model.param.descr('DNACluster_StartRadius', 'Starting radius of the DNA clusters.')
model.param.set('DNACluster_Velocity', sprintf('%d', DNACluster_Velocity))
model.param.descr('DNACluster_Velocity', 'Velocity with which the radius of the DNA clusters increases.')
model.param.set('RelativeTime_Adjustment', sprintf('%d', sigmoid_start(1)))
model.param.descr('RelativeTime_Adjustment', 'Time until first amplification event at t=0')
model.param.set('pH', sprintf('%d', pH))
model.param.descr('pH', 'Starting pH of the solution.')
model.param.set('T_Reaction', sprintf('%d', T_Reaction))
model.param.descr('T_Reaction', 'Constant temperature of the reaction.')
model.param.set('D_H', sprintf('%d', D_H))
model.param.descr('D_H', 'Diffusivity of H+ in water.')
model.param.set('Burst', sprintf('%d', Burst));
model.param.descr('Burst', 'Constant of proportionality between the surface area increase of DNA clusters and the released H+. Previous value of 2e14 proved ~50 times too large.');
model.param.set('Time_Step', '5')
model.param.descr('Time_Step', 'Time step of time dependent solver, at which solutions are stored.')
model.param.set('Time_End', sprintf('%d', Reaction_TimeLength - sigmoid_start(1)))
model.param.descr('Time_End', 'End time for the reaction.')
model.param.set('N_A', '6.02E23[mol^(-1)]')
model.param.descr('N_A', 'Avagadro Number')
model.param.set('SurfaceSiteConcentration', sprintf('%d/N_A', SurfaceSiteDensity))
model.param.descr('SurfaceSiteConcentration', 'Concentration of surface sites.')
model.param.set('M_H', '1.01E-3')
model.param.descr('M_H', 'Molar mass of H+.')
model.param.set('M_SiO', '0.0441')
model.param.descr('M_SiO', 'Molar mass of SiO.')
model.param.set('M_SiOH', '0.0451')
model.param.descr('M_SiOH', 'Molar mass of SiOH')
model.param.set('M_SiOH2', '0.0461')
model.param.descr('M_SiOH2', 'Molar mass of SiOH2')
model.param.set('D_HDNA', sprintf('%d', D_HDNA))
model.param.descr('D_HDNA', 'Diffusivity of H+ in DNA.')
model.param.set('H_Max', sprintf('%d', H_Max))
model.param.descr('H_Max', 'Maximum mesh element size')
model.param.set('H_Min', sprintf('%d', H_Min))
model.param.descr('H_Min', 'Minimum mesh element size')
model.param.set('ChipWall_Separation_Xpos', sprintf('%d', ChipWall_Separation_Xpos))
model.param.descr('ChipWall_Separation_Xpos', 'Separation between the edge of the chip and the reaction chamber wall in the positive x direction.')
model.param.set('ChipWall_Separation_Xneg', sprintf('%d', ChipWall_Separation_Xneg))
model.param.descr('ChipWall_Separation_Xneg', 'Separation between the edge of the chip and the reaction chamber wall in the negative x direction.')
model.param.set('ChipWall_Separation_Ypos', sprintf('%d', ChipWall_Separation_Ypos))
model.param.descr('ChipWall_Separation_Ypos', 'Separation between the edge of the chip and the reaction chamber wall in the positive y direction.')
model.param.set('ChipWall_Separation_Yneg', sprintf('%d', ChipWall_Separation_Yneg))
model.param.descr('ChipWall_Separation_Yneg', 'Separation between the edge of the chip and the reaction chamber wall in the negative y direction.')
model.param.set('RIParameter', sprintf('%d', RIParameter))
model.param.descr('RIParameter', 'Reinitialisation parameter for the level set method.')

%Create starting time parameter for number of clusters defined by user
for i = 1:DNACluster_N
    model.param.set(sprintf('DNACluster%d_Start', int16(i)), sprintf('%d', sigmoid_start_relative(i)))
    model.param.descr(sprintf('DNACluster%d_Start', int16(i)), sprintf('Start time for DNA cluster %d', int16(i)))
    model.param.set(sprintf('DNACluster%d_CentreX', int16(i)), sprintf('%d', DNACluster_CentreX(i)))
    model.param.descr(sprintf('DNACluster%d_CentreX', int16(i)), sprintf('X co-ordinate of centre of DNA cluster %d', int16(i)))
    model.param.set(sprintf('DNACluster%d_CentreY', int16(i)), sprintf('%d', DNACluster_CentreY(i)))
    model.param.descr(sprintf('DNACluster%d_CentreY', int16(i)), sprintf('Y co-ordinate of centre of DNA cluster %d', int16(i)))
    model.param.set(sprintf('DNACluster%d_CentreZ', int16(i)), sprintf('%d', DNACluster_CentreZ(i)))
    model.param.descr(sprintf('DNACluster%d_CentreZ', int16(i)), sprintf('Z co-ordinate of centre of DNA cluster %d', int16(i)))
end

% Create new geometry
model.component('comp1').geom.create('geom1', 3)

%Create block for reaction solution
model.component('comp1').geom('geom1').create('blk1', 'Block')
model.component('comp1').geom('geom1').feature('blk1').set('size', {'Solution_XSize' 'Solution_YSize' 'SolutionHeight'})
model.component('comp1').geom('geom1').feature('blk1').label('ReactionChamber')

%Create array of rectangles to signify ISFET array
model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('unite', true);
model.component('comp1').geom('geom1').feature('wp1').set('selresult', true);

%Create each reactangle individually so that it is easier to arrange
%results later
for i = 1:ISFET_Nx
    for j = 1:ISFET_Ny
        
        model.component('comp1').geom('geom1').feature('wp1').geom.create(sprintf('r%d', (j-1)*ISFET_Nx + i), 'Rectangle');
        model.component('comp1').geom('geom1').feature('wp1').geom.feature(sprintf('r%d', (j-1)*ISFET_Nx + i)).set('size', {'ISFET_Width' '1'});
        model.component('comp1').geom('geom1').feature('wp1').geom.feature(sprintf('r%d', (j-1)*ISFET_Nx + i)).setIndex('size', 'ISFET_Length', 1);
        model.component('comp1').geom('geom1').feature('wp1').geom.feature(sprintf('r%d', (j-1)*ISFET_Nx + i)).set('pos', {sprintf('ChipWall_Separation_Xneg + Electrode_Width + ISFETElectrode_Separation + %d*(ISFET_Width + 2*ISFETElectrode_Separation + Electrode_Width)', i-1) '0'});
        model.component('comp1').geom('geom1').feature('wp1').geom.feature(sprintf('r%d', (j-1)*ISFET_Nx + i)).setIndex('pos', sprintf('ChipWall_Separation_Yneg + ISFET_StartSeparation + %d*(ISFET_Length + ISFETISFET_Separation)', j-1), 1); 
        model.component('comp1').geom('geom1').feature('wp1').geom.feature(sprintf('r%d', (j-1)*ISFET_Nx + i)).set('selresultshow', 'all')
    end
end

%Create cumulative selection for clusters
model.component('comp1').geom('geom1').selection.create('csel1', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('csel1').label('Cumulative Selection 1');

%Store variables for reaction rates
model.component('comp1').variable.create('var1');
    
for i = 1:DNACluster_N
    
    % Create each DNA cluster in the geometry
    model.component('comp1').geom('geom1').create(sprintf('sph%d', int16(i)), 'Sphere')
    model.component('comp1').geom('geom1').feature(sprintf('sph%d', int16(i))).label(sprintf('DNACluster%d', int16(i)))
    model.component('comp1').geom('geom1').feature(sprintf('sph%d', int16(i))).set('pos', {sprintf('DNACluster%d_CentreX', int16(i)) sprintf('DNACluster%d_CentreY', int16(i)) sprintf('DNACluster%d_CentreZ', int16(i))})
    model.component('comp1').geom('geom1').feature(sprintf('sph%d', int16(i))).set('r', 'DNACluster_StartRadius')
    
    % Set selection level to all
    model.component('comp1').geom('geom1').feature(sprintf('sph%d', int16(i))).set('selresult', true)
    model.component('comp1').geom('geom1').feature(sprintf('sph%d', int16(i))).set('selresultshow', 'all')
    
    % Create step functions for each DNA cluster transitioning at the
    % cluster's start time for amplification
    model.component('comp1').func.create(sprintf('step%d', int16(i)), 'Step')
    model.component('comp1').func(sprintf('step%d', int16(i))).label(sprintf('DNACluster%d_Step', int16(i)))
    model.component('comp1').func(sprintf('step%d', int16(i))).set('location', sprintf('DNACluster%d_Start', int16(i)))
  
    %Add cluster to selection
    model.component('comp1').geom('geom1').feature(sprintf('sph%d', int16(i))).set('contributeto', 'csel1');
    
end

%Subtract DNA clusters from reaction solution
model.component('comp1').geom('geom1').create('dif1', 'Difference')
model.component('comp1').geom('geom1').feature('dif1').selection('input').set('blk1')
model.component('comp1').geom('geom1').feature('dif1').selection('input2').named('csel1')
model.component('comp1').geom('geom1').feature('dif1').set('selresult', true);

% Set model transparency so the reaction chamber can be seen through
model.component('comp1').view('view1').set('transparency', true)

% Build the geometry
model.geom('geom1').run

% Level Set Physics
model.component('comp1').physics.create('ls', 'LevelSet', 'geom1')
model.component('comp1').physics('ls').feature('lsm1').set('gamma', 'RIParameter')
model.component('comp1').physics('ls').feature('lsm1').set('u', {'-ls.intnormx*DNACluster_Velocity' '-ls.intnormy*DNACluster_Velocity' '-ls.intnormz*DNACluster_Velocity'})

% Create inlets for each DNA cluster
for i = 1:DNACluster_N
    
    model.component('comp1').physics('ls').create(sprintf('inl%d', int16(i)), 'InletBoundary', 2)
    model.component('comp1').physics('ls').feature(sprintf('inl%d', int16(i))).selection.named(sprintf('geom1_sph%d_bnd', int16(i)))
    model.component('comp1').physics('ls').feature(sprintf('inl%d', int16(i))).set('lscond', 'SpecifyLevelSetFunctionExplicitly')
    model.component('comp1').physics('ls').feature(sprintf('inl%d', int16(i))).set('ls0', sprintf('step%d(t[1/s])', int16(i)))
end

% Transport of Diluted Species Physics
model.component('comp1').physics.create('tds', 'DilutedSpecies', {'c_H'})
model.component('comp1').physics('tds').create('reac1', 'Reactions', 3)
model.component('comp1').physics('tds').feature('init1').setIndex('initc', '10^(-pH)*1e3', 0) %Note conversion from mol/L to mol/m^3
model.component('comp1').physics('tds').feature('cdm1').set('D_c_H', {'D_H*ls.Vf1+ D_HDNA*ls.Vf2' '0' '0' '0' 'D_H*ls.Vf1+ D_HDNA*ls.Vf2' '0' '0' '0' 'D_H*ls.Vf1+ D_HDNA*ls.Vf2'})

% Reaction rate string concatenation
model.component('comp1').physics('tds').feature('reac1').setIndex('R_c_H','ls.delta*2*Burst/N_A' , 0)
model.component('comp1').physics('tds').feature('reac1').selection.set(1);


% Add Chemistry and Surface Reactions
if Simulation_Flag > 0
    
    % Chemistry
    model.component('comp1').physics.create('chem', 'Chemistry', 'geom1');
    
    % Create species: H, SiO, SiOH, SiOH2
    model.component('comp1').physics('chem').create('spec1', 'SpeciesChem', 3);
    model.component('comp1').physics('chem').feature('spec1').set('specName', 'H');
    model.component('comp1').physics('chem').create('spec1', 'SpeciesChem', 3);
    model.component('comp1').physics('chem').feature('spec1').set('specName', 'SiO');
    model.component('comp1').physics('chem').create('spec1', 'SpeciesChem', 3);
    model.component('comp1').physics('chem').feature('spec1').set('specName', 'SiOH');
    model.component('comp1').physics('chem').create('spec1', 'SpeciesChem', 3);
    model.component('comp1').physics('chem').feature('spec1').set('specName', 'SiOH2');
    model.component('comp1').physics('chem').feature('SiO').set('sType', 'surface');
    model.component('comp1').physics('chem').feature('SiOH').set('sType', 'surface');
    model.component('comp1').physics('chem').feature('SiOH2').set('sType', 'surface');
    model.component('comp1').physics('chem').feature('SiOH2_surf').set('M', 'M_SiOH2[kg/mol]');
    model.component('comp1').physics('chem').feature('SiOH2_surf').set('z', 1);
    model.component('comp1').physics('chem').feature('SiOH_surf').set('M', 'M_SiOH[kg/mol]');
    model.component('comp1').physics('chem').feature('SiO_surf').set('M', 'M_SiO[kg/mol]');
    model.component('comp1').physics('chem').feature('SiO_surf').set('z', -1);
    model.component('comp1').physics('chem').feature('H').set('M', 'M_H[kg/mol]');
    model.component('comp1').physics('chem').feature('H').set('z', 1);
    model.component('comp1').physics('chem').prop('ChemistryModelInputParameter').setIndex('csurf', 'c_SiO', 0, 0);
    model.component('comp1').physics('chem').prop('ChemistryModelInputParameter').setIndex('csurf', 'c_SiOH', 1, 0);
    model.component('comp1').physics('chem').prop('ChemistryModelInputParameter').setIndex('csurf', 'c_SiOH2', 2, 0);
    
    % Create reactions
    model.component('comp1').physics('chem').create('rch1', 'ReactionChem', 3);
    model.component('comp1').physics('chem').feature('rch1').set('formula', 'H+SiO(ads) = SiOH2(ads)');
    model.component('comp1').physics('chem').feature('rch1').set('rType', 'rev');
    model.component('comp1').physics('chem').feature('rch1').set('setKeq0', true);
    model.component('comp1').physics('chem').feature('rch1').set('Keq0', '10^(pK_SiO)');
    model.component('comp1').physics('chem').create('rch2', 'ReactionChem', 3);
    model.component('comp1').physics('chem').feature('rch2').set('formula', 'H + SiOH(ads) = SiOH2(ads)');
    model.component('comp1').physics('chem').feature('rch2').set('rType', 'rev');
    model.component('comp1').physics('chem').feature('rch2').set('setKeq0', true);
    model.component('comp1').physics('chem').feature('rch2').set('Keq0', '10^(pK_SiOH)');
    
    % Surface Reactions
    model.component('comp1').physics.create('sr', 'SurfaceReactions', {'cs1'; ''});
    ChipSurface_Entities = mphselectbox(model, 'geom1', [(-1e-6 - ChipWall_Separation_Xneg) (-1e-6 - ChipWall_Separation_Yneg) -1e-6;  (Chip_XSize + ChipWall_Separation_Xpos + 1e-6) (Chip_YSize + ChipWall_Separation_Ypos + 1e-6) 1e-6]', 'boundary');
    model.component('comp1').physics('sr').selection.set([ChipSurface_Entities]);
    model.component('comp1').physics('sr').feature('sp1').set('Gamma', 'SurfaceSiteConcentration[mol/m^2]');
    model.component('comp1').physics('chem').prop('ChemistryModelInputParameter').setIndex('ConcentrationInput', 'c_H', 0, 0);
    model.component('comp1').physics('sr').field('surfaceconcentration').component(1, 'c_SiO');
    model.component('comp1').physics('sr').field('surfaceconcentration').component(2, 'c_SiOH');
    model.component('comp1').physics('sr').field('surfaceconcentration').component(3, 'c_SiOH2');
    model.component('comp1').physics('sr').feature('init1').setIndex('initcs', 'SurfaceSiteConcentration', 1);
    model.component('comp1').physics('sr').create('reac1', 'Reactions', 2);
    model.component('comp1').physics('sr').feature('reac1').setIndex('Rs_c_SiO_src', 'root.comp1.chem.R_SiO_surf', 0);
    model.component('comp1').physics('sr').feature('reac1').setIndex('Rs_c_SiOH_src', 'root.comp1.chem.R_SiOH_surf', 0);
    model.component('comp1').physics('sr').feature('reac1').setIndex('Rs_c_SiOH2_src', 'root.comp1.chem.R_SiOH2_surf', 0);
    model.component('comp1').physics('sr').feature('reac1').selection.set([ChipSurface_Entities]);

    % Ammend Transport of Diluted Species
    model.component('comp1').physics('tds').create('srf1', 'SurfaceReactionsFlux', 2);
    model.component('comp1').physics('tds').feature('srf1').setIndex('J0_c_H_src', 'root.comp1.chem.R_SiO_surf', 0);
    model.component('comp1').physics('tds').feature('srf1').selection.set([ChipSurface_Entities]);
    
end

%Set mesh
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').automatic(false);
model.component('comp1').mesh('mesh1').feature('size').set('custom', true);
model.component('comp1').mesh('mesh1').feature('size').set('hmin', sprintf('%d', H_Min));
model.component('comp1').mesh('mesh1').feature('size').set('hmax', sprintf('%d', H_Max));

%Create time dependent solver
model.study.create('std1');
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').activate('ls', true);
model.study('std1').feature('time').activate('tds', true);
model.study('std1').feature('time').set('tlist', 'range(0,Time_Step,Time_End)');

model.sol.create('sol1');
model.sol('sol1').study('std1');

model.study('std1').feature('time').set('notlistsolnum', 1);
model.study('std1').feature('time').set('notsolnum', '1');
model.study('std1').feature('time').set('listsolnum', 1);
model.study('std1').feature('time').set('solnum', '1');

model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'time');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,Time_Step,Time_End)');
model.sol('sol1').feature('t1').set('plot', 'off');
model.sol('sol1').feature('t1').set('plotgroup', 'Default');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('rtol', Tolerance);
model.sol('sol1').feature('t1').set('atolglobalvaluemethod', 'factor');
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').create('se1', 'Segregated');
model.sol('sol1').feature('t1').feature('se1').feature.remove('ssDef');
model.sol('sol1').feature('t1').feature('se1').create('ss1', 'SegregatedStep');
model.sol('sol1').feature('t1').feature('se1').feature('ss1').set('segvar', 'comp1_phils');
model.sol('sol1').feature('t1').feature('se1').feature('ss1').set('subdamp', 0.8);
model.sol('sol1').feature('t1').feature('se1').feature('ss1').set('subjtech', 'once');
model.sol('sol1').feature('t1').create('i1', 'Iterative');
model.sol('sol1').feature('t1').feature('i1').set('linsolver', 'gmres');
model.sol('sol1').feature('t1').feature('i1').set('prefuntype', 'left');
model.sol('sol1').feature('t1').feature('i1').set('rhob', 40);
model.sol('sol1').feature('t1').feature('i1').set('itrestart', 50);
model.sol('sol1').feature('t1').feature('i1').set('maxlinit', 400);
model.sol('sol1').feature('t1').feature('i1').set('nlinnormuse', 'on');
model.sol('sol1').feature('t1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').set('prefun', 'gmg');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').set('mcasegen', 'any');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('pr').create('sl1', 'SORLine');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('iter', 2);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.2);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('seconditer', 1);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('relax', 0.4);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('po').create('sl1', 'SORLine');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').set('iter', 2);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.2);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').set('seconditer', 2);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('po').feature('sl1').set('relax', 0.4);
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('t1').feature('se1').feature('ss1').set('linsolver', 'i1');
model.sol('sol1').feature('t1').feature('se1').feature('ss1').label('Level Set');
model.sol('sol1').feature('t1').feature('se1').create('ss2', 'SegregatedStep');
%model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('segvar', 'comp1_c');
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('subdamp', 0.8);
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('subjtech', 'once');
model.sol('sol1').feature('t1').create('i2', 'Iterative');
model.sol('sol1').feature('t1').feature('i2').set('linsolver', 'gmres');
model.sol('sol1').feature('t1').feature('i2').set('prefuntype', 'left');
model.sol('sol1').feature('t1').feature('i2').set('rhob', 40);
model.sol('sol1').feature('t1').feature('i2').set('itrestart', 50);
model.sol('sol1').feature('t1').feature('i2').set('maxlinit', 400);
model.sol('sol1').feature('t1').feature('i2').set('nlinnormuse', 'on');
model.sol('sol1').feature('t1').feature('i2').create('mg1', 'Multigrid');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').set('prefun', 'gmg');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').set('mcasegen', 'any');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('pr').create('sl1', 'SORLine');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('iter', 2);
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.2);
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('seconditer', 1);
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('pr').feature('sl1').set('relax', 0.4);
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('po').create('sl1', 'SORLine');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('po').feature('sl1').set('iter', 2);
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.2);
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('po').feature('sl1').set('seconditer', 2);
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('po').feature('sl1').set('relax', 0.4);
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('t1').feature('i2').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('linsolver', 'i2');
model.sol('sol1').feature('t1').feature('se1').feature('ss2').label('Concentration c');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').feature('t1').set('timemethod', 'genalpha');
model.sol('sol1').feature('t1').set('estrat', 'exclude');
model.study('std1').feature('time').set('plot', false);
model.study('std1').feature('time').set('probesel', 'none');

mphlaunch

out = model