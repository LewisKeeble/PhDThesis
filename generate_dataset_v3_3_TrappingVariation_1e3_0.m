%% SimulatingLAMPKunal
% Moniri + Lewis with changes
% Castellations not considered
% For simplicity protons can pass through DNA clusters
% No Electrodes
% Trapping regions are used as before
% v3.3 removes clusters that don't initiate before being enveloped by other
% clusters to reduce runtime

function generate_dataset_v3_3_TrappingVariation_1e3_0(run)

    rng(run)

    %% Reaction
    num_of_mol = 1e3; % # of initial molecules
    speed_of_front = 1.86E-6 ; % metres per second
    r_initial = 1.25E-6; % in metres

    % Starting time to each cluster is Normal distribution. Values taken from
    % fit to real data
    mu_start_time = 13; % in minutes
    sigma_start_time = 3; % in minutes

    % Number of protons released is proportional to area of propagation front.
    % 'burst' is the proportionality constant. Can use estimateburst script to
    % estimate the parameter based on experimental observations.
    burst = 1.8e13;
    BurstFactor = 0.01;

    % Create starting time dist (and sample from it)
    pd = makedist('Normal', 'mu', mu_start_time, 'sigma', sigma_start_time);
    sigmoid_start = round(60 * sort(random(pd, 1, num_of_mol)) / 1); % notice conversion to index 

    total_time = 40; % time in minutes
    time_step = 1; % in seconds
    t = 0:time_step:total_time*60; % in seconds

    trapping_percentage = 0;
    %% Diffusion

    %Effective 'rate constant' for H ion diffusion, k_diffusion.
    %Value from http://omh.umeche.maine.edu/pdfs/JChemPhys_135_124505.01pdf.pdf
    %for 50 degrees C. Units m^2/s.
    diffusivity_H = 1E-8;

    %% Chip design
    
    %Design parameters
    %Line electrodes have width in x and length in y. Length is equal to the
    %sensor_ysize. For no electrodes, best to put electrode_width as zero and halve the
    %electrodeisfet_separation
    chip.N_x = 100;
    chip.N_y = 100; %N_x x N_y array
    chip.electrode_width = 52.1E-6;
    chip.electrode_thickness = 0.05E-6;
    chip.N_electrodes = chip.N_x+1; %Number of electrodes on chip

    chip.electrodeisfet_separation = (95.9E-6-2*35e-6-2e-6)/2; %Separation between the outermost point of electrodes and the edges of ISFETs along x axis
    chip.isfetisfet_separation = 2E-6; %Separation between the edges of neighbouring ISFETs in the y axis
    chip.isfet_width = 35E-6;
    chip.isfet_length = 35E-6;

    %Distance between first and last ISFETs and edge of chip in y-axis
    chip.sensor_startdisplacement = 0E-6;
    chip.sensor_enddisplacement = 0E-6;

    %Width of individual castellations in the y-axis. Extension in x-axis
    chip.castellation_width = 0;
    chip.castellation_extension = 0;

    %If only 1 castellation per pixel, castellation separation will be set to 0
    chip.ncastellation_perpixel = 0;
    chip.castellation_separation = 0;

    %Separation between the edge of the chip and the wall of the reaction
    %chamber
    chip.wall_separation_xpos = 0;%(1.75e-3 + 2*0.262e-3 - (chip.N_x*(2*chip.electrodeisfet_separation+chip.isfet_width) + chip.N_electrodes*chip.electrode_width))/2;
    chip.wall_separation_xneg = 0;%chip.wall_separation_xpos;
    chip.wall_separation_ypos = 0;
    chip.wall_separation_yneg = 0;%1.5e-3 - (chip.sensor_startdisplacement + chip.N_y*(chip.isfet_length) + (chip.N_y-1)*(chip.isfetisfet_separation) + chip.sensor_enddisplacement);

    %Calculate the total x and y size of the sensor array
    sensor_xsize = chip.N_x*(2*chip.electrodeisfet_separation+chip.isfet_width) + chip.N_electrodes*chip.electrode_width;
    sensor_ysize = chip.sensor_startdisplacement + chip.N_y*(chip.isfet_length) + (chip.N_y-1)*(chip.isfetisfet_separation) + chip.sensor_enddisplacement;

    reaction_xsize = 3.1397e-3;%sensor_xsize + chip.wall_separation_xpos + chip.wall_separation_xneg;
    reaction_ysize = 2e-3;%sensor_ysize + chip.wall_separation_ypos + chip.wall_separation_yneg;
    solution_height = 2.5e-3;
    
    %Calculate the total separation between adjacent electrodes
    chip.electrode_separation = 2*chip.electrodeisfet_separation + chip.isfet_width + 2*chip.castellation_extension;

    %% Sensor Locations
    
    image_resolution = 1e-6;
    
    %Coordinates for rgb image
    image_x = image_resolution/2:image_resolution:reaction_xsize-image_resolution/2;
    image_y = image_resolution/2:image_resolution:reaction_ysize-image_resolution/2;

    %Construct the chip given the user input parameters, creating a RGB image -
    %chip_layout - as well as giving the position of the sensors and electrodes
    [~, column_location, row_location, ~] = BuildChip_QuickDiffusionv3(image_x, image_y, chip, image_resolution);

    %% End of Inputs

    %% Initialise Clusters

    %If DNA exists in solution, initialise starting positions of DNA clusters
    if num_of_mol ~= 0
        clusters = InitialiseClustersRandom(num_of_mol, reaction_xsize, reaction_ysize, solution_height, r_initial);
        %clusters = InitializeClusterHeight(num_of_mol, sensor_xsize, sensor_ysize, solution_height, r_initial, 0.9);
    end

    %Set DNA with lowest percentile of z-coords to locations at electrode edges
    trapped_index = find(clusters.centre_z < prctile( clusters.centre_z, trapping_percentage));
    
    %Select a random electrode for cluster origin. Then randomly select
    %one of the two edges of the electrode for cluster origin.
    electrodes_inwell = floor(reaction_xsize/(chip.electrode_width + chip.electrode_separation));
    electrode_origin = randi(electrodes_inwell, [1,length(trapped_index)]) - 1;
    edge_origin = randi(2, [1,length(trapped_index)]) - 1;

    %Compute x position. Chip layout goes in repeating units of
    %electrode_width, electrodeisfet_separation, isfet_width,
    %electrodeisfet_separation
    clusters.centre_x(trapped_index) = electrode_origin*(chip.electrode_width + chip.electrode_separation) + edge_origin*chip.electrode_width;

    %Clusters can take any y position
    clusters.centre_y(trapped_index) = reaction_ysize*rand(1,length(trapped_index));

    %Clusters start on top of electrode
    clusters.centre_z(trapped_index) = chip.electrode_thickness;

    %{
    %If desired, create a scatter plot to show the distribution of DNA cluster
    %centres and the chip layout, as a good reality check
    PlotClustersInitial(num_of_mol, clusters, reaction_xsize, reaction_ysize, solution_height, chip_layout)
    xlim([0,reaction_xsize]);ylim([0,reaction_ysize]);zlim([0,solution_height]);
    %}             
    %Set start time of simulation. signmoid_start(1) to start at the first
    %amplification. Enter 1 to start at t=0.
    t_start = sigmoid_start(1);
    
    % Remove clusters that don't amplify before being enveloped by other
    % clusters

    i=0;
    while i<num_of_mol
        
        i=i+1; 
        j=i;
        while j<num_of_mol
            
            j=j+1;
            
            % Find radius of cluster i at initiation time of cluster j
            radius_at_initiation = (sigmoid_start(j) - sigmoid_start(i))*speed_of_front + r_initial;

            % Determine distance from cluster i to cluster j
            cluster_separation = sqrt((clusters.centre_x(i) - clusters.centre_x(j))^2 + (clusters.centre_y(i) - clusters.centre_y(j))^2 + (clusters.centre_z(i) - clusters.centre_z(j))^2);

            % If cluster j is within cluster i at initiation, remove
            % cluster j from the simulation
            if cluster_separation < radius_at_initiation

                clusters.centre_x(j) = [];
                clusters.centre_y(j) = [];
                clusters.centre_z(j) = [];
                sigmoid_start(j) = [];
                clusters.radius_prev(j) = [];
                clusters.radius(j) = [];
                % Save new number of 'active' DNA clusters
                num_of_mol = num_of_mol - 1;
                j=j-1;
            end
        end
    end

    %% Initialise Sensors
    %Array to store the number of protons at each sensor over time. col =
    %sensor, row = time
    sensor_nH = zeros(chip.N_x*chip.N_y, length(t));
    protons_released = zeros(1, length(t));

    %% Simulation

    for i = t_start:length(t)-1
        %disp(i)
        protons = [];
        
        %Carry out next step of LAMP reaction
        num_of_active_sig = sum(i>sigmoid_start);

        for j = 1:num_of_active_sig

            % If cluster is active then update propogation front
            clusters.radius_prev(j) = clusters.radius(j);
            clusters.radius(j) = speed_of_front*(t(i)-t(sigmoid_start(j))); 

            % Compute position of randomly distributed protons around propogation front
            numtorelease = round(burst*(volume_of_sphere(clusters.radius(j))  - volume_of_sphere(clusters.radius_prev(j)))*BurstFactor/speed_of_front);
            %numtorelease = round(burst*(clusters.radius(j)^2));

            c = 2*rand(numtorelease,1)-1;
            lon=2*pi*rand(numtorelease,1);
            lat=acos(c);
            a=cos(lon).*sin(lat); %random points on a sphere
            b=sin(lon).*sin(lat);

            % These are the produced protons
            additions = [clusters.radius(j)*a'+clusters.centre_x(j);clusters.radius(j)*b'+clusters.centre_y(j); clusters.radius(j)*c'+clusters.centre_z(j)];

            % Remove the computed proton positions outside the sensor
            additions(:,additions(1,:)>reaction_xsize) = []; %this should be xneg
            additions(:,additions(1,:)<0) = [];
            additions(:,additions(2,:)>reaction_ysize) = [];
            additions(:,additions(2,:)<0) = [];
            additions(:,additions(3,:)>solution_height) = [];
            additions(:,additions(3,:)<0) = [];

            % Remove protons that fall within other clusters
            for k = 1:num_of_active_sig
                if(k~=j)
                    additions(:,(additions(1,:) - clusters.centre_x(k)).^2 ...
                        + (additions(2,:) - clusters.centre_y(k)).^2 ...
                        + (additions(3,:) - clusters.centre_z(k)).^2 ...
                        <= clusters.radius(k)^2) = [];
                end
            end

            protons = [protons additions];
        end

        if size(protons, 2) > 0

            protons_released(i) = protons_released(i-1) + size(protons, 2);

            % Diffusion Equation: Z(t+del_t) = Z(t) + sqrt(2.D.del_t).N(0,1)
            % Rearrange to get del_t, time taken to reach detection height
            del_t = (erfinv(1-rand(1, size(protons,2)))*2*sqrt(diffusivity_H)./(protons(3,:))).^(-2);
            %del_t = (((detection_distance*ones(1, size(protons,2)) - protons(3,:))./normrnd(0,1, [1, size(protons,2)])).^2)./(2*diffusivity_H);

            % Use X, Y diffusion equations to compute final locations of protons
            % protons_final has 4 rows: x,y,z,t
            protons_final = zeros(size(protons, 1), size(protons, 2));

            % Calculate upper and lower limits of the normal function CDF for the x and y dimensions 
            CDF_xLL = 0.5 + 0.5*erf(-protons(1,:)./(2*sqrt(diffusivity_H*del_t)));
            CDF_xUL = 0.5 + 0.5*erf((reaction_xsize - protons(1,:))./(2*sqrt(diffusivity_H*del_t)));

            CDF_yLL = 0.5 + 0.5*erf(-protons(2,:)./(2*sqrt(diffusivity_H*del_t)));
            CDF_yUL = 0.5 + 0.5*erf((reaction_ysize - protons(2,:))./(2*sqrt(diffusivity_H*del_t)));

            % Determine diffusion displacement in x and y dimensions given the
            % time to reach chip surface, del_t
            protons_final(1,:) = protons(1,:) + erfinv(2.*(rand(1, size(protons,2)).*(CDF_xUL-CDF_xLL)+CDF_xLL) - 1).*2.*sqrt(diffusivity_H.*del_t);
            protons_final(2,:) = protons(2,:) + erfinv(2.*(rand(1, size(protons,2)).*(CDF_yUL-CDF_yLL)+CDF_yLL) - 1).*2.*sqrt(diffusivity_H.*del_t);

            % Determine diffusion displacement in x and y dimensions given the
            % time to reach chip surface, del_t
            %protons_final(1,:) = protons(1,:) + erfinv(2.*rand(1, size(protons,2)) - 1).*2.*sqrt(diffusivity_H.*del_t);
            %protons_final(2,:) = protons(2,:) + erfinv(2.*rand(1, size(protons,2)) - 1).*2.*sqrt(diffusivity_H.*del_t);

            protons_final(3,:) = zeros(1, size(protons,2));
            protons_final(4,:) = ceil((del_t+t(i))./time_step);
    %         disp('protons_final before removal:')
    %         disp(protons_final)

            % Remove protons outside set reaction time
            protons_final(:, (protons_final(4,:)>length(t))) = [];

            % Remove protons outside sensor array
            protons_final(:,protons_final(1,:)>=reaction_xsize - chip.wall_separation_xpos) = [];
            protons_final(:,protons_final(1,:)<=chip.wall_separation_xneg) = [];
            protons_final(:,protons_final(2,:)>=reaction_ysize - chip.wall_separation_ypos-chip.sensor_enddisplacement) = [];
            protons_final(:,protons_final(2,:)<=chip.wall_separation_yneg + chip.sensor_startdisplacement) = [];

            % Find nearest sensor to intersection point. Find minimum,
            % positive value for difference between x and y coords of
            % intersection point with x and y starting points of ISFET
            % columns and rows
            column_intersection = column_location(:,1) - protons_final(1,:); 
            row_intersection = row_location(:,1) - protons_final(2,:);

            % Set initially to 1. If the protons come before the first
            % sensor this will be picked up in the later check and avoid
            % throwing index errors
            column_nearest = ones(1,size(protons_final,2));
            row_nearest = ones(1,size(protons_final,2));

            for j = 1:size(protons_final,2)

                column_nearest_temp = find( column_intersection(:,j) < 0, 1, 'last');
                if ~isempty(column_nearest_temp)
                    column_nearest(j) = column_nearest_temp;
                end

                row_nearest_temp = find( row_intersection(:,j) < 0, 1, 'last');
                if ~isempty(row_nearest_temp)
                    row_nearest(j) = row_nearest_temp;
                end
            end

            incolumn_flag = logical(protons_final(1,:) <= column_location( column_nearest, 2)');
            inrow_flag = logical(protons_final(2,:) <= row_location( row_nearest, 2)');

            insensor_flag = incolumn_flag & inrow_flag;

            column_nearest = column_nearest(insensor_flag);
            row_nearest = row_nearest(insensor_flag);       
            protons_final = protons_final(:,insensor_flag);

            for j = 1:size(protons_final,2)

                % Add proton to sensor
                sensor_nH((column_nearest(j) - 1)*chip.N_y + row_nearest(j), protons_final(4,j)) = sensor_nH((column_nearest(j) - 1)*chip.N_y + row_nearest(j), protons_final(4,j)) + 1;
            end

        else
            protons_released(i) = protons_released(i-1);
            %If (1) there is a positive number of protons in the solution
            %& No new protons were released, this is a signal to stop rxn
            if max(protons_released) > 0 
               break
            end
        end

    end

    save(sprintf('SensornH_Titan_Trapping_1e3Copy_BurstFactor1e1_Run%d.txt',run),'sensor_nH','-ASCII','-append')
    save(sprintf('ProtonsReleased_Titan_Trapping_1e3Copy_BurstFactor1e1_Run%d',run),'protons_released','-ASCII','-append')
end
