% Random walk vs First Passage with Boundary Conditions

time_step_RW = 0.01; % in seconds

total_time_RW = 2400;

t_H = 0;

diffusivity_H = 1E-8;

protons_nH_initial = 1e5;

reaction_xsize = 2.2e-3;
reaction_ysize = 1.5e-3;
solution_height = 2.5e-3;

for height = [0.1]
%{
protons = ones(3,protons_nH_initial).*[0.5*reaction_xsize;0.5*reaction_ysize;height*2.5e-3];

% Diffusion Equation: Z(t+del_t) = Z(t) + sqrt(22D.del_t).N(0,1)
% Rearrange to get del_t, time taken to reach detection height
del_t = (erfinv(1-rand(1, size(protons,2)))*2*sqrt(diffusivity_H)./(protons(3,:))).^(-2);

% Use X, Y diffusion equations to compute final locations of protons
% protons_final has 4 rows: x,y,z,t
protons_final = zeros(4, size(protons, 2));

% Determine diffusion displacement in x and y dimensions given the
% time to reach chip surface, del_t
protons_final(1,:) = erfinv(2.*rand(1, size(protons,2)) - 1).*2.*sqrt(diffusivity_H.*del_t);
protons_final(2,:) = erfinv(2.*rand(1, size(protons,2)) - 1).*2.*sqrt(diffusivity_H.*del_t);
protons_final(3,:) = zeros(1, size(protons,2));
protons_final(4,:) = del_t;

save(sprintf('RandomWalkVsFirstPassage_FP_Abs_%d',height), 'protons_final')
%}
   
protons = ones(3,protons_nH_initial).*[0.5*reaction_xsize;0.5*reaction_ysize;height*2.5e-3];

protons_final = [];

% Carry out diffusion until the next step in the LAMP reaction
while t_H < total_time_RW

    % Add random diffusion in each of x, y, z to every proton
    j = 1;
    while j <= size(protons, 2)

        %Save vector of diffusion movement
        diffusion_vector = [sqrt(2*diffusivity_H*time_step_RW)*normrnd(0,1), sqrt(2*diffusivity_H*time_step_RW)*normrnd(0,1), sqrt(2*diffusivity_H*time_step_RW)*normrnd(0,1)];

        %Calculate magnitude of difusion vector
        diffusion_mag = norm(diffusion_vector);

        %Normalise diffusion vector
        diffusion_vector = diffusion_vector./diffusion_mag;

        %Save distance along diffusion vector to nearest collision
        %with DNA cluster
        collisiondistance = diffusion_mag;

        % Determine if path of proton takes it within detection region
        % of ISFET

        % Find intersection of path with z = detection_distance plane
        [d_zintersection, p_zintersection] = LinePlane_Intersection( [0,0,1], [ reaction_xsize/2, reaction_ysize/2, 0], diffusion_vector, protons(:,j)'); %distance along diffusion vector to intersection. = inf if no intersection            
        
        new_point = protons(:,j) + collisiondistance.*diffusion_vector';

        % If distance of intersection along diffusion vector is greater
        % than magnitude of diffusion, then proton doesn't reach
        % detection z plane
        if d_zintersection >= 0 && d_zintersection <= diffusion_mag
            
            collisiondistance = d_zintersection;

            protons_final = [protons_final, [p_zintersection';t_H]];

            % Even if not within sensing region, remove proton from
            % simulation as it is now trapped at ISFET surface
            protons(:,j) = [];

            % Subtract from index j to accommodate deletion of
            % proton
            j = j-1;

        elseif new_point(1) < 0 || new_point(1) > reaction_xsize || new_point(2) < 0 || new_point(2) > reaction_ysize || new_point(3) > solution_height
            
            % Even if not within sensing region, remove proton from
            % simulation as it is now trapped at ISFET surface
            protons(:,j) = [];

            % Subtract from index j to accommodate deletion of
            % proton
            j = j-1;
        else
        
            %Update xyz coords of proton
            protons(:,j) = protons(:,j) + collisiondistance.*diffusion_vector';

        end
        j = j+1;
    end

    %Update time
    t_H = t_H + time_step_RW;
end

save(sprintf('RandomWalkVsFirstPassage_RW_Abs_%d',height), 'protons_final')
%}
end