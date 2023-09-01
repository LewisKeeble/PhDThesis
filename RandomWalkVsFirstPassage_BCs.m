% Random walk vs First Passage with Boundary Conditions
%clear
time_step_RW = 0.01; % in seconds

total_time_RW = 2400;

%t_H = 0;

diffusivity_H = 1E-8;

protons_nH_initial = 1e5;

reaction_xsize = 0.0022;
reaction_ysize = 0.0015;
solution_height = 0.0025;

noceiling_flag = 0;
reflectingwalls_flag = 1;
reflectingceiling_flag = 0;

start_z = [0.1,0.5,0.9];
%{
for i = 1:3

    
    protons = ones(3,protons_nH_initial).*[reaction_xsize/2;reaction_ysize/2;start_z(i)*solution_height];

    % Diffusion Equation: Z(t+del_t) = Z(t) + sqrt(2.D.del_t).N(0,1)
    % Rearrange to get del_t, time taken to reach detection height
    del_t = (erfinv(1-rand(1, size(protons,2)))*2*sqrt(diffusivity_H)./(protons(3,:))).^(-2);

    % Use X, Y diffusion equations to compute final locations of protons
    % protons_final has 4 rows: x,y,z,t
    protons_final = zeros(4, size(protons, 2));

    %% Absorbing

    % Determine diffusion displacement in x and y dimensions given the
    % time to reach chip surface, del_t
    protons_final(1,:) = protons(1,:) + normrnd(0,1,[1, size(protons,2)]).*sqrt(2*diffusivity_H.*del_t);
    protons_final(2,:) = protons(2,:) + normrnd(0,1,[1, size(protons,2)]).*sqrt(2*diffusivity_H.*del_t);

    protons_final(3,:) = zeros(1, size(protons,2));
    protons_final(4,:) = del_t;

    % Remove protons outside sensor array
    protons_final(:,protons_final(1,:)>=reaction_xsize) = [];
    protons_final(:,protons_final(1,:)<=0) = [];
    protons_final(:,protons_final(2,:)>=reaction_ysize) = [];
    protons_final(:,protons_final(2,:)<=0) = [];
    
    save(sprintf('RandomWalkVsFirstPassage_BCs_FP_Abs_%dz_Test',start_z(i)), 'protons_final')
    
end
%}
for i = 2:3

    t_H = 0;
    
    protons = ones(3,protons_nH_initial).*[reaction_xsize/2;reaction_ysize/2;start_z(i)*2.5e-3];%solution_height];
    
    protons_final = [];
    
    plane_n = [[1,0,0];[1,0,0];[0,1,0];[0,1,0];[0,0,1];[0,0,1]];
    plane_p = [[0, reaction_ysize/2, solution_height/2];[reaction_xsize, reaction_ysize/2, solution_height/2];[reaction_xsize/2, 0, solution_height/2];[reaction_xsize/2, reaction_ysize, solution_height/2];[reaction_xsize/2, reaction_ysize/2, 0];[reaction_xsize/2, reaction_ysize/2, solution_height]];
    if noceiling_flag == 1
       plane_n(end,:)=[];
       plane_p(end,:)=[];
    end
    % Carry out diffusion until the next step in the LAMP reaction
    while t_H < total_time_RW

            %Save vector of diffusion movement
            diffusion_vector = [sqrt(2*diffusivity_H*time_step_RW)*normrnd(0,1,1,size(protons,2)); sqrt(2*diffusivity_H*time_step_RW)*normrnd(0,1,1,size(protons,2)); sqrt(2*diffusivity_H*time_step_RW)*normrnd(0,1,1,size(protons,2))];
            d_interesct = zeros(6,size(protons,2));
            for j = 1:6
                plane_p(j,:)'-protons
            end
            % Find intersection of path with z = detection_distance plane
            [d_intersection(1), p_intersection(1,:)] =  LinePlane_Intersection( [1,0,0], [0, reaction_ysize/2, solution_height/2], diffusion_vector_norm, protons(:,j)');
            [d_intersection(2), p_intersection(2,:)] =  LinePlane_Intersection( [1,0,0], [ reaction_xsize, reaction_ysize/2, solution_height/2], diffusion_vector_norm, protons(:,j)');
            [d_intersection(3), p_intersection(3,:)] =  LinePlane_Intersection( [0,1,0], [ reaction_xsize/2, 0, solution_height/2], diffusion_vector_norm, protons(:,j)');
            [d_intersection(4), p_intersection(4,:)] =  LinePlane_Intersection( [0,1,0], [ reaction_xsize/2, reaction_ysize, solution_height/2], diffusion_vector_norm, protons(:,j)');
            [d_intersection(5), p_intersection(5,:)] =  LinePlane_Intersection( [0,0,1], [ reaction_xsize/2, reaction_ysize/2, 0], diffusion_vector_norm, protons(:,j)'); %distance along diffusion vector to intersection. = inf if no intersection            
%{
function [d_intersection, p_intersection] =  LinePlane_Intersection_v2( plane_normalvector, plane_point, line_vector, line_point)
    d_intersection = dot((plane_point - line_point), plane_normalvector)/dot( line_vector, plane_normalvector);%distance along diffusion vector to intersection. = inf if no intersection
    p_intersection = line_point + d_intersection*line_vector; %point of intersection
end
 %}
           d_intersection(d_intersection<=0) = inf;
            
            [min_val,min_idx] = min(d_intersection);
            if min_val < diffusion_mag         
                if min_idx == 5
                    protons_final = [protons_final, [p_intersection(min_idx,:)';t_H]];
                    % Even if not within sensing region, remove proton from
                    % simulation as it is now trapped at ISFET surface
                    protons(:,j) = [];
                    % Subtract from index j to accommodate deletion of
                    % proton
                    j = j-1;
                elseif min_idx == 6
                    if reflectingceiling_flag == 0
                        % Even if not within sensing region, remove proton from
                        % simulation as it is now trapped at ISFET surface
                        protons(:,j) = [];

                        % Subtract from index j to accommodate deletion of
                        % proton
                        j = j-1;
                    else
                        protons(:,j) = p_intersection(min_idx,:)';
                    end
                else
                    if reflectingwalls_flag == 0
                        % Even if not within sensing region, remove proton from
                        % simulation as it is now trapped at ISFET surface
                        protons(:,j) = [];

                        % Subtract from index j to accommodate deletion of
                        % proton
                        j = j-1;
                    else
                        protons(:,j) = p_intersection(min_idx,:)';
                    end
                end
            else
               protons(:,j) = protons(:,j) + diffusion_vector'; 
            end
            j = j+1;
        end

        %Update time
        t_H = t_H + time_step_RW;
    end

    save(sprintf('RWVsFP_RW_Ref_AbsCeiling_%dz',start_z(i)), 'protons_final')
    
%}