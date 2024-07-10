classdef trial_prediction_class < handle

    properties

        g = 9.806
        m

        trial_data
        % Horizon length in ms
        horizon_time

        % Sample time in ms
        sample_time_ms

        n_force_samples
        n_marker_samples

        com_pos_ref
        com_vel_ref

        com_pos_kalman
        com_vel_kalman
        com_acc_kalman
        res_f_ext_kalman

        com_pos_integration
        com_vel_integration

        com_pos_opensim
        com_vel_opensim
        com_acc_opensim
        res_f_ext_opensim

        com_pos_pelvis
        com_vel_pelvis
        com_acc_pelvis
        res_f_ext_pelvis

        com_pos_init
        com_vel_init

        n_samples
        n_horizons
        n_samples_horizon

        f_right
        f_left
        f_net_ext
        t_h

        % Instantenous acceleration
        inst_acc

        horizon_start
        acc_pred
        vel_pred
        pos_pred

        horizon_end
        dist_err
        rms_err
        max_err
        pred_traj
        ref_traj
        max_pred_traj
        min_pred_traj

        pred_err_fig

        movie_out

        shank_ang_r
        shank_ang_l

        thigh_ang_r
        thigh_ang_l

        hat_ang

        A_c
        B_c
        C_c

        A_d
        B_d
        C_d

        trial_part

    end

    methods 
        function obj = trial_prediction_class(trial_data, trial_part, duration_sec)
            % Constructor
            obj.trial_data = trial_data;
            obj.m = trial_data.m;

            len_m = length(trial_data.markerStruct.RIAS);
            len_f = length(trial_data.f_right);

            % marker_range = NaN(1,2);
            % force_range = NaN(1,2);

            obj.trial_part = trial_part;

            if (strcmp(trial_part, 'start'))

                marker_range_start = 1;
                marker_range_end = duration_sec*trial_data.rMakers + 1;

                force_range_start = 1;
                force_range_end = duration_sec*trial_data.rForces + trial_data.rForces/trial_data.rMakers;
                
            elseif (strcmp(trial_part, 'end'))

                marker_range_start = len_m - duration_sec*trial_data.rMakers - 1;
                marker_range_end = len_m;

                force_range_start = len_f - duration_sec*trial_data.rForces - trial_data.rForces/trial_data.rMakers;
                force_range_end = len_f;
        
            elseif (strcmp(trial_part, 'entire'))

               marker_range_start = 1;
               marker_range_end = len_m;

               force_range_start = 1;
               force_range_end = len_f;

            else

                error('Invalid part of the trial selected.')

            end

            marker_range = marker_range_start:marker_range_end;
            force_range = force_range_start:force_range_end;

            % All 
            obj.trial_data.f_right = trial_data.f_right(force_range,:);
            obj.trial_data.f_left = trial_data.f_left(force_range,:);

            obj.trial_data.com_pos = trial_data.com_pos(marker_range, :);
            obj.trial_data.com_vel = trial_data.com_vel(marker_range, :);
            
            obj.trial_data.markerStruct.RIAS = trial_data.markerStruct.RIAS(marker_range,:);
            obj.trial_data.markerStruct.LIAS = trial_data.markerStruct.LIAS(marker_range,:);
            obj.trial_data.markerStruct.RIPS = trial_data.markerStruct.RIPS(marker_range,:);
            obj.trial_data.markerStruct.LIPS = trial_data.markerStruct.LIPS(marker_range,:);
            
            obj.trial_data.ik_results = trial_data.ik_results(marker_range,:);


            obj.n_force_samples = length(force_range);
              
        end

        function init_pred_params(this)

            this.n_samples = length(1:this.sample_time_ms:this.n_force_samples);
            this.n_horizons = this.n_samples - 1;
            this.n_samples_horizon = this.horizon_time/this.sample_time_ms + 1;
            this.t_h = 0:this.sample_time_ms/1000:this.horizon_time/1000;

            this.set_external_forces();
            this.set_ref_com_state();
            this.set_init_com_state();
            
        end

        function set_external_forces(this)

            % External forces
            this.f_right = this.trial_data.f_right(1:this.sample_time_ms:this.n_force_samples,:)';
            this.f_left = this.trial_data.f_left(1:this.sample_time_ms:this.n_force_samples,:)';
            this.f_net_ext = this.f_right + this.f_left;
            this.f_net_ext(2,:) = this.f_net_ext(2,:) - this.m*this.g;

            temp = this.f_net_ext;
            this.f_net_ext(1,:) = this.f_net_ext(3,:);
            this.f_net_ext(3,:) = temp(1,:);

        end

        function estimate_com_opensim(this)

            this.com_pos_opensim(1,:) = this.trial_data.com_pos.center_of_mass_Z';
            this.com_pos_opensim(2,:) = this.trial_data.com_pos.center_of_mass_Y';
            this.com_pos_opensim(3,:) = this.trial_data.com_pos.center_of_mass_X';

            this.com_pos_opensim(1,:) = this.com_pos_opensim(1,:) - this.com_pos_opensim(1,1);
            this.com_pos_opensim(2,:) = this.com_pos_opensim(2,:) - this.com_pos_opensim(2,1);
            this.com_pos_opensim(3,:) = this.com_pos_opensim(3,:) - this.com_pos_opensim(3,1);

            this.com_vel_opensim(1,:) = this.trial_data.com_vel.center_of_mass_Z';
            this.com_vel_opensim(2,:) = this.trial_data.com_vel.center_of_mass_Y';
            this.com_vel_opensim(3,:) = this.trial_data.com_vel.center_of_mass_X';

            %smoothed_vel_x = smoothdata(this.trial_data.com_vel.center_of_mass_X', 'movmean', 5);
            %smoothed_vel_y = smoothdata(this.trial_data.com_vel.center_of_mass_Y', 'movmean', 5);
            %smoothed_vel_z = smoothdata(this.trial_data.com_vel.center_of_mass_Z', 'movmean', 5);

            smoothed_vel_x = smooth(this.trial_data.com_vel.center_of_mass_X', 'loess', 0.2);
            smoothed_vel_y = smooth(this.trial_data.com_vel.center_of_mass_Y', 'loess', 0.2);
            smoothed_vel_z = smooth(this.trial_data.com_vel.center_of_mass_Z', 'loess', 0.2);

            %sgolay_order = 1; 
            %sgolay_window = 11;
            %smoothed_vel_x = sgolayfilt(this.trial_data.com_vel.center_of_mass_X', sgolay_order, sgolay_window);
            %smoothed_vel_y = sgolayfilt(this.trial_data.com_vel.center_of_mass_Y', sgolay_order, sgolay_window);
            %smoothed_vel_z = sgolayfilt(this.trial_data.com_vel.center_of_mass_Z', sgolay_order, sgolay_window);

            %smoothed_vel_x = smoothdata(this.trial_data.com_vel.center_of_mass_X', 'gaussian', 30);
            %smoothed_vel_y = smoothdata(this.trial_data.com_vel.center_of_mass_Y', 'gaussian', 30);
            %smoothed_vel_z = smoothdata(this.trial_data.com_vel.center_of_mass_Z', 'gaussian', 30);

            %smoothed_vel_x = smooth(this.trial_data.com_vel.center_of_mass_X', 'lowess', 0.5);
            %smoothed_vel_y = smooth(this.trial_data.com_vel.center_of_mass_Y', 'lowess', 0.5);
            %smoothed_vel_z = smooth(this.trial_data.com_vel.center_of_mass_Z', 'lowess', 0.5);

            this.com_vel_opensim(1,:) = smoothed_vel_z;
            this.com_vel_opensim(2,:) = smoothed_vel_y;
            this.com_vel_opensim(3,:) = smoothed_vel_x;

            this.com_acc_opensim(1,:) = gradient(this.com_vel_opensim(1,:), this.sample_time_ms/1000);
            this.com_acc_opensim(2,:) = gradient(this.com_vel_opensim(2,:), this.sample_time_ms/1000);
            this.com_acc_opensim(3,:) = gradient(this.com_vel_opensim(3,:), this.sample_time_ms/1000);

            this.res_f_ext_opensim(1,:) = this.f_net_ext(1,:) - this.m*this.com_acc_opensim(1,:);
            this.res_f_ext_opensim(2,:) = this.f_net_ext(2,:) - this.m*this.com_acc_opensim(2,:);
            this.res_f_ext_opensim(3,:) = this.f_net_ext(3,:) - this.m*this.com_acc_opensim(3,:);

        end

        function estimate_com_pelvis(this)

            RIAS = this.trial_data.markerStruct.RIAS;
            LIAS = this.trial_data.markerStruct.LIAS;

            RIPS = this.trial_data.markerStruct.RIPS;
            LIPS = this.trial_data.markerStruct.LIPS;

            com_pos = 0.25*(RIAS + LIAS + RIPS + LIPS);
            com_pos(:,1) = com_pos(:,1) - com_pos(1,1);
            com_pos(:,2) = com_pos(:,2) - com_pos(1,2);
            com_pos(:,3) = com_pos(:,3) - com_pos(1,3);

            this.com_pos_pelvis(1,:) = com_pos(:,3)';
            this.com_pos_pelvis(2,:) = com_pos(:,2)';
            this.com_pos_pelvis(3,:) = com_pos(:,1)';

            this.com_vel_pelvis(1,:) = gradient(this.com_pos_pelvis(1,:), this.sample_time_ms/1000);
            this.com_vel_pelvis(2,:) = gradient(this.com_pos_pelvis(2,:), this.sample_time_ms/1000);
            this.com_vel_pelvis(3,:) = gradient(this.com_pos_pelvis(3,:), this.sample_time_ms/1000);

            %smoothed_vel_1 = smoothdata(this.com_vel_pelvis(1,:), 'movmean', 5);
            %smoothed_vel_2 = smoothdata(this.com_vel_pelvis(2,:), 'movmean', 5);
            %smoothed_vel_3 = smoothdata(this.com_vel_pelvis(3,:), 'movmean', 5);

            smoothed_vel_1 = smooth(this.com_vel_pelvis(1,:), 'loess', 0.2);
            smoothed_vel_2 = smooth(this.com_vel_pelvis(2,:), 'loess', 0.2);
            smoothed_vel_3 = smooth(this.com_vel_pelvis(3,:), 'loess', 0.2);

            %sgolay_order = 1; 
            %sgolay_window = 11;
            %smoothed_vel_1 = sgolayfilt(this.com_vel_pelvis(1,:), sgolay_order, sgolay_window);
            %smoothed_vel_2 = sgolayfilt(this.com_vel_pelvis(2,:), sgolay_order, sgolay_window);
            %smoothed_vel_3 = sgolayfilt(this.com_vel_pelvis(3,:), sgolay_order, sgolay_window);

            %smoothed_vel_1 = smoothdata(this.com_vel_pelvis(1,:), 'gaussian', 30);
            %smoothed_vel_2 = smoothdata(this.com_vel_pelvis(2,:), 'gaussian', 30);
            %smoothed_vel_3 = smoothdata(this.com_vel_pelvis(3,:), 'gaussian', 30);

            %smoothed_vel_1 = smooth(this.com_vel_pelvis(1,:), 'lowess', 0.5);
            %smoothed_vel_2 = smooth(this.com_vel_pelvis(2,:), 'lowess', 0.5);
            %smoothed_vel_3 = smooth(this.com_vel_pelvis(3,:), 'lowess', 0.5);

            this.com_vel_pelvis(1,:) = smoothed_vel_1;
            this.com_vel_pelvis(2,:) = smoothed_vel_2;
            this.com_vel_pelvis(3,:) = smoothed_vel_3;

            this.com_acc_pelvis(1,:) = gradient(this.com_vel_pelvis(1,:), this.sample_time_ms/1000);
            this.com_acc_pelvis(2,:) = gradient(this.com_vel_pelvis(2,:), this.sample_time_ms/1000);
            this.com_acc_pelvis(3,:) = gradient(this.com_vel_pelvis(3,:), this.sample_time_ms/1000);

            this.res_f_ext_pelvis(1,:) = this.f_net_ext(1,:) - this.m*this.com_acc_pelvis(1,:);
            this.res_f_ext_pelvis(2,:) = this.f_net_ext(2,:) - this.m*this.com_acc_pelvis(2,:);
            this.res_f_ext_pelvis(3,:) = this.f_net_ext(3,:) - this.m*this.com_acc_pelvis(3,:);

        end

        function set_ref_com_state(this)

            this.estimate_com_opensim();
            this.com_pos_ref = this.com_pos_opensim;
            this.com_vel_ref = this.com_vel_opensim;

        end

        function calc_segment_angles(this)

            this.shank_ang_r = this.trial_data.ik_results.ankle_angle_r;
            this.shank_ang_r = 90 + this.shank_ang_r;
            
            this.shank_ang_l = this.trial_data.ik_results.ankle_angle_l;
            this.shank_ang_l = 90 + this.shank_ang_l;
            
            figure;
            plot(this.shank_ang_r,'r');
            hold on;
            plot(this.shank_ang_l,'b');
            legend('Right shank', 'Left shank');
            xlabel(' time (samples)');
            ylabel('Angle (deg)');
            
            knee_ang_r = -this.trial_data.ik_results.knee_angle_r;
            knee_ang_l = -this.trial_data.ik_results.knee_angle_l;

            this.thigh_ang_r = knee_ang_r + this.shank_ang_r;
            this.thigh_ang_l = knee_ang_l + this.shank_ang_l;
            
            figure;
            plot(this.thigh_ang_r,'r');
            hold on;
            plot(this.thigh_ang_l,'b');
            legend('Right thigh', 'Left thigh');
            xlabel(' time (samples)');
            ylabel('Angle (deg)');
            
            this.hat_ang = this.trial_data.ik_results.lumbar_extension + this.trial_data.ik_results.pelvis_tilt;
            this.hat_ang = -this.hat_ang;
            this.hat_ang = 90 + this.hat_ang;
            
            figure;
            plot(this.hat_ang,'r');
            legend('HAT angle');
            xlabel(' time (samples)');
            ylabel('Angle (deg)');

            this.shank_ang_r = this.shank_ang_r*pi/180;
            this.shank_ang_l = this.shank_ang_l*pi/180;

            this.thigh_ang_r = this.thigh_ang_r*pi/180;
            this.thigh_ang_l = this.thigh_ang_l*pi/180;

            this.hat_ang = this.hat_ang*pi/180;
            
        end

        function y = calc_measurement(this, k, par)

            % q1 = this.shank_ang_r(k);
            % q2 = this.thigh_ang_r(k);
            % 
            % x_s_r = par.d_s*cos(q1);
            % y_s_r = par.d_s*sin(q1);
            % 
            % x_t_r = par.L_s*cos(q1) + par.d_t*cos(q2);
            % y_t_r = par.L_s*sin(q1) + par.d_t*sin(q2);
            % 
            % q1 = this.shank_ang_l(k);
            % q2 = this.thigh_ang_l(k);
            % 
            % x_s_l = par.d_s*cos(q1);
            % y_s_l = par.d_s*sin(q1);
            % 
            % x_t_l = par.L_s*cos(q1) + par.d_t*cos(q2);
            % y_t_l = par.L_s*sin(q1) + par.d_t*sin(q2);
            % 
            % % How to account for asymmetry?
            % q1 = (this.shank_ang_r(k) + this.shank_ang_l(k))/2;
            % q2 = (this.thigh_ang_r(k) + this.thigh_ang_l(k))/2;
            % q3 = this.hat_ang(k)*pi/180;
            % 
            % x_u = par.L_s*cos(q1) + par.L_t*cos(q2) + par.d_u*cos(q3);
            % y_u = par.L_s*sin(q1) + par.L_t*sin(q2) + par.d_u*sin(q3);
            % 
            y = NaN(3,1);
            % 
            % y(1) = (par.m_s*x_s_r + par.m_s*x_s_l + par.m_t*x_t_r + par.m_t*x_t_l + par.m_u*x_u)/par.m;
            % y(2) = (par.m_s*y_s_r + par.m_s*y_s_l + par.m_t*y_t_r + par.m_t*y_t_l + par.m_u*y_u)/par.m;
            % y(3) = 0;

            y(1) = this.com_pos_opensim(1,k);
            y(2) = this.com_pos_opensim(2,k);
            y(3) = this.com_pos_opensim(3,k);
            
        end

        function estimate_com_kalman(this)

            dt = 0.005;
            this.n_samples = length(this.shank_ang_l);

             this.A_c = [0, 1, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 0;
                    0, 0, 0, 1, 0, 0;
                    0, 0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0, 1;
                    0, 0, 0, 0, 0, 0]; 

            this.B_c = 1/this.m*[0, 0, 0;
                1, 0, 0;
                0, 0, 0;
                0, 1, 0;
                0, 0, 0;
                0, 0, 1];             

            this.C_c = [1,0,0,0,0,0;
                   0,0,1,0,0,0;
                   0,0,0,0,1,0]; 

            sys = ss(this.A_c,this.B_c,this.C_c,0);
            d_sys = c2d(sys,dt);

            this.A_d = d_sys.A;
            this.B_d = d_sys.B;
            this.C_d = d_sys.C;

            A = this.A_d;
            B = this.B_d;
            C = this.C_d;

            H = 1.80;

            par.m = this.m;

            L_f = 0.152*H;
            par.L_t = (0.53-0.285)*H;
            par.L_s = (0.285-0.039)*H;
            par.L_u = (0.818-0.53)*H;

            par.d_u = 0.626*par.L_u;
            par.d_t = 0.433*par.L_t;
            par.d_s = 0.433*par.L_s;

            par.h_f = 0.039*H;

            par.m_u = 0.678*par.m;
            par.m_t = 0.1*par.m;
            par.m_s = 0.0465*par.m;

            x_init = [0;0;0;0;0;0];
            y = this.calc_measurement(1, par);

            x_init(1) = y(1);
            x_init(3) = y(2);
            x_init(5) = y(3);

            x_init(2) = this.com_vel_ref(1,1);
            x_init(4) = this.com_vel_ref(2,1);
            x_init(6) = this.com_vel_ref(3,1);

            P_init = diag(0.1*[1, 1, 1, 1, 1, 1]);

            x_prev = x_init;
            P_prev = P_init;

            x_est = NaN(6,this.n_samples);

            n_states = length(x_init);

            P = cell(this.n_samples);

            % V_p: process noise covariance matrix
            % V_m: measurement noise covariance matrix
            V_p = diag([1,5,1,5,1,5]);
            V_m = diag([10,10,10]);

            K_cell = cell(this.n_samples);

            u(1,:) = this.f_net_ext(1,:);
            u(2,:) = this.f_net_ext(2,:);
            u(3,:) = this.f_net_ext(3,:);

            x_est(:,1) = x_init;

            for k = 2:this.n_samples

                % Estimate the state at k from the dynamics and control at k
                x_pred = A*x_prev + B*u(:,k-1);
                % Estimate the covariance matrix
                P_pred = A*P_prev*A.' + V_p;
                % Kalman gain
                K_cell{k} = P_pred*C.'*inv(C*P_pred*C.' + V_m);
                K = K_cell{k};
                % Obtain the measurements y
                y = this.calc_measurement(k, par);

                % Update the estimates.
                x_est(:,k) = x_pred + K*(y - C*x_pred);
                P{k} = (eye(n_states)-K*C)*P_pred;
                x_prev = x_est(:,k);
                P_prev = P{k};

            end

            this.com_pos_kalman(1,:) = x_est(1,:) - x_est(1,1);
            this.com_pos_kalman(2,:) = x_est(3,:) - x_est(3,1);
            this.com_pos_kalman(3,:) = x_est(5,:) - x_est(5,1);

            this.com_vel_kalman(1,:) = x_est(2,:);
            this.com_vel_kalman(2,:) = x_est(4,:);
            this.com_vel_kalman(3,:) = x_est(6,:);

            this.com_acc_kalman(1,:) = gradient(this.com_vel_kalman(1,:), this.sample_time_ms/1000);
            this.com_acc_kalman(2,:) = gradient(this.com_vel_kalman(2,:), this.sample_time_ms/1000);
            this.com_acc_kalman(3,:) = gradient(this.com_vel_kalman(3,:), this.sample_time_ms/1000);

            this.res_f_ext_kalman(1,:) = this.f_net_ext(1,:) - this.m*this.com_acc_kalman(1,:);
            this.res_f_ext_kalman(2,:) = this.f_net_ext(2,:) - this.m*this.com_acc_kalman(2,:);
            this.res_f_ext_kalman(3,:) = this.f_net_ext(3,:) - this.m*this.com_acc_kalman(3,:);

        end

        function estimate_com_integration(this)
            
            % t = this.trial_data.ik_results.time;
            % a = this.f_net_ext ./ this.m;
            % this.com_vel_integration = this.com_vel_ref(:,1) + cumtrapz(t, a, 2);
            % this.com_pos_integration = cumtrapz(t, this.com_vel_integration,2);

            x_init = [0;0;0;0;0;0];
            % y = this.calc_measurement(1, par);

            % x_init(1) = y(1);
            % x_init(3) = y(2);
            % x_init(5) = y(3);

            x_init(2) = this.com_vel_ref(1,1);
            x_init(4) = this.com_vel_ref(2,1);
            x_init(6) = this.com_vel_ref(3,1);

            u(1,:) = this.f_net_ext(1,:);
            u(2,:) = this.f_net_ext(2,:);
            u(3,:) = this.f_net_ext(3,:);

            x_est = NaN(6,this.n_samples);
            x_prev = x_init;

            x_est(:,1) = x_init;

            for k=2:this.n_samples

                x_est(:,k) = this.A_d*x_prev + this.B_d*u(:,k-1);
                x_prev = x_est(:,k);

            end

            this.com_pos_integration(1,:) = x_est(1,:) - x_est(1,1);
            this.com_pos_integration(2,:) = x_est(3,:) - x_est(3,1);
            this.com_pos_integration(3,:) = x_est(5,:) - x_est(5,1);

            this.com_vel_integration(1,:) = x_est(2,:);
            this.com_vel_integration(2,:) = x_est(4,:);
            this.com_vel_integration(3,:) = x_est(6,:);

        end

        function set_init_com_state(this)
            
            this.estimate_com_kalman();
            this.com_pos_init = this.com_pos_kalman;
            this.com_vel_init = this.com_vel_kalman;

        end

        % function pred_com_acc(this, inst_acc)
        % 
        % 
        % 
        % end

        function predict_com(this)

            this.horizon_start = zeros(this.n_horizons,1);

            this.acc_pred = cell(this.n_horizons,1);
            this.vel_pred = cell(this.n_horizons,1);
            this.pos_pred = cell(this.n_horizons,1);
            
            for i=1:this.n_horizons
            
                this.horizon_start(i) = i;
            
                %this.inst_acc = this.f_net_ext(:,i)/this.m;
                %this.inst_acc = this.m*this.com_acc_opensim(:,i)/this.m;
                this.inst_acc = this.m*this.com_acc_pelvis(:,i)/this.m;
            
                % this.acc_pred{i} = this.predic_com_acc
                this.acc_pred{i} = this.inst_acc .* ones(3, this.n_samples_horizon);
                this.vel_pred{i} = this.com_vel_init(:,i) + cumtrapz(this.t_h, this.acc_pred{i},2);
                this.pos_pred{i} = this.com_pos_init(:,i) + cumtrapz(this.t_h, this.vel_pred{i},2);
            
            end

        end


        function calc_pred_error(this)
            
            this.horizon_end = zeros(this.n_horizons,1);

            this.dist_err = cell(this.n_horizons,1);
            this.rms_err = zeros(this.n_horizons,1);
            this.max_err = zeros(this.n_horizons,1);
            
            this.pred_traj = cell(this.n_horizons,1);
            this.ref_traj = cell(this.n_horizons,1);
            this.max_pred_traj = zeros(3, this.n_horizons);
            this.min_pred_traj = zeros(3, this.n_horizons);
            
            
            for i=1:this.n_horizons
            
            if (this.horizon_start(i) + this.n_samples_horizon > this.n_samples)
                this.horizon_end(i) = this.n_samples;
            else
                this.horizon_end(i) = this.horizon_start(i) + this.n_samples_horizon - 1;
            end
            
            % horizon samples range
            range = this.horizon_start(i):this.horizon_end(i);
            % Needed for the shorter horizons at the end.
            horizon_length = length(range);
            this.ref_traj{i} = this.com_pos_ref(:,range);
            this.pred_traj{i} = this.pos_pred{i}(:,1:horizon_length);
            this.dist_err{i} = sqrt(sum((this.pred_traj{i} - this.ref_traj{i}).^2));
            this.rms_err(i) = rms(this.dist_err{i});
            this.max_err(i) = max(this.dist_err{i});
            this.max_pred_traj(:,i) = max(this.pred_traj{i},[],2);
            this.min_pred_traj(:,i) = min(this.pred_traj{i},[],2);
            
            end


        end


        function plot_pred_results(this)

            f = figure;
            f.Position = [1007         318        1198         898];
            % yyaxis left;
            box on; grid on; hold on;
            plot(this.rms_err, 'r');
            xlabel('Horizon number');
            % ylabel('RMS prediction error (m)');
            set(gca,'FontSize', 18);
            % set(gca, 'YColor', 'r'); 
            
            % f = figure;
            % yyaxis right;
            box on; grid on;
            plot(this.max_err, 'b');
            xlabel('Horizon number');
            % ylabel('Maximum prediction error (m)');
            set(gca,'FontSize', 22);
            % set(gca, 'YColor', 'b'); 
            ylabel('Prediction error (m)');

            legend('RMS error', 'Maximum error');

            file_name = strcat(this.trial_data.subject_id, '_', this.trial_data.trial_name, '_', ...
            this.trial_part, '_', num2str(this.horizon_time), 'ms', '_', 'const_', 'err'); 

            exportgraphics(f, strcat(file_name,'.png'));

        end

        function create_pred_animation(this, plot_title, plot_view)

            close all;
            fig = figure; hold on; grid on;
            fig.Position = [500 100 900 900];
            set(gcf,'color','w');
            marker_size = 12;
            line_width = 1.6;
            
            font_size = 24;
            
            xlabel('X (walkway length) (m)');
            ylabel('Z (walkway width) (m)');
            zlabel('Y (vertical) (m)');
            title(plot_title)
            set(gca,'Fontsize',font_size);

            % Initialize the movie structure.
            this.movie_out = struct('cdata', cell(1,this.n_horizons), 'colormap', cell(1,this.n_horizons));
            frame_counter = 1;
            
            % Object handles.
            h = zeros(1,1);
            h_counter = 1;
            
            delta = 0.1;

            % x-axis of the figure : walkway length (mocap x-axis) 
            x_plot_ind = 1;
            % y-axis of the figure : walkway width (mocap z-axis)
            y_plot_ind = 3;
            % z-axis of the figure : vertical (mocap y-axis)
            z_plot_ind = 2;
            
            x_lim_min = min([min(this.com_pos_ref(x_plot_ind,:)), ...
                min(this.min_pred_traj(x_plot_ind,:))]) - delta;
            x_lim_max = max([max(this.com_pos_ref(x_plot_ind,:)), ...
                max(this.max_pred_traj(x_plot_ind,:))]) + delta;
            
            y_lim_min = min([min(this.com_pos_ref(y_plot_ind,:)), ...
                min(this.min_pred_traj(y_plot_ind,:))]) - delta;
            y_lim_max = max([max(this.com_pos_ref(y_plot_ind,:)), ...
                max(this.max_pred_traj(y_plot_ind,:))]) + delta;
            
            z_lim_min = min([min(this.com_pos_ref(z_plot_ind,:)), ...
                min(this.min_pred_traj(z_plot_ind,:))]) - delta;
            z_lim_max = max([max(this.com_pos_ref(z_plot_ind,:)), ...
                max(this.max_pred_traj(z_plot_ind,:))]) + delta;
            
            for i = 1:this.n_horizons

                % i
            
                h(h_counter) = plot3(this.ref_traj{i}(x_plot_ind,1), ...
                    this.ref_traj{i}(y_plot_ind,1), ...
                    this.ref_traj{i}(z_plot_ind,1),'bo', 'MarkerSize',marker_size, 'LineWidth',line_width); 
                h_counter = h_counter + 1;

                h(h_counter) = plot3(this.ref_traj{i}(x_plot_ind,:), ...
                    this.ref_traj{i}(y_plot_ind,:), ...
                    this.ref_traj{i}(z_plot_ind,:), ...
                    'b--','LineWidth',line_width); h_counter_ref = h_counter;
                h_counter = h_counter + 1;

                h(h_counter) = plot3(this.ref_traj{i}(x_plot_ind,end), ...
                    this.ref_traj{i}(y_plot_ind,end), ...
                    this.ref_traj{i}(z_plot_ind,end), ...
                    'bd', 'MarkerSize',marker_size, 'LineWidth',line_width); 
                h_counter = h_counter + 1; 
            
                h(h_counter) = plot3(this.pred_traj{i}(x_plot_ind,1), ...
                    this.pred_traj{i}(y_plot_ind,1), ...
                    this.pred_traj{i}(z_plot_ind,1), ...
                    'ro', 'MarkerSize',marker_size, 'LineWidth',line_width); 
                h_counter = h_counter + 1;

                h(h_counter) = plot3(this.pred_traj{i}(x_plot_ind,:), ...
                    this.pred_traj{i}(y_plot_ind,:), ...
                    this.pred_traj{i}(z_plot_ind,:), ...
                    'r--','LineWidth',line_width); h_counter_pred = h_counter;
                h_counter = h_counter + 1;

                h(h_counter) = plot3(this.pred_traj{i}(x_plot_ind,end), ...
                    this.pred_traj{i}(y_plot_ind,end), ...
                    this.pred_traj{i}(z_plot_ind,end), ...
                    'rd', 'MarkerSize',marker_size, 'LineWidth',line_width); 
            
                axis([x_lim_min x_lim_max y_lim_min y_lim_max z_lim_min z_lim_max]);
                view(plot_view);
                legend([h(h_counter_pred),h(h_counter_ref)],'Prediction','OpenSim','Location','northeast','FontSize',24);
            
                this.movie_out(frame_counter) = getframe(fig);
                frame_counter = frame_counter + 1;
            
                % Delete the graphical objects of this frame. 
                delete(h);
                % Reinitialize h because the number of objects in the next frame coud be different.
                h = zeros(1,1);
                h_counter = 1;
            
            end

            sample_rate = 200;
            slowdown = 1;
            fps = sample_rate/slowdown;
            movie(this.movie_out,1,fps)

        end

        function save_pred_animation(this, view_name, frame_rate, quality)

            % AVI supports high frame rates
            file_name = strcat(this.trial_data.subject_id, '_', this.trial_data.trial_name, '_', ...
            this.trial_part, '_', num2str(this.horizon_time), 'ms', '_', 'const_', view_name); 

            v = VideoWriter(strcat(file_name,'.avi'));
            v.FrameRate = frame_rate;
            v.Quality = quality;
            open(v)
            writeVideo(v,this.movie_out)
            close(v);

        end

        function plot_com_estimates(this)

            f = figure;
            f.Position = [936         128        1437        1123];
            tcl = tiledlayout(3,1);
            nexttile(tcl);
            p1 = plot(this.com_pos_kalman(1,:),'r'); 
            hold on; grid on;
            p2 = plot(this.com_pos_opensim(1,:),'b');
            p3 = plot(this.com_pos_pelvis(1,:),'k'); 
            ylabel('X pos. (walkway length) (m)');
            % legend('Kalman', 'OpenSim', 'Integration','Location','best');
            set(gca,'FontSize', 18);
            xlim([1,this.n_samples]);
            
            nexttile(tcl);
            plot(this.com_pos_kalman(2,:),'r'); 
            hold on; grid on;
            plot(this.com_pos_opensim(2,:),'b');
            plot(this.com_pos_pelvis(2,:),'k'); 
            ylabel('Y pos. (vertical) (m)');
            % legend('Kalman', 'OpenSim', 'Integration','Location','best');
            set(gca,'FontSize', 18);
            xlim([1,this.n_samples]);
            
            nexttile(tcl);
            plot(this.com_pos_kalman(3,:),'r'); 
            hold on; grid on;
            plot(this.com_pos_opensim(3,:),'b');
            plot(this.com_pos_pelvis(3,:),'k'); 
            ylabel('Z pos. (walkway width) (m)');
            % legend('Kalman', 'OpenSim', 'Integration','Location','best');
            set(gca,'FontSize', 18);
            xlabel('Time (sample)');
            xlim([1,this.n_samples]);
            
            hL = legend([p1,p2,p3],'Kalman', 'OpenSim', 'Pelvis', 'Location','Best','Orientation','horizontal','FontSize',16); 
            % Move the legend to the right side of the figure
            hL.Layout.Tile = 'South';

            exportgraphics(f, 'CoM_pos_estimation_comparison.png');
            
            f = figure;
            f.Position = [936         128        1437        1123];
            tcl = tiledlayout(3,1);
            nexttile(tcl);
            p1 = plot(this.com_vel_kalman(1,:),'r'); 
            hold on; grid on;
            p2 = plot(this.com_vel_opensim(1,:),'b');
            p3 = plot(this.com_vel_pelvis(1,:),'k'); 
            ylabel('X vel. (walkway length) (m/s)');
            % legend('Kalman', 'OpenSim', 'Integration','Location','best');
            set(gca,'FontSize', 18);
            xlim([1,this.n_samples]);
            
            nexttile(tcl);
            plot(this.com_vel_kalman(2,:),'r'); 
            hold on; grid on;
            plot(this.com_vel_opensim(2,:),'b');
            plot(this.com_vel_pelvis(2,:),'k'); 
            ylabel('Y vel. (vertical) (m/s)');
            % legend('Kalman', 'OpenSim', 'Integration','Location','best');
            set(gca,'FontSize', 18);
            xlim([1,this.n_samples]);
            
            nexttile(tcl);
            plot(this.com_vel_kalman(3,:),'r'); 
            hold on; grid on;
            plot(this.com_vel_opensim(3,:),'b');
            plot(this.com_vel_pelvis(3,:),'k'); 
            ylabel('Z vel. (walkway width) (m/s)');
            % legend('Kalman', 'OpenSim', 'Integration','Location','best');
            set(gca,'FontSize', 18);
            xlabel('Time (sample)');
            xlim([1,this.n_samples]);
            
            hL = legend([p1,p2,p3],'Kalman', 'OpenSim', 'Pelvis', 'Location','Best','Orientation','horizontal','FontSize',16); 
            % Move the legend to the right side of the figure
            hL.Layout.Tile = 'South';

            exportgraphics(f, 'CoM_vel_estimation_comparison.png');

        end

        function plot_f_ext_res(this)

            f = figure;
            f.Position = [936         128        1437        1123];
            tcl = tiledlayout(3,1);
            nexttile(tcl);
            p1 = plot(this.res_f_ext_kalman(1,:),'r'); 
            hold on; grid on;
            p2 = plot(this.res_f_ext_opensim(1,:),'b');
            p3 = plot(this.res_f_ext_pelvis(1,:),'k'); 
            ylabel('X pos. (walkway length) (m)');
            % legend('Kalman', 'OpenSim', 'Integration','Location','best');
            set(gca,'FontSize', 18);
            xlim([1,this.n_samples]);
            
            nexttile(tcl);
            plot(this.res_f_ext_kalman(2,:),'r'); 
            hold on; grid on;
            plot(this.res_f_ext_opensim(2,:),'b');
            plot(this.res_f_ext_pelvis(2,:),'k'); 
            ylabel('Y pos. (vertical) (m)');
            % legend('Kalman', 'OpenSim', 'Integration','Location','best');
            set(gca,'FontSize', 18);
            xlim([1,this.n_samples]);
            
            nexttile(tcl);
            plot(this.res_f_ext_kalman(3,:),'r'); 
            hold on; grid on;
            plot(this.res_f_ext_opensim(3,:),'b');
            plot(this.res_f_ext_pelvis(3,:),'k'); 
            ylabel('Z pos. (walkway width) (m)');
            % legend('Kalman', 'OpenSim', 'Integration','Location','best');
            set(gca,'FontSize', 18);
            xlabel('Time (sample)');
            xlim([1,this.n_samples]);
            
            hL = legend([p1,p2,p3],'Kalman', 'OpenSim', 'Pelvis', 'Location','Best','Orientation','horizontal','FontSize',16); 
            % Move the legend to the right side of the figure
            hL.Layout.Tile = 'South';

            exportgraphics(f, 'f_ext_res_comparison.png');

        end

    end

end