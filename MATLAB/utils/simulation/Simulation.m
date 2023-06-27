classdef Simulation
    properties
        Nsim
        prob_info
        rand_seed = 42
    end
    methods
        function self = Simulation(Nsim, prob_info)
            % constructor for Simulation object
            self.Nsim = Nsim;
            self.prob_info = prob_info;
            if isfield(self.prob_info, 'rand_seed')
                self.rand_seed = self.prob_info.rand_seed;
            end
        end
        
        function [self, sim_data] = run_closed_loop(self, c, tt, observer, CEM, multistage, offset)
            % check controller type
            if isa(c,'MPC')
                mpc_controller = true;
                hw_controller = false;
                disp('Using a MPC.')
            elseif isa(c,'DNN')
                mpc_controller = false;
                hw_controller = false;
                disp('Using an approximate controller.')
            elseif isa(c, 'HWController')
                mpc_controller = false;
                hw_controller = true;
                disp('Using an approximate controller on hardware.')
            else
                error('Unsupported controller type.')
            end
            % extract relevant information from prob_info
            nu = self.prob_info.nu;
            nx = self.prob_info.nx;
            ny = self.prob_info.ny;
            nyc = self.prob_info.nyc;
            nd = self.prob_info.nd;
            nv = self.prob_info.nv;
            nw = self.prob_info.nw;
            
            x0 = self.prob_info.x0;
            myref = self.prob_info.myref;
            ts = self.prob_info.ts;
            
            v_min = self.prob_info.v_min;
            v_max = self.prob_info.v_max;
            w_min = self.prob_info.w_min;
            w_max = self.prob_info.w_max;
            if multistage
                Wset = self.prob_info.Wset;
            end
            
            % check if a plant model was provided
            plant_provided = false;
            if isfield(self.prob_info, 'fp')
                plant_provided = true;
                disp('Using user-provided plant...')
                fp = self.prob_info.fp;
                hp = self.prob_info.hp;
            end
            if ~plant_provided
                error('no plant provided')
            end
            
            % set random seed (for reproducible results)
            rng(self.rand_seed)
            % generate noise/disturbance variables
            Vsim = v_min + (v_max-v_min) .* rand(nv,self.Nsim+1);
            Wsim = w_min + (w_max-w_min) .* rand(nw,self.Nsim+1);
            
            % instantiate container variables for storing simulation data
            Xsim = zeros(nx,self.Nsim+1); % state trajectories (plant)
            Ysim = zeros(ny,self.Nsim+1); % output trajectories (plant)
            Usim = zeros(nu,self.Nsim);   % input trajectories (plant)
            Xhat = zeros(size(Xsim));  % state estimation
            Dhat = zeros(nd, self.Nsim+1);  % disturbance estimation
            Ypred = zeros(ny,self.prob_info.Np,self.Nsim);
            
            if offset
                Xss_sim = zeros(size(Xsim));
                Yss_sim = zeros(size(Ysim));
                Uss_sim = zeros(size(Usim));
            end

            ctime = zeros(1,self.Nsim);   % computation time
            Jsim = zeros(1,self.Nsim);    % cost/optimal objective value (controller)
            Yrefsim = zeros(nyc,self.Nsim);  % output reference/target (as sent to controller)

            if CEM
                CEMsim = zeros(1,self.Nsim+1);
                CEMadd = self.prob_info.CEMadd;
            end
            
            % set initial states and reset controller for consistency
            Xsim(:,1) = x0;
            Xhat(:,1) = Xsim(:,1);
            if ~isempty(observer)
                observer.xhat = Xhat(:,1);
                if offset
                    observer.dhat = zeros(nd,1);
                end
            end
            if plant_provided
                Ysim(:,1) = full(hp(Xhat(:,1),Vsim(:,1)));
            else
                Ysim(:,1) = 37;
            end
            if mpc_controller
                c = c.reset_initial_guesses();
                if offset && strcmpi(c.offset_style, 'rawlings')
                    tt.reset_initial_guesses();
                end
            end

            % loop over simulation time
            for k = 1:self.Nsim
                startTime = tic;

                Yrefsim(:,k) = myref(k*ts);
                if mpc_controller
                    % get optimal input from a MPC
                    if CEM
                        if multistage
                            % multistage version
                            c = c.set_parameters({Xhat(:,k), Yrefsim(:,k), CEMsim(:,k), Wset}, false);
                        else
                            % economic version
                            c = c.set_parameters({Xhat(:,k), Yrefsim(:,k), CEMsim(:,k)}, false);
                        end
                    else
                        if offset
                            % offset-free version
                            if strcmpi(c.offset_style, 'limon')
                                c = c.set_parameters({Xhat(:,k), Dhat(:,k), Yrefsim(:,k)}, false);
                            elseif strcmpi(c.offset_style, 'rawlings')
                                % target tracker
                                tt = tt.set_parameters({Yrefsim(:,k), Dhat(:,k)}, false);
                                [tt, Xss_sim(:,k), Uss_sim(:,k)] = tt.solve_ocp();
                                
                                c = c.set_parameters({Xhat(:,k), Dhat(:,k), Xss_sim(:,k), Uss_sim(:,k)}, false);
                            else
                                error('Controller type not supported.')
                            end
                        else
                            % nominal version
                            c = c.set_parameters({Xhat(:,k), Yrefsim(:,k)}, false);
                        end
                    end
                    [c, res, feas] = c.solve_mpc();
                    Uopt = res.U;
                    Jopt = res.J;
                else
                    % get optimal input from some other controller type
                    % (i.e., an approximate controller (DNN))
                    Jopt = NaN;
                    if CEM
                        in_vec = [Xhat(:,k); CEMsim(:,k)];
                    else
                        in_vec = [Xhat(:,k); Yrefsim(:,k)];
                    end
                    if hw_controller
                        [Uopt, c] = c.getControlInput(in_vec);
                    else
                        if isempty(c.netca)
                            Uopt = c.net(in_vec);
                        else
                            Uopt = full(c.netca(in_vec));
                        end
                    end
                end

                ctime(k) = toc(startTime);
                if mpc_controller
                    if ~feas
                        fprintf('Was not feasible on iteration %d of simulation\n',k)
                    end
                end

                Usim(:,k) = Uopt;
                Jsim(k) = Jopt;

                % send optimal input to plant/real system
                if plant_provided
                    Xsim(:,k+1) = full(fp(Xsim(:,k),Usim(:,k),Wsim(:,k)));
                    Ysim(:,k+1) = full(hp(Xsim(:,k+1),Vsim(:,k+1)));
                end
                if CEM
                    CEMsim(:,k+1) = CEMsim(:,k) + full(CEMadd(Ysim(:,k+1)));
                end

                % get estimates
                if ~isempty(observer)
                    if offset
                        [observer, xhat, dhat] = observer.update_observer(Usim(:,k), Ysim(:,k+1));
                        Xhat(:,k+1) = xhat;
                        Dhat(:,k+1) = dhat;
                    else
                        [observer, xhat] = observer.update_observer(Usim(:,k), Ysim(:,k+1));
                        Xhat(:,k+1) = xhat;
                    end
                else
                    Xhat(:,k+1) = Xsim(:,k+1);
                end
            end

            % create dictionary of simulation data
            sim_data = struct();
            sim_data.Xsim = Xsim;
            sim_data.Ysim = Ysim;
            sim_data.Usim = Usim;
            sim_data.Jsim = Jsim;
            sim_data.Yrefsim = Yrefsim;
            sim_data.ctime = ctime;
            sim_data.Xhat = Xhat;
            sim_data.Ypred = Ypred;
            if offset
                sim_data.Xss_sim = Xss_sim;
                sim_data.Uss_sim = Uss_sim;
                sim_data.Yss_sim = Yss_sim;
            end
            if CEM
                sim_data.CEMsim = CEMsim;
            end
        end
    end
end