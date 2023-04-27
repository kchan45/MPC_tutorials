classdef MultistageMPC < MPC
    % This class defines an MPC implementation for the thermal dose
    % delivery of plasma jets. This class utilizes the Opti stack interface
    % of CasADi.
    properties
        opti;
        opti_vars;
        opti_params;
        mpc_type = 'multistage';
    end
    methods
        function self = get_mpc(self)
            % This method creates the optimization problem for the MPC. All
            % information necessary for the creation of this controller is
            % passed upon instantiation of this object within the prob_info
            % dict. For more details on the optimization problem, the user
            % is referred to the paper associated with the release of this
            % code.
            % 
            % This code uses IPOPT for the NLP solver which is distributed
            % with CasADi. Users are referred to IPOPT
            % [https://coin-or.github.io/Ipopt/] and the associated paper
            % for more information on this solver.
            
            import casadi.*
            % unpack problem information
            Np = self.prob_info.Np;
            N_robust = self.prob_info.N_robust;
            case_idx = self.prob_info.case_idx;
            Wset = self.prob_info.Wset;

            nu = self.prob_info.nu;
            nx = self.prob_info.nx;
            ny = self.prob_info.ny;
            nyc = self.prob_info.nyc;
            nw = self.prob_info.nw;

            u_min = self.prob_info.u_min;
            u_max = self.prob_info.u_max;
            x_min = self.prob_info.x_min;
            x_max = self.prob_info.x_max;
            y_min = self.prob_info.y_min;
            y_max = self.prob_info.y_max;

            u_init = self.prob_info.u_init;
            x_init = self.prob_info.x_init;
            y_init = self.prob_info.y_init;

            f = self.prob_info.f;
            h = self.prob_info.h;
            r = self.prob_info.r;
            lstage = self.prob_info.lstage;
            
            % Build the scenario matrix with all the combinations.
            % Dimension (Ncases^Nrobust, Nrobust)
            scenario_mat = combvec(case_idx, case_idx)';
            N_scenarios = length(scenario_mat);

            % Weights for cost function
            w_i = (1/N_scenarios)*ones(N_scenarios,1);
            
            % create NLP opti object
            self.opti = Opti();
            
            % Initialize container lists for all states, inputs, outputs, and
            % predicted noise over horizon
            X = cell(Np+1, N_scenarios);
            Y = cell(Np+1, N_scenarios);
            U = cell(Np, N_scenarios);
            wPred = cell(Np, N_scenarios);

            J = 0; % initialize cost/objective function
            n_disc = 0;
            
            % define parameter(s), variable(s), and problem
            CEMref = self.opti.parameter(nyc); % target/reference output
            self.opti.set_value(CEMref, zeros(nyc,1))

            CEM0 = self.opti.parameter(nyc); % initial CEM
            self.opti.set_value(CEM0, zeros(nyc,1))
            
            for n_sc = 1:N_scenarios
                sc_vec = scenario_mat(n_sc,:)';
                w_idx = sc_vec - min(min(scenario_mat)) + 1;
            
                X{1,n_sc} = self.opti.parameter(nx); % initial state as a parameter
                self.opti.set_value(X{1,n_sc}, zeros(nx,1))

                Y{1,n_sc} = self.opti.variable(ny); % initial output variable
                self.opti.subject_to(Y{1,n_sc} == h(X{1,n_sc}))
                self.opti.set_initial(Y{1,n_sc}, y_init)
                n_disc = n_disc + ny;

                % the loop below systematically defines the variables of the
                % optimal control problem (OCP) over the prediction horizon
                for k = 1:Np
                    % variable and constraints for u_{k}
                    U{k,n_sc} = self.opti.variable(nu);
                    self.opti.subject_to(u_min <= U{k,n_sc} <= u_max)
                    self.opti.set_initial(U{k,n_sc}, u_init)
                    n_disc = n_disc + nu;
                    
                    % noise prediction parameter
                    wPred{k,n_sc} = self.opti.parameter(nw);
                    if k <= N_robust
                        self.opti.set_value(wPred{k,n_sc}, Wset(:,w_idx(k)));
                    elseif k < Np
                        self.opti.set_value(wPred{k,n_sc}, Wset(:,w_idx(N_robust)));
                    else
                        self.opti.set_value(wPred{k,n_sc}, zeros(nw,1));
                    end

                    Jstage = lstage(X{k,n_sc}, wPred{k,n_sc});
                    J = J + w_i(n_sc)*Jstage; % add to the cost (construction of the objective)

                    % variable x_{k+1}
                    X{k+1,n_sc} = self.opti.variable(nx);
                    self.opti.subject_to(x_min <= X{k+1,n_sc} <= x_max)
                    self.opti.set_initial(X{k+1,n_sc}, x_init)
                    n_disc = n_disc + nx;

                    % variable y_{k+1}
                    Y{k+1,n_sc} = self.opti.variable(ny);
                    self.opti.subject_to(y_min <= Y{k+1,n_sc} <= y_max)
                    self.opti.set_initial(Y{k+1,n_sc}, y_init)
                    n_disc = n_disc + ny;

                    % dynamics constraint
                    self.opti.subject_to(X{k+1,n_sc} == f(X{k,n_sc},U{k,n_sc},wPred{k,n_sc}));

                    % output equality constraint
                    self.opti.subject_to(Y{k+1,n_sc} == h(X{k+1,n_sc}));
                end

                % terminal cost and constraints
                J_end = lstage(X{end,n_sc}, zeros(nw,1));
                J = J+w_i(n_sc)*J_end;
            end

            Jcon = J + CEM0;
            J = (Jcon-CEMref)^2;
            
            % Non-anticipativity constraint (first split)
            for jj = 1:N_scenarios-1
                self.opti.subject_to(U{1,jj} == U{1,jj+1});
            end
            
            % Non-anticipativity constraint (second split)
            if N_robust > 1
                n_usc = length(case_idx); % number of unique scenarios at each node
                scenario_vec = 1:N_scenarios;
                for ii = 1:n_usc
                    % get the scenarios that correspond to the same parent node
                    same_nodes = scenario_vec(scenario_mat(:,1) == case_idx(ii));

                    % add constraint
                    for jj = 1:length(same_nodes)-1
                        self.opti.subject_to(U{2,same_nodes(jj)} == U{2,same_nodes(jj+1)});
                    end
                end
            end
            
            % set to minimize objective/cost
            self.opti.minimize( J )
            
            % set solver and problem options
            discrete = zeros(n_disc,1);
            p_opts = struct('verbose', 0, ...
                            'expand', true, ...
                            'discrete', discrete, ...
                            'print_time', 0); % problem options
            s_opts = struct('max_iter', 1000, ...
                            'print_level', 0, ...
                            'tol', 1e-6);     % solver options
            self.opti.solver('ipopt', p_opts, s_opts) % add the solver to the opti object

            % save list containers of variables/parameters into a dict for
            % portability
            self.opti_vars = struct;
            self.opti_vars.U = U;
            self.opti_vars.X = X(2:end,:);
            self.opti_vars.Y = Y;
            self.opti_vars.J = J;

            self.opti_params  = struct;
            self.opti_params.X0 = X(1,:);
            self.opti_params.CEMref = CEMref;
            self.opti_params.CEM0 = CEM0;
            self.opti_params.wPred = wPred;

        end
        
        function self = reset_initial_guesses(self, varargin)
            % This method resets the intial guesses for the variables of
            % the optimization problem back to those defined in the
            % problem_info dict provided upon instantiation of the
            % NonlinearMPC object.
        
            % unpack relevant information from the prob_info dict
            Np = self.prob_info.Np;
            N_robust = self.prob_info.N_robust;
            case_idx = self.prob_info.case_idx;
            u_init = self.prob_info.u_init;
            x_init = self.prob_info.x_init;
            y_init = self.prob_info.y_init;

            % unpack relevant variable containers from problem creation
            U = self.opti_vars.U;
            X = self.opti_vars.X;
            Y = self.opti_vars.Y;
            
            N_scenarios = length(case_idx)^N_robust;

            for n_sc = 1:N_scenarios
                self.opti.set_initial(Y{1,n_sc}, y_init)
                for k = 1:Np
                    self.opti.set_initial(U{k,n_sc}, u_init)
                    self.opti.set_initial(X{k,n_sc}, x_init)
                    self.opti.set_initial(Y{k+1,n_sc}, y_init)
                end
            end
        end
        
        function self = set_parameters(self, params_list, warn)
            % This method sets the values of the parameters of the
            % optimization problem to those provided as arguments to this
            % method.
            
            Np = self.prob_info.Np;
            nw = self.prob_info.nw;
            N_robust = self.prob_info.N_robust;
            case_idx = self.prob_info.case_idx;
            % Build the scenario matrix with all the combinations.
            % Dimension (Ncases^Nrobust, Nrobust)
            scenario_mat = combvec(case_idx, case_idx)';
            N_scenarios = length(scenario_mat);
            
            % unpack params_list argument
            if isequal(class(params_list), 'struct')
                x0 = params_list.X0;
                cemref = params_list.CEMref;
                cem0 = params_list.CEM0;
                Wset = params_list.Wset;
            elseif isequal(class(params_list), 'cell')
                if warn
                    warning('You passed in a cell array. Make sure the cell array order is the same as the following: ')
                    disp(fieldnames(self.opti_params))
                end
                x0 = params_list{1};
                cemref = params_list{2};
                cem0 = params_list{3};
                Wset = params_list{4};
            else
                warning('Unsupported params_list class... Parameters were not set!')
            end
            self.opti.set_value(self.opti_params.CEMref, cemref)
            self.opti.set_value(self.opti_params.CEM0, cem0)
            
            % update intial state and uncertainty bounds
            for n_sc = 1:N_scenarios
                % reset initial conditions
                self.opti.set_value(self.opti_params.X0{n_sc}, x0)
                
                sc_vec = scenario_mat(n_sc,:)';
                w_idx = sc_vec - min(min(scenario_mat)) + 1;
                for k = 1:Np
                    % change bounds depending on robust horizon
                    if k <= N_robust
                        self.opti.set_value(self.opti_params.wPred{k,n_sc}, Wset(:,w_idx(k)));
                    elseif k < Np
                        self.opti.set_value(self.opti_params.wPred{k,n_sc}, Wset(:,w_idx(N_robust)));
                    else
                        self.opti.set_value(self.opti_params.wPred{k,n_sc}, zeros(nw,1));
                    end
                end
            end
        end
        
        function [self, res, feas] = solve_mpc(self)
            % This method solves the MPC as defined by the get_mpc() method
            % of this class. This method can only be called after the the
            % get_mpc() method has been called (i.e., the optimization
            % problem must be defined before it can be solved).
            
            % unpack relevant information from problem creation
            u_min = self.prob_info.u_min;
            u_max = self.prob_info.u_max;

            feas = 1;
            res = struct();
            try
                sol = self.opti.solve();
                res.U = sol.value(self.opti_vars.U{1,1});
                res.J = sol.value(self.opti_vars.J);

                if self.prob_info.warm_start
                    self.opti.set_initial(sol.value_variables())
                end

            catch
                % if solve fails, get the last value
                feas = 0;
                
                u = self.opti.debug.value(self.opti_vars.U{1,1});
                res.U = max(min(u,u_max),u_min);
                res.J = self.opti.debug.value(self.opti_vars.J);

                % disp('U_0:'); disp(res.U)
                % disp('J:'); disp(res.J)
            end
        end
    end
end