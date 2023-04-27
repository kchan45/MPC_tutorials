classdef EconomicMPC < MPC
    % This class defines a MPC implementation for the thermal dose delivery
    % of plasma jets. This class utilizes the Opti stack interface of
    % CasADi.
    properties
        opti;
        opti_vars;
        opti_params;
        mpc_type = 'economic';
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

            nu = self.prob_info.nu;
            nx = self.prob_info.nx;
            ny = self.prob_info.ny;
            nyc = self.prob_info.nyc;

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
%             r = self.prob_info.r;
            lstage = self.prob_info.lstage;
            
            % create NLP opti object
            self.opti = Opti();
            p_opts = struct('verbose', 0, ...
                            'expand', true, ...
                            'print_time', 0); % problem options
            s_opts = struct('max_iter', 1000, ...
                            'print_level', 0, ...
                            'tol', 1e-6);     % solver options
            self.opti.solver('ipopt', p_opts, s_opts) % add the solver to the opti object

            % Initialize container lists for all states, inputs, outputs, and
            % predicted noise over horizon
            X = cell(Np+1,1);
            Y = cell(Np+1,1);
            U = cell(Np,1);

            J = 0; % initialize cost/objective function

            % define parameter(s), variable(s), and problem
            X{1} = self.opti.parameter(nx); % initial state as a parameter
            self.opti.set_value(X{1}, zeros(nx,1))

            CEMref = self.opti.parameter(nyc); % target/reference output
            self.opti.set_value(CEMref, zeros(nyc,1))

            CEM0 = self.opti.parameter(nyc); % initial CEM
            self.opti.set_value(CEM0, zeros(nyc,1))

            Y{1} = self.opti.variable(ny); % initial output variable
            self.opti.subject_to(Y{1} == h(X{1}))
            self.opti.set_initial(Y{1}, y_init)

            % the loop below systematically defines the variables of the
            % optimal control problem (OCP) over the prediction horizon
            for k = 1:Np
                % variable and constraints for u_{k}
                U{k} = self.opti.variable(nu);
                self.opti.subject_to(u_min <= U{k} <= u_max)
                self.opti.set_initial(U{k}, u_init)

                Jstage = lstage(X{k});
                J = J + Jstage; % add to the cost (construction of the objective)

                % variable x_{k+1}
                X{k+1} = self.opti.variable(nx);
                self.opti.subject_to(x_min <= X{k+1} <= x_max)
                self.opti.set_initial(X{k+1}, x_init)

                % variable y_{k+1}
                Y{k+1} = self.opti.variable(ny);
                self.opti.subject_to(y_min <= Y{k+1} <= y_max)
                self.opti.set_initial(Y{k+1}, y_init)

                % dynamics constraint
                self.opti.subject_to(X{k+1} == f(X{k},U{k}))

                % output equality constraint
                self.opti.subject_to(Y{k+1} == h(X{k+1}))
            end

            % terminal cost and constraints
            J_end = lstage(X{end});
            J = J+J_end;

            Jcon = J + CEM0;
            J = (Jcon-CEMref)^2;

            % set to minimize objective/cost
            self.opti.minimize( J )

            % save list containers of variables/parameters into a dict for
            % portability
            self.opti_vars = struct;
            self.opti_vars.U = U;
            self.opti_vars.X = X(2:end);
            self.opti_vars.Y = Y;
            self.opti_vars.J = J;

            self.opti_params  = struct;
            self.opti_params.X0 = X{1};
            self.opti_params.CEMref = CEMref;
            self.opti_params.CEM0 = CEM0;

        end
        
        function self = reset_initial_guesses(self, varargin)
            % This method resets the intial guesses for the variables of
            % the optimization problem back to those defined in the
            % problem_info dict provided upon instantiation of the
            % NonlinearMPC object.
        
            % unpack relevant information from the prob_info dict
            Np = self.prob_info.Np;
            u_init = self.prob_info.u_init;
            x_init = self.prob_info.x_init;
            y_init = self.prob_info.y_init;

            % unpack relevant variable containers from problem creation
            U = self.opti_vars.U;
            X = self.opti_vars.X;
            Y = self.opti_vars.Y;

            self.opti.set_initial(Y{1}, y_init)
            for k = 1:Np
                self.opti.set_initial(U{k}, u_init)
                self.opti.set_initial(X{k}, x_init)
                self.opti.set_initial(Y{k+1}, y_init)
            end
        end
        
        function self = set_parameters(self, params_list, warn)
            % This method sets the values of the parameters of the
            % optimization problem to those provided as arguments to this
            % method.
            
            if isequal(class(params_list), 'struct')
                self.opti.set_value(self.opti_params.X0, params_list.X0)
                self.opti.set_value(self.opti_params.CEMref, params_list.CEMref)
                self.opti.set_value(self.opti_params.CEM0, params_list.CEM0)
            elseif isequal(class(params_list), 'cell')
                if warn
                    warning('You passed in a cell array. Make sure the cell array order is the same as the following: ')
                    disp(fieldnames(self.opti_params))
                end
                self.opti.set_value(self.opti_params.X0, params_list{1})
                self.opti.set_value(self.opti_params.CEMref, params_list{2})
                self.opti.set_value(self.opti_params.CEM0, params_list{3})
            else
                warning('Unsupported params_list class... Parameters were not set!')
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
                res.U = sol.value(self.opti_vars.U{1});
                res.J = sol.value(self.opti_vars.J);

                if self.prob_info.warm_start
                    self.opti.set_initial(sol.value_variables())
                end

            catch
                % if solve fails, get the last value
                feas = 0;
                
                u = self.opti.debug.value(self.opti_vars.U{1});
                res.U = max(min(u,u_max),u_min);
                res.J = self.opti.debug.value(self.opti_vars.J);

                % disp('U_0:'); disp(res.U)
                % disp('J:'); disp(res.J)
            end
        end
    end
end