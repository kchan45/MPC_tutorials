classdef MinimumTimeOCP < OCP
    properties
        opti
        opti_vars
        opti_params
    end
    methods
        function [self] = get_ocp(self, N)
            % create the OCP for the minimum time problem for the plasma
            % system
            nx = self.prob_info.nx;
            nu = self.prob_info.nu;
            f = self.prob_info.f;
            Tmax = self.prob_info.Tmax;
            u_min = self.prob_info.u_min;
            u_max = self.prob_info.u_max;
            x0 = self.prob_info.x0;
            xN = self.prob_info.xN;
            x_init = self.prob_info.x_init;
            
            % requires casadi
            self.opti = casadi.Opti();
            X = self.opti.variable(nx,N+1);  % state trajectory
            U = self.opti.variable(nu,N);    % input trajectory
            tf = self.opti.variable();       % final time
            
            X0 = self.opti.parameter(nx);    % initial state parameter
            self.opti.set_value(X0, x0);
            XN = self.opti.parameter(nx);    % terminal state parameter
            self.opti.set_value(XN, xN);
            

            self.opti.minimize(tf);      % objective; treatment time
            
            % formulate discrete dynamics using fixed step 
            dt = tf/N; % length of a control interval
            for k = 1:N % loop over control intervals
               % Runge-Kutta 4 integration
               k1 = f(X(:,k),         U(:,k));
               k2 = f(X(:,k)+dt/2*k1, U(:,k));
               k3 = f(X(:,k)+dt/2*k2, U(:,k));
               k4 = f(X(:,k)+dt*k3,   U(:,k));
               x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4); 
               self.opti.subject_to(X(:,k+1)==x_next); % close the gaps
            end
            
            % ---- path constraints -----------
            self.opti.subject_to(X(1,:) <= Tmax);     % path constraint
            self.opti.subject_to(u_min <= U <= u_max);        % control is limited

            % ---- boundary conditions --------
            self.opti.subject_to(X(:,1) == X0);       % start at initial states 
%             self.opti.subject_to(abs(X(:,N+1) - XN) <= 1e-3);     % finish at desired terminal states (within tolerance)
            self.opti.subject_to(X(:,N+1) == XN);

            % ---- misc. constraints  ----------
            self.opti.subject_to(tf>=0); % Time must be positive

            % ---- initial values for solver ---
            self.opti.set_initial(X(1,:), x_init(1));
            self.opti.set_initial(tf, 100);

            % ---- solve NLP              ------
            p_opts = struct('verbose', 0, ...
                            'print_time', 0); % problem options
            s_opts = struct('print_level', 0, ...
                            'max_iter', 500);     % solver options
            self.opti.solver('ipopt', p_opts, s_opts); % set numerical backend
            
            self.opti_vars = struct();
            self.opti_vars.X = X;
            self.opti_vars.U = U;
            self.opti_vars.tf = tf;
            
            self.opti_params = struct();
            self.opti_params.X0 = X0;
            self.opti_params.XN = XN;
        end
        
        function [self] = reset_initial_guesses(self)
            % reset the initial guesses of the variables (for
            % reproducibility between runs if re-using object)
            x_init = self.prob_info.x_init;
            
            self.opti.set_initial(self.opti_vars.X(1), x_init(1));
            self.opti.set_initial(self.opti_vars.tf, 100);
        end
        
        function [self] = set_parameters(self, params_list, warn)
            % set the parameters of the problem
            if isequal(class(params_list), 'struct')
                self.opti.set_value(self.opti_params.X0, params_list.X0)
                self.opti.set_value(self.opti_params.XN, params_list.XN)
            elseif isequal(class(params_list), 'cell')
                if warn
                    warning('You passed in a cell array. Make sure the cell array order is the same as the following: ')
                    disp(fieldnames(self.opti_params))
                end
                self.opti.set_value(self.opti_params.X0, params_list{1})
                self.opti.set_value(self.opti_params.XN, params_list{2})
            else
                warning('Unsupported params_list class... Parameters were not set!')
            end            
        end
        
        function [self, res, feas, sol] = solve_ocp(self)
            % solve the OCP and extract the solution
            % unpack relevant information from problem creation
            u_min = self.prob_info.u_min;
            u_max = self.prob_info.u_max;

            feas = 1;
            res = struct();
            try
                sol = self.opti.solve();
                res.U = sol.value(self.opti_vars.U(1));
                res.J = sol.value(self.opti_vars.tf);

                if self.prob_info.warm_start
                    self.opti.set_initial(sol.value_variables())
                end

            catch
                % if solve fails, get the last value
                feas = 0;
                sol = [];
                
                u = self.opti.debug.value(self.opti_vars.U(1));
                res.U = max(min(u,u_max),u_min);
                res.J = self.opti.debug.value(self.opti_vars.tf);

                % disp('U_0:'); disp(res.U)
                % disp('J:'); disp(res.J)
            end
        end
    end
end