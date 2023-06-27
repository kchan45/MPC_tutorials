classdef OffsetFreeMPC < MPC
    % This class defines an offset-free MPC implementation for general
    % systems. This class utilizes the Opti stack interface of CasADi.
    properties
        opti;
        opti_vars;
        opti_params;
        mpc_type = 'offsetfree';
        offset_style = 'Limon';
    end
    methods
        function self = set_offset_style(self, style)
            % Setter function for offset free implementation style.
            % Accepted strings include 'Limon' (default) and 'Rawlings'
            if strcmpi(style, 'limon') || strcmpi(style,'rawlings')
                self.offset_style = style;
            else
                error('Unsupported offset-free implementation')
            end
        end
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
            nd = self.prob_info.nd;

            u_min = self.prob_info.u_min;
            u_max = self.prob_info.u_max;
            x_min = self.prob_info.x_min;
            x_max = self.prob_info.x_max;
            y_min = self.prob_info.y_min;
            y_max = self.prob_info.y_max;
            yc_min = self.prob_info.yc_min;
            yc_max = self.prob_info.yc_max;

            u_init = self.prob_info.u_init;
            x_init = self.prob_info.x_init;
            y_init = self.prob_info.y_init;

            f = self.prob_info.f;
            h = self.prob_info.h;
            r = self.prob_info.r;
            lstage = self.prob_info.lstage;
            lterm = self.prob_info.lterm;
            
            % create NLP opti object
            self.opti = Opti();
            p_opts = struct('verbose', 0, ...
                            'expand', true, ...
                            'print_time', 0); % problem options
            s_opts = struct('max_iter', 1000, ...
                            'print_level', 0, ...
                            'tol', 1e-6);     % solver options
            self.opti.solver('ipopt', p_opts, s_opts) % add the solver to the opti object

            % Initialize container lists for all states, inputs, and
            % predicted noise over horizon
            X = cell(Np+1,1);
            U = cell(Np,1);

            J = 0; % initialize cost/objective function
            
            % define parameter(s), variable(s), and problem
            
            X{1} = self.opti.parameter(nx); % initial state as a parameter
            self.opti.set_value(X{1}, zeros(nx,1))
            
            D = self.opti.parameter(nd); % offset free disturbances
            self.opti.set_value(D, zeros(nyc,1))
            
            if strcmpi(self.offset_style, 'limon')
                % limon style implementation
                Yref = self.opti.parameter(nyc); % target/reference output
                self.opti.set_value(Yref, zeros(nyc,1))

                Xss = self.opti.variable(nx);
                self.opti.subject_to(x_min <= Xss <= x_max)
                self.opti.set_initial(Xss, x_init)

                Uss = self.opti.variable(nu);
                self.opti.subject_to(u_min <= Uss <= u_max)
                self.opti.set_initial(Uss, u_init)

                Yss = self.opti.variable(ny);
                self.opti.subject_to(y_min <= Yss <= y_max)
                self.opti.set_initial(Yss, y_init)

                Ycss = self.opti.variable(nyc);
                self.opti.subject_to(yc_min <= Ycss <= yc_max)
                self.opti.set_initial(Ycss, zeros(nyc,1))

                self.opti.subject_to(Xss == f(Xss, Uss, D))
                self.opti.subject_to(Yss == h(Xss, D))
                self.opti.subject_to(Ycss == r(Yss))
                
            elseif strcmpi(self.offset_style, 'rawlings')
                % rawlings style implementation
                Xss = self.opti.parameter(nx);
                self.opti.set_value(Xss, x_init)

                Uss = self.opti.parameter(nu);
                self.opti.set_value(Uss, u_init)
                
            else
                error('Unsupported offset-free implementation')
            end

            % the loop below systematically defines the variables and
            % objective of the optimal control problem (OCP) over the
            % prediction horizon
            for k = 1:Np
                % variable and constraints for u_{k}
                U{k} = self.opti.variable(nu);
                self.opti.subject_to(u_min <= U{k} <= u_max)
                self.opti.set_initial(U{k}, u_init)

                Jstage = lstage(X{k}, U{k}, Xss, Uss);
                J = J + Jstage; % add to the cost (construction of the objective)

                % variable x_{k+1}
                X{k+1} = self.opti.variable(nx);
                self.opti.subject_to(x_min <= X{k+1} <= x_max)
                self.opti.set_initial(X{k+1}, x_init)
                
                self.opti.subject_to(y_min <= h(X{k+1},D) <= y_max)

                % dynamics constraint
                self.opti.subject_to(X{k+1} == f(X{k},U{k},D))
            end

            % terminal cost and constraints
            J_end = lterm(X{end}, Xss);
            J = J + J_end;
            
            % target penalty
            if strcmpi(self.offset_style, 'limon')
                if isfield(self.prob_info, 'target_penalty')
                    target_penalty = self.prob_info.target_penalty;
                else
                    warning('No target penalty defined. Assuming 0.')
                    target_penalty = 0;
                end
                J = J + target_penalty * sum((Yref - Ycss).^2);
            end
            
            % terminal equality constraint
            if isfield(self.prob_info, 'term_eq_cons')
                if self.prob_info.term_eq_cons
                    self.opti.subject_to(X{end} == Xss)
                end
            end

            % set to minimize objective/cost
            self.opti.minimize( J )

            % save containers of variables/parameters into a struct for
            % portability
            self.opti_vars = struct;
            self.opti_vars.U = U;
            self.opti_vars.X = X(2:end);
            self.opti_vars.J = J;

            self.opti_params  = struct;
            self.opti_params.X0 = X{1};
            self.opti_params.D = D;
            
            if strcmpi(self.offset_style, 'limon')
                self.opti_vars.Xss = Xss;
                self.opti_vars.Uss = Uss;
                self.opti_vars.Yss = Yss;
                self.opti_vars.Ycss = Ycss;
                
                self.opti_params.Yref = Yref;
            else
                self.opti_params.Xss = Xss;
                self.opti_params.Uss = Uss;
            end
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
            yc_init = self.prob_info.yc_init;

            % unpack relevant variable containers from problem creation
            U = self.opti_vars.U;
            X = self.opti_vars.X;
            if strcmpi(self.offset_style, 'limon')
                Uss = self.opti_vars.Uss;
                Xss = self.opti_vars.Xss;
                Yss = self.opti_vars.Yss;
                Ycss = self.opti_vars.Ycss;
                
                self.opti.set_initial(Uss, u_init)
                self.opti.set_initial(Xss, x_init)
                self.opti.set_initial(Yss, y_init)
                self.opti.set_initial(Ycss, yc_init)
            end

            for k = 1:Np
                self.opti.set_initial(U{k}, u_init)
                self.opti.set_initial(X{k}, x_init)
            end
        end
        
        function self = set_parameters(self, params_list, warn)
            % This method sets the values of the parameters of the
            % optimization problem to those provided as arguments to this
            % method.
            
            if isequal(class(params_list), 'struct')
                self.opti.set_value(self.opti_params.D, params_list.D)
                self.opti.set_value(self.opti_params.X0, params_list.X0)
                if strcmpi(self.offset_style, 'limon')
                    self.opti.set_value(self.opti_params.Yref, params_list.Yref)
                else
                    self.opti.set_value(self.opti_params.Xss, params_list.Xss)
                    self.opti.set_value(self.opti_params.Uss, params_list.Uss)
                end
            elseif isequal(class(params_list), 'cell')
                if warn
                    warning('You passed in a cell array. Make sure the cell array order is the same as the following: ')
                    disp(fieldnames(self.opti_params))
                end
                self.opti.set_value(self.opti_params.X0, params_list{1})
                self.opti.set_value(self.opti_params.D, params_list{2})
                if strcmpi(self.offset_style, 'limon')
                    self.opti.set_value(self.opti_params.Yref, params_list{3})
                else
                    self.opti.set_value(self.opti_params.Xss, params_list{3})
                    self.opti.set_value(self.opti_params.Uss, params_list{4})
                end
                
            else
                warning('Unsupported params_list class... Parameters were not set!')
            end            
        end
        
        function [self, res, feas] = solve_mpc(self)
            % This method solves the MPC as defined by the get_mpc() method
            % of this class. This method can only be called after the the
            % get_mpc() method has been called (i.e., the optimization
            % problem must be defined before it can be solved).
            
            % get fields of the opti_vars struct
            fldnames = fieldnames(self.opti_vars);
            
            % unpack relevant information from problem creation
            u_min = self.prob_info.u_min;
            u_max = self.prob_info.u_max;

            feas = 1;
            res = struct;
            try
                sol = self.opti.solve();
                for i = 1:length(fldnames)
                    if strcmpi(fldnames{i}, 'U')
                        res.(fldnames{i}) = sol.value(self.opti_vars.(fldnames{i}){1});
                    elseif isa(self.opti_vars.(fldnames{i}),'cell')
                        res.(fldnames{i}) = cell2mat(cellfun(@sol.value,self.opti_vars.(fldnames{i}),'UniformOutput',false)');
                    else
                        res.(fldnames{i}) = sol.value(self.opti_vars.(fldnames{i}));
                    end
                end

                if self.prob_info.warm_start
                    self.opti.set_initial(sol.value_variables())
                end

            catch e
                % if solve fails, get the last value
                feas = 0;
                disp(e.message)
                
                for i = 1:length(fldnames)
                    if strcmpi(fldnames{i}, 'U')
                        u = self.opti.debug.value(self.opti_vars.(fldnames{i}){1});
                        res.(fldnames{i}) = max(min(u,u_max),u_min);
                    elseif isa(self.opti_vars.(fldnames{i}),'cell')
                        res.(fldnames{i}) = cell2mat(cellfun(@(x)self.opti.debug.value(x),self.opti_vars.(fldnames{i}),'UniformOutput',false)');
                    else
                        res.(fldnames{i}) = self.opti.debug.value(self.opti_vars.(fldnames{i}));
                    end
                end

                % disp('U_0:'); disp(res.U)
                % disp('J:'); disp(res.J)
            end
        end
    end
end