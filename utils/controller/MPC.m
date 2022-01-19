classdef MPC
    % MPC is a super class designed to be a template for particular
    % implementations of model predictive controllers (MPCs). Users should
    % develop their own MPCs by using the general structure/methods
    % provided below. Upon initialization of this class or any of its child
    % classes, users should provide a MATLAB struct that contains all of
    % the relevant problem information. This class is designed to be used
    % with CasADi to generate the optimization problems associated with
    % MPC. Users are referred to the CasADi documentation for more
    % information on CasADi.
    properties
        prob_info
    end
    methods
        function self = MPC(prob_info)
            % constructor for MPC object
            self.prob_info = prob_info;
        end
        function self = get_mpc(self, varargin)
            % This method should generate the MPC problem by unpacking
            % relevant information from the prob_info dict defined upon
            % instantiation of the class. Upon defining the optimization
            % problem, this method should save and/or return the
            % appropriate objects such that this object may be called upon
            % later. (e.g. if using the Opti stack interface of CasADi,
            % save/return the objects that reference the optimization
            % object (typically named opti) and the optimization variable
            % references)
        end
        function self = reset_initial_guesses(self, varargin)
            % This method should reset any initial guesses of the decision
            % variables passed into the optimization problem. This method
            % provides a way to simulate repeated solves of the
            % optimization problem without re-defining an entirely new
            % problem. If using the Opti stack interface of CasADi, this
            % method mainly involves using the set_initial() method of the
            % Opti object.
        end
        function self = set_parameters(self, varargin)
            % This method should (re)set any parameters in the optimization
            % problem. This method provides a way to simulated consistent
            % and repeated solves of the optimization problem without
            % re-defining an entirely new problem. If using the Opti stack
            % interface of CasADi, this method mainly involves using the
            % set_value() method of the Opti object.
        end
        function self = solve_mpc(self, varargin)
            % This method should solve the optimization problem and
            % return/save the relavent optimal variables. For MPC, this is
            % typically the first optimal input determined by the solver.
            % Users may also wish to return other values of the
            % optimization problem and/or the entire solution of the
            % optimization problem. This method should also handle any
            % Exceptions that may occur upon a call to solve the
            % optimization problem in the form of a try/except clause. If
            % using the Opti stack interface of CasADi, this method mainly
            % involves the call to the solve() method of the Opti object,
            % as well as calls to the value() method of OptiSolution/Opti
            % objects.
        end
    end
end