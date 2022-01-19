classdef OCP
    % base definition of an OCP class
    properties
        prob_info
    end
    methods
        function [self] = OCP(prob_info)
            % OCP object constructor
            self.prob_info = prob_info;
        end
        
        function [self] = get_ocp(self)
            % method to build OCP
            error('An OCP must be built using CasADi')
        end
        
        function [self] = reset_initial_guesses(self)
            % method to reset the initial guesses of each variable of the
            % OCP
            error('Method has not been defined for this type of OCP!')
        end
        
        function [self] = set_parameters(self, varargin)
            % method to set values of the parameters of the OCP
            error('Method has not been defined for this type of OCP!')
        end
        
        function [self] = solve_ocp(self)
            % method to solve the OCP
            error('Method has not been defined for this type of OCP!')
        end
    end
end