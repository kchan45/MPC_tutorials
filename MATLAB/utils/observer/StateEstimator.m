classdef StateEstimator
    properties
        prob_info
        P
        Q
        R
        xhat
    end
    
    methods
        function [self] = get_observer(self, prob_info)
            % function which gets the specfied observer object using
            % information loaded from the prob_info structure
        end
        function [self] = update_observer(self, u, ymeas)
            % function which updates the observer in simulation
        end
    end
end