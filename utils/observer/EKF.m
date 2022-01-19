classdef EKF < StateEstimator
    properties
        dhat
        A1
        H1
    end
    
    methods
        function [self] = get_observer(self, prob_info)
            % function which gets the specfied observer object using
            % information loaded from the prob_info structure
            self.prob_info = prob_info;
            self.Q = prob_info.Qobs;
            self.R = prob_info.Robs;
            
            nx = prob_info.nx;
            nu = prob_info.nu;
            nd = prob_info.nd;
            
            self.xhat = zeros(nx,1);
            self.dhat = zeros(nd,1);
            self.P = eye(nx+nd);
            
            import casadi.*
            % create casadi variables
            x = SX.sym('x',nx);
            u = SX.sym('u',nu);
            d = SX.sym('d',nd);
            
            % get linearized system matrices
            xdot = prob_info.f(x,u,d);
            ymeas = prob_info.h(x,d);
            ddot = d;
            self.A1 = Function('A1', {x,u,d}, {jacobian(vertcat(xdot,ddot),vertcat(x,d))});
            self.H1 = Function('H1', {x,d}, {jacobian(ymeas,vertcat(x,d))});
        end
        function [self, xhat, dhat] = update_observer(self, u, ymeas)
            % function which updates the observer in simulation
            
            % get predicted states
            x_next = full(self.prob_info.f(self.xhat, u, self.dhat));
            d_next = self.dhat;
            
            % get predicted measurement
            y_next = full(self.prob_info.h(x_next, d_next));
            
            % get predicted covariance
            A = full(self.A1(self.xhat, u, self.dhat));
            H = full(self.H1(self.xhat, self.dhat));
            phi = A;
            P_aug_next = phi * self.P * phi' + self.Q;
            
            % get kalman gain
            S_next = pinv(H * P_aug_next * H' + self.R);
            K_aug_next = P_aug_next * H' * S_next;
            
            % update state estimation and covariance
            x_aug_next = [x_next; d_next];
            x_aug_update = x_aug_next + K_aug_next * (ymeas - y_next);
            P_aug_update = (eye(self.prob_info.nx+self.prob_info.nd) - K_aug_next * H) * P_aug_next;
            
            % update observer properties
            self.xhat = x_aug_update(1:self.prob_info.nx);
            self.dhat = x_aug_update(self.prob_info.nx+1:end);
            self.P = P_aug_update;
            
            % define return values
            xhat = self.xhat;
            dhat = self.dhat;
        end
    end
end            

