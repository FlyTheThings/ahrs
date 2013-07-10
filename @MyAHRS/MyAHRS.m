classdef MyAHRS < handle
 
    %% Public properties
    properties (Access = public)
        SamplePeriod = 1/50;
        % Sensor characteristics
        % err covariance of accelerometers: users must change these values
        % depending on the environment under the vehicle is in operation
        var_az  = 0.962361;     % (0.1*g)^2
        var_ax	= 0.962361;
        var_ay  = 0.962361;
        % err covariance of magnetometer heading					      */
        var_psi = 0.014924;     % (7*d2r)^2
        % normal state inital value for flight build
        State=[1,0,0,0,0,0,0]; % four quaternions, 3 rate gyro biases
        Quaternion = [1,0,0,0];
        the = 0;
        phi = 0;
        psi = 0;
    end
    
    %% Private properties
    properties (Access = private)
        magCheck = 0;
        K_G = eye(3);
        be = zeros(3,1);
        % initialization of err, measurement, and process cov. matrices
        aP= eye(7)*1.0e-1;
        tprev = 0;
   end    
 
    %% Public methods
    methods (Access = public)
        function obj = MyAHRS(varargin)
            for i = 1:2:nargin
                if  strcmp(varargin{i}, 'SamplePeriod'), obj.SamplePeriod = varargin{i+1};
                else error('Invalid argument');
                end
            end;
        end
        function obj = Update(obj, Time, Gyroscope, Accelerometer, Magnetometer)
            g = 9.81;
            g2 = 2*g;
            Hdt = 0.5*(Time - obj.tprev);
            obj.tprev = Time;
            %Hdt = 0.5*obj.SamplePeriod;
            xs = obj.State;
   
            %assign new variables			*/
            pc = Gyroscope(1)*Hdt;
            qc = Gyroscope(2)*Hdt;
            rc = Gyroscope(3)*Hdt;

            %initialization of state transition matrix
            Fsys = eye(7); %mat_creat(7,7,UNIT_MATRIX);
            
            %state transition matrix			*/
            Fsys(1,2) = -pc; Fsys(1,3) = -qc; Fsys(1,4) = -rc;
            Fsys(2,1) =  pc; Fsys(2,3) =  rc; Fsys(2,4) = -qc;
            Fsys(3,1) =  qc; Fsys(3,2) = -rc; Fsys(3,4) =  pc;
            Fsys(4,1) =  rc; Fsys(4,2) =  qc; Fsys(4,3) = -pc;

            Fsys(1,5) = xs(2)*Hdt;  Fsys(1,6) = xs(3)*Hdt;  Fsys(1,7) = xs(4)*Hdt;
            Fsys(2,5) =-xs(1)*Hdt;  Fsys(2,6) = xs(4)*Hdt;  Fsys(2,7) =-Fsys(1,6);
            Fsys(3,5) =-Fsys(2,6);  Fsys(3,6) = Fsys(2,5);  Fsys(3,7) = Fsys(1,5);
            Fsys(4,5) = Fsys(1,6);  Fsys(4,6) =-Fsys(1,5);  Fsys(4,7) = Fsys(2,5);

            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %Extended Kalman filter: prediction step
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %propagation of quaternion using gyro measurement at a given sampling interval dt
            xs(1) = xs(1) - pc*xs(2) - qc*xs(3) - rc*xs(4);
            xs(2) = xs(2) + pc*xs(1) - qc*xs(4) + rc*xs(3);
            xs(3) = xs(3) + pc*xs(4) + qc*xs(1) - rc*xs(2);
            xs(4) = xs(4) - pc*xs(3) + qc*xs(2) + rc*xs(1);

            aQ = eye(7)*1.0e-7;
            aQ(5,5) = 1.0e-11;
            aQ(6,6) = 1.0e-11;
            aQ(7,7) = 1.0e-11;
            
            %error covriance propagation: P = Fsys*P*Fsys' + Q
            % 	mat_mul(Fsys,aP,tmp77);
            % 	mat_transmul(tmp77,Fsys,aP);
            %     for(i=0;i<7;i++) aP(i)(i) += aQ(i)(i);
            obj.aP = Fsys*obj.aP*Fsys' + aQ;

            % Pitch and Roll Update

            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %Extended Kalman filter: correction step for pitch and roll
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %nonlinear measurement equation of h_ahrs(x)
            % Estimated direction of gravity field
            h_ahrs(1)    = -g2*(xs(2)*xs(4)-xs(1)*xs(3));
            h_ahrs(2)    = -g2*(xs(1)*xs(2)+xs(3)*xs(4));
            h_ahrs(3)    = -g*(xs(1)*xs(1)-xs(2)*xs(2)-xs(3)*xs(3)+xs(4)*xs(4));

            %initialization of Jacobian matrix
            Hj = zeros(3,7); %mat_creat(3,7,ZERO_MATRIX);
            %compute Jacobian matrix of h_ahrs(x)
            Hj(1,1) = g2*xs(3); Hj(1,2) =-g2*xs(4); Hj(1,3) = g2*xs(1); Hj(1,4) = -g2*xs(2);
            Hj(2,1) = Hj(1,4);  Hj(2,2) =-Hj(1,3);  Hj(2,3) = Hj(1,2);  Hj(2,4) = -Hj(1,1);
            Hj(3,1) =-Hj(1,3);  Hj(3,2) =-Hj(1,4);  Hj(3,3) = Hj(1,1);  Hj(3,4) =  Hj(1,2);
            aR = eye(3)*obj.var_ax;
            %gain matrix aK = aP*Hj'*(Hj*aP*Hj' + aR)^-1
            aK = obj.aP*Hj'/(Hj*obj.aP*Hj' + aR);
            % 	mat_transmul(aP,Hj,tmp73);
            % 	mat_mul(Hj,tmp73,tmp33);
            % 	for(i=0;i<3;i++) tmp33(i)(i) += aR(i)(i);
            % 	mat_inv(tmp33,Rinv);
            % 	mat_mul(tmp73,Rinv,aK);

            %state update
            for k=1:7
                xs(k) = xs(k)+ aK(k,1)*(Accelerometer(1) - h_ahrs(1))...
                            +  aK(k,2)*(Accelerometer(2) - h_ahrs(2))...
                            +  aK(k,3)*(Accelerometer(3) - h_ahrs(3));
            end

            %error covariance matrix update aP = (I - aK*Hj)*aP
            obj.aP = (eye(7) -aK*Hj)*obj.aP;
            % 	mat_mul(aK,Hj,mat77);
            % 	mat_sub(Iden,mat77,tmpr);
            % 	mat_mul(tmpr,aP,tmp77);
            % 	mat_copy(tmp77,aP);

            obj.magCheck = obj.magCheck +1;
            if(obj.magCheck==5) % Heading update at 10 Hz
                obj.magCheck = 0;
                %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                % second stage kalman filter update to estimate the heading angle
                %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                % corrected body axes magnetic field strength vector
                % h_body = K_G*(h_measured-bias)
                h_b = [Magnetometer(1);Magnetometer(2);Magnetometer(3)];

                %magnetic heading correction due to roll and pitch angle
                cPHI= cos(obj.phi);
                sPHI= sin(obj.phi);
                Bxc = h_b(1,1)*cos(obj.the) + (h_b(2,1)*sPHI + h_b(3,1)*cPHI)*sin(obj.the);
                Byc = h_b(2,1)*cPHI - h_b(3,1)*sPHI;

                %scaling of quertonian,||q||^2 = 1
                norm_ahrs = 1.0/sqrt(xs(1)*xs(1)+xs(2)*xs(2)+xs(3)*xs(3)+xs(4)*xs(4));
                xs(1:4) = xs(1:4)*norm_ahrs;
                
                %Jacobian
                coeff1(1)= 2*(xs(2)*xs(3)+xs(1)*xs(4));
                coeff1(2)= 1 - 2*(xs(3)*xs(3)+xs(4)*xs(4));
                coeff1(3)= 2/(coeff1(1)*coeff1(1)+coeff1(2)*coeff1(2));

                temp(1) = coeff1(2)*coeff1(3);
                temp(2) = coeff1(1)*coeff1(3);
                
                Hpsi  = zeros(1,7);
                Hpsi(1,1) = xs(4)*temp(1);
                Hpsi(1,2) = xs(3)*temp(1);
                Hpsi(1,3) = xs(2)*temp(1)+2*xs(3)*temp(2);
                Hpsi(1,4) = xs(1)*temp(1)+2*xs(4)*temp(2);

                %gain matrix Kpsi = aP*Hpsi'*(Hpsi*aP*Hpsi' + Rpsi)^-1
                % 		mat_transmul(aP,Hpsi,tmp71);
                % 		invR = 1/(Hpsi(1)(1)*tmp71(1)(1)+Hpsi(1)(2)*tmp71(2)(1)+Hpsi(1)(3)*tmp71(3)(1)+Hpsi(1)(4)*tmp71(4)(1)+var_psi);
                Kpsi = obj.aP*Hpsi'/(Hpsi*obj.aP*Hpsi' + obj.var_psi);

                %state update
                psi = atan2(coeff1(1),coeff1(2));
                for k=1:7
                    xs(k) =xs(k)+ Kpsi(k,1)*obj.wraparound(atan2(-Byc,Bxc) - psi);
                end

                %error covariance matrix update aP = (I - Kpsi*Hpsi)*aP
                obj.aP = (eye(7) - Kpsi*Hpsi)*obj.aP;
                % 		mat_mul(Kpsi,Hpsi,mat77);
                % 		mat_sub(Iden,mat77,tmpr);
                % 		mat_mul(tmpr,aP,tmp77);
                % 		mat_copy(tmp77,aP);

            end

            %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %scaling of quertonian,||q||^2 = 1
            %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            norm_ahrs = 1.0/sqrt(xs(1)*xs(1)+xs(2)*xs(2)+xs(3)*xs(3)+xs(4)*xs(4));
            xs(1:4) = xs(1:4)*norm_ahrs;

            %obtain euler angles from quaternion
            obj.the = asin(-2*(xs(2)*xs(4)-xs(1)*xs(3)));
            obj.phi = atan2(2*(xs(1)*xs(2)+xs(3)*xs(4)),1-2*(xs(2)*xs(2)+xs(3)*xs(3)));
            obj.psi = atan2(2*(xs(2)*xs(3)+xs(1)*xs(4)),1-2*(xs(3)*xs(3)+xs(4)*xs(4)));
            obj.State = xs;
            obj.Quaternion = xs(1:4);
        end
        function [dta] = wraparound(obj, dta)
            %bound heading angle between -180 and 180
            if(dta >  pi), dta = dta - 2*pi; end
            if(dta < -pi), dta = dta + 2*pi; end
        end
    end
end