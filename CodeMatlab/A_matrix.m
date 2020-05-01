function A = A_matrix(steeringAngle, pedalSpeed,internalStateIn,v_linear,theta)
B=0.8+0.1*0.8*randn(1,1);
       A=[v_linear*cos(theta)
end