
% paramteres for nesterov accelerated gradient descent
classdef Nesterov
   properties
      lambda
      gamma
   end
methods
   function obj = Nesterov
        obj.lambda = [0 0];
        obj.gamma = (1-obj.lambda(1))/obj.lambda(2);    
   end
   
   function obj = nextGamma(obj)
        obj.lambda(1) = obj.lambda(2);
        obj.lambda(2) = (sqrt(4*(obj.lambda(1)^2)+1)+1)/2;
       
        obj.gamma = (1-obj.lambda(1))/obj.lambda(2);    

   end
       
   end
end