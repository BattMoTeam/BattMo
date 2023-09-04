classdef SweepParameter
   properties
       values
       parameter_name
   end
   methods
       function param = SweepParameter(parameter_name,values)
            param.parameter_name=parameter_name;
            param.values = values;
       end
   end
end