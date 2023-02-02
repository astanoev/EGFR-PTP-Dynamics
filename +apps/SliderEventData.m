classdef (ConstructOnLoad) SliderEventData < event.EventData
   properties
      NewState
   end
   
   methods
      function data = SliderEventData(newState)
         data.NewState = newState;
      end
   end
end