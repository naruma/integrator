#!/bin/csh

foreach fa (*.f90)
 cp $fa ${fa:r}.f03 
end

