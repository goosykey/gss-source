clear; clc;

s1 = OTT_MCSequence;

s1.addSegment('type','iburn');
s1.addSegment('type','prop','linecolor','c');
s1.stopConditions{4}{1} = 'apogee';

s1.addSegment('type','fburn','linecolor','y');
s1.stopConditions{5}{2} = 20000;

s1.run;
s1.plotTrajectory;