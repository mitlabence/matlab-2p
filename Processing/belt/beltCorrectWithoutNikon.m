function [belt_struct] = beltCorrectWithoutNikon(belt_struct)
belt_struct = beltCorrectArduinoArtifacts(belt_struct);
belt_struct = beltCorrectLength(belt_struct);
belt_struct  = beltAddRunningProperties(belt_struct);
belt_struct = beltSpeedToMeterPerSecond(belt_struct);
end