# launcher_trajectory
Repository for a space launch vehicle trajectory simulation

This program is not yet where I want it to be, but it's working perfectly fine.
It's is simulating the ascent of a multistage rocket. The technical data must be stored in 
an excel spreadsheet that is in the form of the Saturn5.xls example. The output is a set of graphs that show
altitude vs time, altitude vs ground range, acceleration vs time, velocity vs time, dynamic pressure vs time and rocket mass vs time.

However, at the moment the rocket is ascending completely vertically and will therefore not enter an orbit.
That will be fixed by implementing a function to "steer" the rocket, either pre-programmed as an input in the excel spreadsheet or
as a realtime interactive input.
