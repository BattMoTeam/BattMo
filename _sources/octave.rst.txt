======================
Note on Octave Support
======================

`GNU Octave <https://octave.org/>`_ is an open source scientific programming language. The language syntax is by design
compatible with Matlab and the Octave developpers strive to maintain this compatibility (*drop-in compatibility*). Thus,
Octave can be used as a free alternative to Matlab.

In BattMo, we use only the matlab core functionalities meaning that none of the matlab toolboxes are required. We also
favor standard data structures and functions. It is easier to maintain and likely more robust. Then, it is not hard to
make BattMo compatible with Octave.

Indeed, BattMo solvers can be run with Octave and we **aim at maintaining this compatibility**. It means that scripts
that takes json input, compute the solutions and return the solution for each time step, will be supported. The generic
script :battmo:`runJsonScript` runs in Octave.

**However**, we have met some difficulties with 3D-plotting in Octave and we do not have the resources to solve
those. Octave user may therefore have to develop their own visualization support. Note that 1d-plotting should not be an
issue.
 


