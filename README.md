SIES
====

A Matlab library for shape identification in electro-sensing.

This library deals with the problem of shape identification in some inverse problems with 
electric or acoustic measurements. It allows to reproduce some numerical results of the following papers:

1. Target identification using dictionary matching of generalized polarization tensors, in Foundations of Computational Mathematics, 2014
2. Tracking of a Mobile Target Using Generalized Polarization Tensors, in SIAM Journal on Imaging Sciences, 2013
3. Shape recognition and classification in electro-sensing, Proceedings of the National Academy of Sciences of the United States of America, 2013 
4. Wavelet methods for shape perception in electro-sensing, Mathematics of Computation, 2013 
5. Shape identification and classification in echolocation, Mathematics of Computation, 2013 


The library is written in Matlab object-oriented programming language and is organized in the following manner:

* +acq: classes of acquisition system
* +asymp: asymptotic expansions, CGPT, Scattering coefficients etc.
* +dico: tools for shape identification in a dictionary
* +ops: classes of layer potentials
* +PDE: classes of PDE modeling for the simulation and solution of inverse problems
* +shape: classes for C2 smooth domain
* +tools: functions of utility
* +tracking: functions for tracking of a mobile target
* +wavelet: classes of wavelet representation for shape perception
* examples: demos scripts

For usage of the package, go to the path examples/ and run the demo scripts. Your suggestions and remarks are welcome.
