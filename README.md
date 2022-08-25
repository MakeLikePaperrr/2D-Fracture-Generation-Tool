# 2D-Fracture-Generation-Tool
A MATLAB-based tool used for generating 2D complex and geologically significant fracture networks. Originally developed as part of the master thesis of Andrea Sartori (http://resolver.tudelft.nl/uuid:d20ab3d6-a63d-41b2-b74d-198a3f3f44c5) and later extended by Stephan de Hoop. Several input parameters can be specified, particularly:
- number of fractures;
- max. fractures per set;
- min. length small fractures;
- min. length large fractures;
- angle of small fractures;
- angle of large fractures.
Furthermore, an upper and lower bound on the connectivity can be specified. The connectivity affects the total number of fractures, and it can be fully honored by connecting disconnected clusters or disconnecting connected clusters of fractures such that connectivity agrees with the sampled connectivity from 

See https://doi.org/10.2118/203968-MS for more information on the generation and application of the framework.

Two test cases are provided, one for a single fracture network and one for an ensemble of fracture networks.
