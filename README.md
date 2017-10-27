Receive Combining vs. Multi-Stream Multiplexing in Downlink Systems with Multi-Antenna Users
==================

This is a code package is related to the follow scientific article:

Emil Björnson, Marios Kountouris, Mats Bengtsson, Björn Ottersten, “[Receive Combining vs. Multi-Stream Multiplexing in Downlink Systems with Multi-Antenna Users](http://arxiv.org/pdf/1207.2776),” IEEE Transactions on Signal Processing, vol. 61, no. 13, pp. 3431-3446, July 2013.

The package contains a simulation environment, based on Matlab, that reproduces all the numerical results and figures in the article. *We encourage you to also perform reproducible research!*


## Abstract of Article

In downlink multi-antenna systems with many users, the multiplexing gain is strictly limited by the number of transmit antennas N and the use of these antennas. Assuming that the total number of receive antennas at the multi-antenna users is much larger than N, the maximal multiplexing gain can be achieved with many different transmission/reception strategies. For example, the excess number of receive antennas can be utilized to schedule users with effective channels that are near-orthogonal, for multi-stream multiplexing to users with well-conditioned channels, and/or to enable interference-aware receive combining. In this paper, we try to answer the question if the N data streams should be divided among few users (many streams per user) or many users (few streams per user, enabling receive combining). Analytic results are derived to show how user selection, spatial correlation, heterogeneous user conditions, and imperfect channel acquisition (quantization or estimation errors) affect the performance when sending the maximal number of streams or one stream per scheduled user—the two extremes in data stream allocation.

While contradicting observations on this topic have been reported in prior works, we show that selecting many users and allocating one stream per user (i.e., exploiting receive combining) is the best candidate under realistic conditions. This is explained by the provably stronger resilience towards spatial correlation and the larger benefit from multi-user diversity. This fundamental result has positive implications for the design of downlink systems as it reduces the hardware requirements at the user devices and simplifies the throughput optimization.


## Content of Code Package

The article contains 9 simulation figures. The majority of the figures are generated by the Matlab script simulationFigure4and5and6and7and8and9.m by changing a parameter. The remaining figures are generated by simulationFigure3.m, simulationFigure10.m, and simulationFigure11.m.

The package contains 8 additional scripts with new Matlab functions: functionBlockDiagonalization.m, functionGreedyStreamAllocation.m, functionGreedyZeroForcing.m, functionRVQfeedback_BD.m, functionRVQfeedback_ZF_MESC.m, functionRVQfeedback_ZF_QBC.m, functionSumrateComputation.m, and functionUserSelectionCBSUS.m. These functions are called by the Matlab scripts.

See each file for further documentation. 


## Acknowledgements

The research leading to these results has received funding from the European Research Council under the European Community's Seventh Framework Programme (FP7/2007-2013)/ERC grant agreement number 228044. The work of E. Björnson is funded by the International Postdoc Grant 2012-228 from The Swedish Research Council.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
