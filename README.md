# Tabu-Search Based Pilot Assignment for Cell-Free Massive MIMO Systems

This is a code package is related to the following scientific article:


H. Liu, J. Zhang, X. Zhang, A. Kurniawan, T. Juhana, and B. Ai, ''Tabu-Search-Based pilot assignment for cell-free massive MIMO systems,'' IEEE Trans. Veh. Technol., vol. 69, no. 2, pp. 2286–2290, Feb. 2020

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. *We encourage you to also perform reproducible research!*


## Abstract of Article

In this correspondence, an efficient pilot assignment scheme based on tabu search (TS) is proposed to mitigate the severe pilot contamination effect caused by pilot reuse in cell-free massive multiple-input multiple-output systems. By exploiting the local neighborhood search, we first define the neighborhood and search objective function. Then, the tabu list is introduced to record the prohibited solutions and to avoid being trapped in the local optimum. In addition, we consider a fair policy for the pilot reused times to further improve the system performance. We propose a TS-based pilot assignment scheme iteration to intelligently search the suboptimal pilot assignment result with low complexity. Numerical results verify that the proposed TS based pilot assignment scheme can significantly improve the 95%-likely per-user net throughput with similar complexity compared to the conventional random and greedy schemes.

## Content of Code Package

The package generates the simulation SE results which are used in Figure 2, Figure 3, Figure 4, and Figure 5. To be specific:

- `simulation_main`: Main function;
- `functiontabu_new`: Generate the tabu search based pilot assignment result;

See each file for further documentation.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
