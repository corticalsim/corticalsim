# Introduction
Plant morphogenesis is strongly dependent on the directional growth and the subsequent oriented division of individual cells.
It has been shown that the plant cortical microtubule array plays a key role in controlling both of these processes.
This ordered structure emerges as the collective result of stochastic interactions between large numbers of dynamic microtubules.
A number of analytical and computational approaches to studying the dynamics of cortical microtubules have been proposed in order to elucidate this complex self-organization process.
To date, however, these models have been restricted to two-dimensional planes or geometrically simple surfaces in three dimensions, which strongly limits their applicability as plant cells display a wide variety of shapes.
This limitation is even more acute, as both local as well as global geometrical features of cells are expected to influence the overall organization of the array.
CorticalSim is a framework for efficiently simulating microtubule dynamics on triangulated approximations of arbitrary three-dimensional surfaces.
This allows the study of microtubule array organization on realistic cell surfaces obtained by segmentation of microscopic images.


```{figure} assets/cs-overview.png
:alt: Overview of the CorticalSim software

Overview of the CorticalSim software, from https://doi.org/10.1371/journal.pcbi.1005959
```

For more information about CorticalSim and its use in various applications see the following publications:
- for general use http://dx.doi.org/10.3389/fphy.2014.00019
- for microtubule dependent nucleation http://dx.doi.org/10.1371/journal.pcbi.1013282
- for microtubule dependent nucleation - old algorithm only (nuc_ellipse) http://dx.doi.org/10.1088/1478-3975/8/5/056002
- for microtubule severing http://dx.doi.org/10.1073/pnas.1702650114
- for microtubule deflections https://doi.org/10.1017/qpb.2024.17
- for microtubules on triangulated geometries https://doi.org/10.1371/journal.pcbi.1005959
- for microtubule based cell division orientation https://doi.org/10.1016/j.cub.2018.07.025

