# ICO_thinning_simulation
Algorithm to simulate heterogenous restoration treatments on stem maps to target basal area and target tree group size distributions

Summary: This simulation software uses a tree stand map as input, and simulates restoration treatments with goals to thin the forest to a given basal area (m2 ha-1) while maintaining specified group sizes. 


Description:
In many dry-conifer forests of the western U.S., restoration objectives include thinning forests to reduce potential for wildfire, while also restoring elements of historical forest structure that once characterized dry conifer forests such as large openings, and variability in tree group size. This simulation software is designed to simulate such restoraiton treatments using a stem map. Using a stem map as input data, the algorithm thins trees until a target basal area and group size distribution is acheived. For example, The simulation can be used to simulate a restoration treatment where basal area is reduced to 30 m2 ha-1 while maintaining 50% of trees as individuals, 40% in groups of 1-2, and 10% in groups of 3-5. Basal area targets as well as group size distributions and size ranges can be specified. Although groups are randomly seeded in the simulation, more realistic group seeding can be acheived by adjusting the probability weights for tree selection. Default settings, weight trees according to dbh^1.6.

Citation:
Cannon, JB, 2018. Heterogeneous restoration treatment simulator. Github.
