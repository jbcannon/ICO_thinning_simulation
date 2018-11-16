#############################################################################################
#### Script for simulating stand-scale heterogeneous restoration treatments
#### Author: Jeffery Cannon (Jeffery.Cannon@colostate.edu)
#### Based on Python script from Yvette Dickinson
#### Institution: Colorado Forest Restoration Institute, Colorado State University
#### Date Created: 02/23/2018
#### Last Modified: 10/24/2018
#############################################################################################
# The goal of this script is to simulate a heterogeneous restoration treatment from a stem
# map. Levels of 'clumpiness' can be define based on proportions of trees in various group
# sizes (Tinkham et al. RMRS-GTR 365). The original code was based on a Phython script
# authored by Y. Dickinson, but updated to increase speed and flexibility
#############################################################################################

#######################################START SET UP##########################################
#---> Load libraries
packages <- c('sp', 'plyr', 'rgdal', 'rgeos')
for(package in packages){
  if(suppressMessages(!require(package,character.only=T))){
    install.packages(package,repos='https://cran.mtu.edu/')
    suppressMessages(library(package,character.only=T))
  }
}

#---> Load data (cite sources)
# Stem map data
all_trees = read.csv('dummy_trees.csv')
plot_area_ha = 4  # plot size in ha

#---> Select level of clumping, must match 'clumping' column in clumping_targets.csv
clump_level = 'High'  # 'Low', 'Mod', or 'High'

#--> Set Target basal area
target_ba_m2ha =  9.18  # Target basal area (m2/ha)

#---> Clumping targets (table based on Tinkham et al. GTR-365)
master_targets = read.csv('clumping_targets.csv')
group_size_descriptors = c("unassigned",
                           "single",
                           "2-4 trees",
                           "5-9 trees",
                           "10-15 trees",
                           "16-25 trees")

#---> Reclassification table (matches, group sizes to group names in clumping_targets.csv)
rcl = read.csv('group_size_reclass_table.csv')

#---> Number of simulations to run
run_sims = 2

#---> Group seed probability weight (tree dbh ^ x) 
w = 1.6  # preferentially anchors groups on larger trees; replace with a 1 for equal weighting.

#--> Set intertree distance (i.e., 2 x crown radius, sensu Churchill et al. 2013)
intertree_dist = 6

#---> Plot analysis live to visualize progress? (Printing output slows analysis by ~50%)
plot_progress = FALSE

#---> Load custom functions with explanations
#########################
# FUNCTION: pickGroups()
#
# This function selects a group of trees within a given size range from a list of unassigned
# trees. The inputs to the function include an unassigned tree list (i.e., stem map subset to
# unassigned trees), a distance matrix of all trees, a binary distance matrix (i.e., based on
# treshold intertree distance). In addition, the function requests a group size name (grp_sz)
# which refers 9to the gap size range in a table of tree sizes and proportions (targets).
# This code works by selecting a random (seed) tree, adding neighbors to the group, and
# expanding outward until the correct number of trees are included. If groups cannot be made
# large enough, a second seed tree is selected. If groups are too large, the outermost trees
# from the seed tree are removed. Concepts based on Churchill et al. 2013 and Tinkham et al.
# RMRS-GTR-365
#########################

pickGroup = function(unassigned_trees, dist_mat, dist_mat_b, targets){
  
  # Largest target group size currently needed
  target_group_size = subset(targets, need_more == TRUE)
  target_group_size = target_group_size[nrow(target_group_size),c('grp_min', 'grp_max')]
  target_group_size = as.numeric(target_group_size)
  target_group_size = target_group_size[1]:target_group_size[2]
  
  #Pick random tree to start group
  seed_tree = sample(unassigned_trees$tree_id, size = 1, prob = unassigned_trees$weight)
  curr_group = seed_tree
  
  # Use loop to grow trees until minimum target group size is reached, or growth stops
  while(!length(curr_group) > min(target_group_size)) {
    starting_number = length(curr_group)
    #add all intersecting trees to group
    tmp_dist = dist_mat_b[curr_group,]
    if(class(tmp_dist) == 'logical') {curr_group = which(tmp_dist)} else {
      curr_group = which(apply(tmp_dist, MARGIN = 2, FUN = sum) > 0)
    }
    #catch if trees aren't being added.
    if(starting_number == length(curr_group)) {break}
  }
  
  # Update target group size to largest possible size given group
  target_group_size = subset(targets, need_more == TRUE & grp_min <= length(curr_group))
  target_group_size = tail(target_group_size, 1)
  target_group_size = with(target_group_size, grp_min:grp_max)
  
  #remove excess trees if group size was exceeded
  if(length(curr_group) > max(target_group_size)){
    exact_target_size = sample(target_group_size, size = 1)
    tmp = d[seed_tree, curr_group]
    tmp = sort(tmp)
    tmp = head(tmp, exact_target_size)
    tmp = as.numeric(names(tmp))
    curr_group = tmp
  }
  curr_group = as.numeric(curr_group)
  return(curr_group)
}

#########################
# FUNCTION: cutGroup()
#
# This function identifies all trees that are 'touching' (i.e., within intertree distance of
# a group of trees selected by pickGroup(). The function takes a list of tree ids forming a group
# a stem map of unassigned trees (i.e., stem map subset to unassigned trees), and a binary
# distance matrix of all trees. The output is a list of all trees within an intertree distance
# of the input trees.
#########################

cutGroup = function(tree_group, unassigned_trees, dist_mat_b){
  tmp_dist = dist_mat_b[tree_group,]
  if(class(tmp_dist) == 'logical') {cut_group = which(tmp_dist)} else {
    cut_group = which(apply(tmp_dist, MARGIN = 2, FUN = sum) > 0)
  }
  cut_group = cut_group[! cut_group %in% tree_group]
  cut_group = cut_group[cut_group %in% unassigned_trees$tree_id]
  
  cut_group = as.numeric(cut_group)
  return(cut_group)
}

#########################
# FUNCTION: updateTargets()
#
# This function takes a table containing target proportions of trees in each group size class
# (targets) and an updated tree list containing group membership (trees) and outputs a new
# targets table updated with current distribution of tree size membership. 
#########################

updateTargets = function(targets, trees){
  for(i in levels(trees$grp_nm)){
    targs = ddply(trees, .(grp_nm), summarize, curr_BA = sum(ba_m2), .drop = FALSE)
    targs = tail(targs, -1)
    targets$curr_BA <- NULL
    targets = merge(targets, targs)
    targets$need_more = with(targets, curr_BA < BA)
    return(targets)
  }
}
########################################END SET UP###########################################
  
######################################START ANALYSIS#########################################
#---> Clean up targets table for run
master_targets = subset(master_targets, clumping == clump_level)
master_targets$clumping = NULL
master_targets$grp_nm = as.character(master_targets$grp_nm)

#---> Prepare master tree list
all_trees$ba_m2 = with(all_trees, DBH_cm/200)^2*pi
all_trees$grp_id = 0
all_trees$grp_size = 0
all_trees$grp_nm = "u"
all_trees$grp_nm <- as.factor(all_trees$grp_nm)
levels(all_trees$grp_nm) <- letters[c(21, 1:5)]
all_trees$tree_id <- 1:length(all_trees$ba_m2)

#---> Calculate distance matrix from tree list
master_d = dist(cbind(all_trees$X, all_trees$Y), diag = TRUE)
master_d = as.matrix(master_d)
master_b = master_d <= intertree_dist

#--> Set weights based on dbh
all_trees$weight =  (with(all_trees, DBH_cm)^w) #weight based on dbh

for(run in 1:run_sims){ #include the numbers you'd like to run in a loop... ie, 1:10, change file name base below to (hi, low, med, etc.)

  #---> Set outputs for individual runs
  run = as.character(sprintf('%02d', run))
  output_stem_map = paste0('Thin_Sim', clump_level, run, '.csv')
  plot_title = paste(clump_level, 'clumping')
  trees <- all_trees # create neew tree object so simulation can reset.
  d = master_d # Reset distance matrices
  b = master_b # 
  
  #---> Calculate target basal area for each tree group size
  targets = master_targets
  targets$BA = targets$p * target_ba_m2ha * plot_area_ha
  targets$curr_BA = 0
  targets$need_more = targets$curr_BA < targets$BA
  total_BA = 0
  
  #---> Show plot to display progress
  if(plot_progress) plot(Y ~ X, data = trees, asp = TRUE, pch = 1, cex = 0.1)
  i = 0
  
  #---> Loop through all group sizes (letters 'a' - 'e', and identify groups until targets are met)
  while(any(targets$need_more) == TRUE & total_BA <= target_ba_m2ha){
      i = i + 1
      unassigned_trees = subset(trees, grp_id == 0)
      
      # If there are enough trees, select next group
      if(dim(unassigned_trees)[1] > 1){
        tmp_grp = pickGroup(unassigned_trees, d, b, targets)
        # Find group and assign group membership
        tmp_grp = pickGroup(unassigned_trees, d, b, targets)
        
        lett = rcl[length(tmp_grp),'grp_nm']
        trees[trees$tree_id %in% tmp_grp, "grp_id"] <- paste(lett, sprintf('%04d', i), sep = '')
        trees[trees$tree_id %in% tmp_grp, "grp_size"] <- length(tmp_grp)
        trees[trees$tree_id %in% tmp_grp, "grp_nm"] <- lett
        
        # Cut trees surrounding/touching group
        unassigned_trees = subset(trees, grp_id == 0)
        cutTrees = cutGroup(tmp_grp, unassigned_trees, b)
        trees[trees$tree_id %in% cutTrees, "grp_id"] <- -1
        trees[trees$tree_id %in% cutTrees, "grp_size"] <- -1
        trees[trees$tree_id %in% cutTrees, "grp_nm"] <- NA
        
        # Remove cut trees from consideration by updating distance matrix
        b[cutTrees, ] <- FALSE
        b[, cutTrees] <- FALSE
        
      } else {
        # If there are not enough trees, randomly choose a 'cut' tree to add back in
        add_tree = with(subset(trees, grp_size == -1), sample(x = tree_id, size = 1, prob = weight))
        trees[trees$tree_id == add_tree,]$grp_id = 'addback'
        trees[trees$tree_id == add_tree,]$grp_size = 1
        trees[trees$tree_id == add_tree,]$grp_nm = 'a'
      }
      
      # update targets
      targets = updateTargets(targets, trees)
      total_BA = sum(targets$curr_BA) / plot_area_ha
      
      # Print progress to console and plot device
      if(plot_progress) points(Y ~ X, data = trees, subset = tree_id %in% cutTrees, pch = 16, cex = 0.1, col = 'grey')
      if(plot_progress) points(Y ~ X, data = trees, subset = tree_id %in% tmp_grp, pch = 16, cex = 0.5, col = sample(rainbow(10), size = 1))
      prog = round(sum(targets$curr_BA)/sum(targets$BA), 3)*100
      if(plot_progress) message('Simulation ', run, ' --- ', prog, '% Complete. Group of size ', length(tmp_grp), ' created')
  }
  
  #---> Change level names (letters) to descriptive titles
  levels(trees$grp_nm) <- group_size_descriptors
  
  #---> Subset data to remove all unassigned and cut trees
  trees = subset(trees, grp_nm != 'unassigned')
  #######################################END ANALYSIS##########################################
  
  ######################################START OUTPUTS##########################################
  # Write thinned stem map (tabular) to file
  write.csv(trees, output_stem_map, row.names = FALSE)
  
  # Print final basal area
  BA = sum(targets$curr_BA)/plot_area_ha
  BA = round(BA,2)
  message(paste('Final BA:', BA))
  #######################################END OUTPUTS###########################################
}

