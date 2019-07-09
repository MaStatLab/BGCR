
#' @title Bayesian Graphical Compositional Regression
#' @description Fit the Bayesian Graphical Compositional Regression for comparing two groups of microbiome count data.
#' @param PrJAP a value between 0 and 1. The prior joint alternative probability used for choosing \eqn{\alpha}.
#' @param sum_PrMAP a positive value less than the total number of internal nodes in the phylogenetic tree. The sum of prior marginal alternative probability used for choosing \eqn{\tau}. The default value is set to 2, recommended for 100 OTUs.
#' @param threshold a positive value controls the accuracy used for the empirical Bayes estimate of \eqn{\tau}.
#' @param kappa a non-negative value controls the chaining pattern will occur along the phylogenetic tree. The default value 0 introduces an explaning away effect.
#' @param tree an object of class ``phylo'' defined by the \code{ape} package that summarizes the phylogenetic information of the OTUs.
#' @param otu_group_1 a matrix containing the OTU samples of the first group. Each row represents an OTU, each column represents a sample.
#' @param otu_group_2 a matrix containing the OTU samples of the second group, same as \code{otu_group_1}.
#' @param X_group_1 a matrix containg the covariate of samples in the first group that need to be adjusted in the test. Each row represetns a sample, each column represents a covariate. The default setting adjusts for no covariate.
#' @param X_group_2 a matrix containg the covariate of samples in the second group, same as \code{X_group_2}.
#' @param nu a vector contains a sequence of the dispersion parameter \eqn{\nu} used for numerically computing the marginal likelihood.
#' @param sigma a positive value. The standard deviation of the normal prior on the regression coefficient.
#' @param verbose logicals. If true, print the number of internal node processed while computing the marginal likelihoods.
#' @return \code{BGCR} returns an object of class \code{"BGCR"}. A \code{"BGCR"} object is a list containing the following components:
#'  \itemize{
#'   \item \code{tree} - the phylogenetic tree. An object of class \code{"phylo"} containing the phylogenetic information of the OTUs.
#'   \item \code{PJAP} - the posterior joint alternative probability.
#'   \item \code{PMAP} - a vector containing the posterior marginal alternative probability at each internal node of the phylogenetic tree, ordered according to the label of the node in the phylogenetic tree. Note that the first half of the vector contains zeros since they represent the leaves of the tree.
#'   \item \code{PJAP_ind} - the posterior joint alternative probability returned by BCR, similar to \code{PJAP}.
#'   \item \code{PMAP_ind} - the posterior marginal alternative probabilities returned by BCR, similar to \code{PMAP}.
#'   \item \code{BF} - a vector containing the Bayes factor comparing the local hypotheses at each internal node. The orders are same as the orders in \code{PMAP}.
#'   \item \code{alpha} - the value of \eqn{\alpha} used for BGCR.
#'   \item \code{tau} - the value of \eqn{\tau} used for BGCR estimated by the empirical Bayes procedure.
#' }
#' @examples
#' library(ape)
#' exp_tree <- rtree(4) #create a phylogenetic tree with 4 leaves
#'
#' #generate the OTU table for the first group
#' otu1 <- matrix(c(2,3,4,6,3,5,5,8), nrow = 4, byrow = TRUE)
#'
#' #generate the OTU table for the second group
#' otu2 <- matrix(c(1,1,1,1,2,2,2,2), nrow = 4, byrow = TRUE)
#'
#' rownames(otu1) <- exp_tree$tip.label #name the OTUs in the OTU table
#' rownames(otu2) <- exp_tree$tip.label #name the OTUs in the OTU table
#'
#' result <- BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005,
#'                kappa = 0, tree = exp_tree, otu_group_1 = otu1, otu_group_2 = otu2)
#'
#' plot(x = result)
#'
#' @useDynLib BGCR
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics axis par plot polygon
#' @importFrom stats optimize reorder
#' @importFrom ape nodelabels
#' @export
BGCR <- function(PrJAP = 0.5,
                sum_PrMAP = "default",
                threshold = 0.005,
                kappa = 0,
                tree,
                otu_group_1,
                otu_group_2,
                X_group_1 = "default",
                X_group_2 = "default",
                nu = 10 ^ (seq(-1, 4)),
                sigma = 4,
                verbose = FALSE){


  ################################################################################################
  ####----Functions deal with the index----####

  ####

  get_child_index = function(x, map_for_child, which = 0){

    num_leaf = dim(map_for_child)[1]/2 + 1
    if(x <= num_leaf){
      print("Warning: node has no child.")
      return(NA)
    }else if(x > 2 * num_leaf - 1){
      print("Error: node index out of bound.")
      return(NA)
    }else{
      if(which == 0){
        return(map_for_child[(x - num_leaf) * 2, 2])
      }else if(which == 1){
        return(map_for_child[(x - num_leaf) * 2 - 1, 2])
      }else{
        print("Error: 'which' should be 0/1.")
        return(NULL)
      }
    }
  }


  ####

  get_parent_index = function(x, map_for_parent){

    num_leaf = dim(map_for_parent)[1]/2 + 1
    if(0 <= x && x <= num_leaf){
      return(map_for_parent[x, 1])
    }else if(num_leaf + 1 < x && x< 2 * num_leaf){
      return(map_for_parent[x - 1, 1])
    }else if(x == num_leaf + 1){
      print("Warning: root node has no parent.")
      return(NA)
    }else{
      print("Error: index out of bound.")
      return(NA)
    }
  }


  ####

  get_which_child = function(x, map_for_parent, map_for_child, root = 0){

    num_leaf = dim(map_for_parent)[1]/2 + 1
    if(x == num_leaf + 1){
      return(root)
    }else{
      parent = get_parent_index(x, map_for_parent)
      if(x == get_child_index(parent, map_for_child, which = 0)) return(0)
      else return(1)
    }
  }


  ################################################################################################
  ####----Functions deal with the data/BF----####

  ####

  get_data = function(tree, otu){

    data = matrix(0, nrow = num_leaf + num_node, ncol = dim(otu)[2])
    label = tree$tip.label
    for(i in 1:num_leaf){
      data[i, ] = otu[label[i], ]
    }

    for(i in end : start){
      child_left = get_child_index(i, map_for_child, 0)
      child_right = get_child_index(i, map_for_child, 1)
      data[i, ] = data[child_left, ] + data[child_right, ]
    }
    return(data)
  }


  ####

  get_cov = function(X_group_1 = "default", X_group_2 = "default", otu_group_1, otu_group_2){

    if(X_group_1 == "default" && X_group_2 == "default"){
      X_null = matrix(1, nrow = dim(otu_group_1)[2] + dim(otu_group_2)[2], ncol = 1)
      X_alt = cbind(X_null, c(rep(0, dim(otu_group_1)[2]), rep(1, dim(otu_group_2)[2])))
      return(list(X_null = X_null, X_alt = X_alt))
    }else if(X_group_1 == "default" && X_group_2 != "default"){
      print("Error: covariates for the two group should have same columns.")
      return(0)
    }else if(X_group_1 != "default" && X_group_2 == "default"){
      print("Error: covariates for the two group should have same columns.")
      return(0)
    }else if(dim(X_group_1)[2] != dim(X_group_2)[2]){
      print("Error: covariates for the two group should have same columns.")
      return(0)
    }else{
      X = rbind(X_group_1, X_group_2)
      X_null = cbind(rep(1, dim(X)[1]), X)
      X_alt = cbind(X_null, c(rep(0, dim(X_group_1)[1]), rep(1, dim(X_group_2)[1])))
      return(list(X_null = X_null, X_alt = X_alt))
    }
  }


  ####

  get_BF = function(otu_group_1, otu_group_2, X_group_1, X_group_2, tree, nu = 10 ^ (seq(-1, 4)), sigma = 4, verbose){

    data_group_1 = get_data(tree, otu_group_1)
    data_group_2 = get_data(tree, otu_group_2)

    BF = rep(0, num_leaf + num_node)
    cov = get_cov(X_group_1, X_group_2, otu_group_1, otu_group_2)
    X_null = cov$X_null
    X_alt = cov$X_alt

    for(i in start : end){
      left = get_child_index(i, map_for_child, 0)
      right = get_child_index(i, map_for_child, 1)
      n1 = cbind(data_group_1[left, ], data_group_1[right, ])
      n2 = cbind(data_group_2[left, ], data_group_2[right, ])

      BF[i] =  beta_binom_cov(n1, n2, X_null, X_alt, nu, sigma)
      if(verbose == TRUE){
        print(i)
      }
    }
    return(BF)
  }


  ################################################################################################
  ####----Functions deal with the prior----####

  ####

  get_prior = function(alpha, beta, kappa){

    prior = matrix(0, ncol = 2, nrow = 4)
    prior[1, 2] = exp(alpha)/(1 + exp(alpha))
    prior[2, 2] = exp(alpha + beta)/(1 + exp(alpha + beta))
    prior[3, 2] = prior[2, 2]
    prior[4, 2] = exp(alpha + beta + kappa)/(1 + exp(alpha + beta + kappa))
    prior[ , 1] = 1 - prior[ , 2]
    return(prior)
  }


  ####

  get_prior_updown = function(alpha, beta, kappa, only_prior = TRUE, root = 0){

    prior_updown = array(0, dim = c(num_node + num_leaf + 1, 2, 4))
    prior = get_prior(alpha, beta, kappa)
    clique_margin = matrix(0, nrow = num_node + num_leaf, ncol = 8)
    node_margin = matrix(0, nrow = num_node + num_leaf, ncol = 2)

    for(i in 1:num_leaf){
      node_margin[i, ] = prior[1, ]
      clique_margin[i, ] = c(rep(node_margin[i, 1], 4), rep(node_margin[i, 2], 4))/4
      node_margin[i, ] = c(1, 0)
    }

    for(i in end : start){
      left = get_child_index(i, map_for_child, 0)
      right = get_child_index(i, map_for_child, 1)
      left_margin = node_margin[left, ]
      right_margin = node_margin[right, ]

      margin = prior
      margin[1, ] = margin[1, ] * left_margin[1] * right_margin[1]
      margin[2, ] = margin[2, ] * left_margin[1] * right_margin[2]
      margin[3, ] = margin[3, ] * left_margin[2] * right_margin[1]
      margin[4, ] = margin[4, ] * left_margin[2] * right_margin[2]

      clique_margin[i, 1:4] = margin[, 1]
      clique_margin[i, 5:8] = margin[, 2]

      node_margin[i, ] = c(sum( clique_margin[i, 1:4]),  sum(clique_margin[i, 5:8]))
      prior_updown[i, 1, ] = clique_margin[i, 1:4]/node_margin[i, 1]
      prior_updown[i, 2, ] = clique_margin[i, 5:8]/node_margin[i, 2]
    }

    r = num_leaf + 1
    if(root == 0){
      prior_updown[end + 1, 1, ] = c(node_margin[r, 1], node_margin[r, 1], node_margin[r, 2], node_margin[r, 2])/2
    }else{
      prior_updown[end + 1, 1, ] = c(node_margin[r, 1], node_margin[r, 2], node_margin[r, 1], node_margin[r, 2])/2
    }
    prior_updown[end + 1, 2, ] = prior_updown[end + 1, 1, ]

    if(only_prior){
      return(prior_updown)
    }else{
      return(list(prior_updown = prior_updown, clique_margin = clique_margin, node_margin = node_margin) )
    }
  }


  ####

  compute_PrJAP = function(alpha, beta, kappa = 0, root = 0){

    res = get_prior_updown(alpha, beta, kappa, FALSE, root)
    prior_updown = res$prior_updown
    node_margin = res$node_margin
    pr_jap = node_margin[num_leaf + 1, 1]
    for(i in (num_leaf + 1):(num_leaf + num_node)){
      left = get_child_index(i, map_for_child, 0)
      right = get_child_index(i, map_for_child, 1)
      pr_jap = pr_jap * prior_updown[i, 1, 1]
    }
    pr_sum_pmap = sum(node_margin[ , 2])
    return(list(pr_jap = 1 - pr_jap, pr_sum_pmap = pr_sum_pmap))
  }


  ####

  get_alpha_beta = function(PrJAP = 0.5, sum_PrMAP = 2, threshold, max_iter = 50){

    alpha_b = -20
    alpha_t = 20
    alpha_mid = (alpha_b + alpha_t)/2
    beta = 0
    prjap_mid = compute_PrJAP(alpha_mid, beta, kappa = 0)$pr_jap
    diff = abs(prjap_mid - PrJAP)

    iter = 0
    while(iter <= max_iter && diff >= threshold){
      if(prjap_mid <= PrJAP){
        alpha_b = alpha_mid
        alpha_mid = (alpha_b + alpha_t)/2
        prjap_mid = compute_PrJAP(alpha_mid, beta, kappa = 0)$pr_jap
        diff = abs(prjap_mid - PrJAP)
      }else{
        alpha_t = alpha_mid
        alpha_mid = (alpha_b + alpha_t)/2
        prjap_mid = compute_PrJAP(alpha_mid, beta, kappa = 0)$pr_jap
        diff = abs(prjap_mid - PrJAP)
      }
      iter = iter + 1
    }

    alpha = alpha_mid
    beta_b = -20
    beta_t = 20
    beta_mid = (beta_b + beta_t)/2
    pr_sum_pmap_mid = compute_PrJAP(alpha, beta_mid, kappa = 0)$pr_sum_pmap
    diff = abs(pr_sum_pmap_mid - sum_PrMAP)

    while(iter <= max_iter && diff >= threshold){
      if(pr_sum_pmap_mid <= sum_PrMAP){
        beta_b = beta_mid
        beta_mid = (beta_b + beta_t)/2
        pr_sum_pmap_mid = compute_PrJAP(alpha, beta_mid, kappa = 0)$pr_sum_pmap
        diff = abs(pr_sum_pmap_mid - sum_PrMAP)
      }else{
        beta_t = beta_mid
        beta_mid = (beta_b + beta_t)/2
        pr_sum_pmap_mid = compute_PrJAP(alpha, beta_mid, kappa = 0)$pr_sum_pmap
        diff = abs(pr_sum_pmap_mid - sum_PrMAP)
      }
      iter = iter + 1
    }
    return(list(alpha = alpha, beta = beta_mid))
  }

  ##################################################################################################
  ####----Functions for the posterior----####

  ####

  get_xi = function(prior_updown, root = 0){

    xi = array(0, dim = c(num_node + num_leaf + 1, 8, 8))
    for(i in start : end){
      updown_node = prior_updown[i, , ]
      updown_clique = matrix(0, nrow = 8, ncol = 8)
      if(get_which_child(i, map_for_parent, map_for_child, root) == 0){
        updown_clique[c(1, 2, 5, 6), 1:4] = rep(updown_node[1, ], each = 4)
        updown_clique[c(3, 4, 7, 8), 5:8] = rep(updown_node[2, ], each = 4)
      }else{
        updown_clique[c(1, 3, 5, 7), 1:4] = rep(updown_node[1, ], each = 4)
        updown_clique[c(2, 4, 6, 8), 5:8] = rep(updown_node[2, ], each = 4)
      }
      xi[i, , ] = updown_clique
    }

    i = end + 1
    updown_node = prior_updown[i, , ]
    updown_clique = matrix(0, nrow = 8, ncol = 8)
    updown_clique[c(1, 2, 5, 6), 1:4] = rep(updown_node[1, ], each = 4)
    updown_clique[c(3, 4, 7, 8), 5:8] = rep(updown_node[2, ], each = 4)
    xi[i, , ] = updown_clique
    return(xi)
  }


  ####

  compute_phi = function(xi, BF){

    phi = matrix(1, ncol = 8, nrow = num_leaf + num_node + 1)
    for(i in end : start){
      mar_likelihood = c(rep(1, 4), rep(BF[i], 4))
      left = get_child_index(i, map_for_child, which = 0)
      right = get_child_index(i, map_for_child, which = 1)
      xi_A = xi[i, , ]
      phi[i, ] = xi_A %*% diag(phi[left, ] * phi[right, ]) %*% mar_likelihood
    }

    i = end + 1
    xi_A = xi[i, , ]
    mar_likelihood = rep(1, 8)
    left = num_leaf + 1
    phi[i, ] = xi_A %*% diag(phi[left, ] * rep(1, 8)) %*% mar_likelihood
    return(phi)
  }


  ####

  compute_posterior = function(prior_updown, BF, root = 0){

    xi = get_xi(prior_updown, root = 0)
    phi = compute_phi(xi, BF)
    xi_tilde = array(0, dim = c(num_node + num_leaf + 1, 8, 8))

    for(i in start : end){
      left = get_child_index(i, map_for_child, which = 0)
      right = get_child_index(i, map_for_child, which = 1)
      mar_likelihood = c(rep(1, 4), rep(BF[i], 4))
      xi_tilde[i, , ] = diag(1/phi[i, ]) %*% xi[i, , ] %*% diag(mar_likelihood) %*% diag(phi[left, ] * phi[right, ])
    }

    i = end + 1
    xi_tilde[i, , ] = diag(1/phi[i, ]) %*% xi[i, , ] %*% diag(phi[num_leaf + 1, ])
    post_updown = array(0, dim = c(num_node + num_leaf + 1, 2, 4))
    for(i in start : (end + 1)){
      post_updown[i, 1, ] = xi_tilde[i, 1, 1:4]
      post_updown[i, 2, ] = xi_tilde[i, 8, 5:8]
    }
    return(post_updown)
  }


  ####

  compute_PJAP = function(post_updown){

    null = 1
    for(i in start : end ){
      null = null * post_updown[i, 1, 1]
    }
    null = null * post_updown[end + 1, 1, 1] * 2
    return(1 - null)
  }


  ####

  compute_PMAP = function(post_updown){

    PMAP = rep(0, num_leaf + num_node)
    already = rep(0, num_leaf + num_node)
    PMAP[num_leaf + 1] = post_updown[end + 1, 2, 4] * 2
    already[num_leaf + 1] = 1

    for(i in (start + 1):end){
      if(already[i] == 0){
        clique_margin = matrix(0, ncol = 4, nrow = 2)
        parent = get_parent_index(i, map_for_parent)
        which = get_which_child(i, map_for_parent, map_for_child, root = 0)
        clique_margin[1, ] = post_updown[parent, 1, ] * (1 - PMAP[parent])
        clique_margin[2, ] = post_updown[parent, 2, ] * PMAP[parent]
        j = get_child_index(parent, map_for_child, 1 - which)
        if(which == 0){
          PMAP[i] = sum(clique_margin[ , 3:4])
          if(j > num_leaf) PMAP[j] = sum(clique_margin[ , c(2, 4)])
          already[i] = already[j] = 1
        }else{
          PMAP[i] = sum(clique_margin[ , c(2, 4)])
          if(j > num_leaf) PMAP[j] = sum(clique_margin[ , 3:4])
          already[i] = already[j] = 1
        }
      }
    }
    return(PMAP)
  }

  ####

  compute_ml = function(beta){

    prior_updown = get_prior_updown(alpha, beta, kappa, only_prior = TRUE, root = 0)
    xi = get_xi(prior_updown)
    phi = compute_phi(xi, BF)
    ml = phi[end + 1, 1]
    return(ml)
  }

  ##################################################################################################
  ####----necessary variables----####

  tree = reorder(tree, order = "cladewise")
  num_node = tree$Nnode
  num_leaf = tree$Nnode + 1

  start = num_leaf + 1
  end = num_leaf + num_node

  tree_postorder = reorder(tree, order = "postorder")
  map_for_child = tree_postorder$edge[dim(tree_postorder$edge)[1]:1, ]  ## map to find child index
  map_for_parent = map_for_child[order(map_for_child[,2]), ]            ## map to find parent index

  ####----get BF----####

  BF = get_BF(otu_group_1, otu_group_2, X_group_1, X_group_2, tree, nu, sigma, verbose)
  BF[BF == "NaN"] = 1

  ####----get parameters----####

  if(sum_PrMAP == "default"){
    parameter = get_alpha_beta(PrJAP, 0, threshold)
    alpha = parameter$alpha
    PrJAP_temp = compute_PrJAP(alpha, -20)
    beta_min = 0
    beta_max = get_alpha_beta(PrJAP_temp$pr_jap, 2, threshold)$beta

    beta = optimize(compute_ml, lower = beta_min, upper = beta_max, maximum  = TRUE)$maximum
    prior_updown = get_prior_updown(alpha, beta, kappa, only_prior = TRUE, root = 0)
    post_updown = compute_posterior(prior_updown, BF, root = 0)
    PJAP = compute_PJAP(post_updown)
    PMAP = compute_PMAP(post_updown)
  }else{
    parameter = get_alpha_beta(PrJAP, sum_PrMAP, threshold)
    alpha = parameter$alpha
    beta = parameter$beta
    prior_updown = get_prior_updown(alpha, beta, kappa, only_prior = TRUE, root = 0)
    post_updown = compute_posterior(prior_updown, BF, root = 0)
    PJAP = compute_PJAP(post_updown)
    PMAP = compute_PMAP(post_updown)
  }

  ####----get posterior----####

  beta_ind = 0
  prior_updown_ind = get_prior_updown(alpha, beta_ind, kappa, only_prior = TRUE, root = 0)
  post_updown_ind = compute_posterior(prior_updown_ind, BF, root = 0)
  PJAP_ind = compute_PJAP(post_updown_ind)
  PMAP_ind = compute_PMAP(post_updown_ind)

  res = list(tree = tree, PJAP = PJAP, PMAP = PMAP, PJAP_ind = PJAP_ind, PMAP_ind = PMAP_ind,
             BF = BF, alpha = alpha, tau = beta)
  class(res) = "BGCR"
  return(res)
}









##################################################################################################
##################################################################################################
####---plot the PMAP----####


#' @title Visualizing the BGCR results on the phylogenetic tree
#' @description This function plots the BGCR results on the phylogenetic tree.
#' @param x a BGCR object returned by the \code{"BGCR"} function.
#' @param BF logicals. If true, the bayes factors are plotted at each internal node.
#' @param ind logicals. If true, plot the PMAPs under BCR, the model without the graphical structure.
#' @param cex a positive value controls the size of the node label.
#' @param main a string contains the title of the plot.
#' @param legend logicals. If true, plot the legend.
#' @param ... further arguments to be passed to plot or to plot.BGCR.
#' @export
plot.BGCR <- function(x, BF = FALSE, ind = FALSE, cex = 0.3, main = "PMAP", legend = TRUE, ...){

  res = x
  tree = res$tree
  num_node = res$tree$Nnode
  num_leaf = res$tree$Nnode + 1
  start = num_leaf + 1
  end = num_leaf + num_node
  for(i in 1:num_node){
    tree$node.label[i] = i + num_leaf
  }

  #if(subtree != "whole"){
  #  tree = subtrees(tree)[[subtree]]
  #  start = min(as.numeric(tree$node.label))
  #  end = max(as.numeric(tree$node.label))
  #}

  if(class(res)!="BGCR"){
    print("ERROR: this function takes a BGCR object.")
    return(0)
  }

  if(ind == FALSE){
    if(BF == FALSE){
      par(mai=c(0.5, 0.4, 0.6, 1.1))
      col_Pal = colorRampPalette(c('white', 'red'))
      node_col = col_Pal(500)[as.numeric(cut(c(res$PMAP[start : end], 0, 1), breaks = 500)) ]
      plot(tree, show.tip.label = FALSE, use.edge.length = FALSE,show.node.label=FALSE,
           direction="downwards", cex = 0.6, main=main)
      nodelabels(text=format(round(res$PMAP[start : end], digits=5), nsmall=2), cex=cex, bg=node_col, frame="circle")

      if(legend == TRUE) {legend_col(col_Pal(500), seq(0, 1)) }
    }else{
      par(mai=c(0.5, 0.4, 0.6 , 1.1))
      pi0 = exp(log(0.5) / num_node)
      pi1 = 1 - exp(log(0.5) / num_node)
      priorp = pi1 * res$BF[start : end] / (pi1 * res$BF[start : end] + pi0)
      col_Pal = colorRampPalette(c('white', 'red'))
      node_col = col_Pal(500)[as.numeric(cut(priorp, breaks = 500)) ]
      plot(tree, show.tip.label = FALSE, use.edge.length = FALSE,show.node.label=FALSE,
           direction="downwards", cex = 0.3, main="BF")
      nodelabels(text=format(round(priorp, digits=2), nsmall=2), cex=0.3, bg=node_col, frame="circle")
      legend_col(col_Pal(500), seq(0, 1))
    }
  }else{
    if(BF == FALSE){
      par(mai=c(0.5, 0.4, 0.6, 1.1))
      col_Pal = colorRampPalette(c('white', 'red'))
      node_col = col_Pal(500)[as.numeric(cut(c(res$PMAP_ind[start : end], 0, 1), breaks = 500)) ]
      plot(tree, show.tip.label = FALSE, use.edge.length = FALSE,show.node.label=FALSE,
           direction="downwards", cex = 0.6, main=main)
      nodelabels(text=format(round(res$PMAP_ind[start : end], digits=5), nsmall=2), cex=cex, bg=node_col, frame="circle")

      if(legend == TRUE) {legend_col(col_Pal(500), seq(0, 1)) }
    }else{
      par(mai=c(0.5, 0.4, 0.6 , 1.1))
      pi0 = exp(log(0.5) / num_node)
      pi1 = 1 - exp(log(0.5) / num_node)
      priorp = pi1 * res$BF[start : end] / (pi1 * res$BF[start : end] + pi0)
      col_Pal = colorRampPalette(c('white', 'red'))
      node_col = col_Pal(500)[as.numeric(cut(priorp, breaks = 500)) ]
      plot(tree, show.tip.label = FALSE, use.edge.length = FALSE,show.node.label=FALSE,
           direction="downwards", cex = 0.3, main="BF")
      nodelabels(text=format(round(priorp, digits=2), nsmall=2), cex=0.3, bg=node_col, frame="circle")
      legend_col(col_Pal(500), seq(0, 1))
    }
  }
}


##################################################################################################
##################################################################################################

