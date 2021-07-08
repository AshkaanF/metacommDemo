## ---
## initialize
## ---
# source('./R/tl_func.R')
set.seed(777)

## load libs
if(!require(dplyr)) install.packages('dplyr'); library(dplyr)
if(!require(magrittr)) install.packages('magrittr'); library(magrittr)
if(!require(ggplot2)) install.packages('ggplot2'); library(ggplot2)
if(!require(reshape2)) install.packages('reshape2'); library(reshape2)
if(!require(Matrix)) install.packages('Matrix'); library(Matrix)
if(!require(igraph)) install.packages('igraph'); library(igraph)
if(!require(expm)) install.packages('expm'); library(expm)
if(!require(NetIndices)) install.packages('NetIndices'); library(NetIndices)
if(!require(RColorBrewer)) install.packages('RColorBrewer'); library(RColorBrewer)
if(!require(e1071)) install.packages('e1071'); library(e1071)
if(!require(doParallel)) install.packages('doParallel'); library(doParallel)
if(!require(parallel)) install.packages('parallel'); library(parallel)
if(!require(foreach)) install.packages('foreach'); library(foreach)

##---
## (1) import and load data in the 'test' directory
##---
out.dir <- '../webs/' ## use this for all jacobian
web.list <- list.dirs(path = out.dir, recursive = F) %>% paste(., '/', sep = '') %>% gsub('//', '/', .)

## empty list to store results
master <- list()

## diffusion and z sweeps
diffs <- c(1e-3, 1e-2, 1e-1)
zs <- c(-0.75, 0, 0.75)
spatials <- expand.grid('dif' = diffs, 'z' = zs)

## number of random permutations to compare with
number.of.dispersal.permutations <- 10

## geometric graph
geo.graph <- function(n = 64, r = 0.32) {
  xy <- cbind(runif(n), runif(n))
  distMat <- as.matrix(dist(xy, method = 'euclidean', upper = T, diag = T))
  adjMat <- matrix(F, nr = n, nc = n)
  adjMat[ distMat < r ] <- T
  diag(adjMat) <- F
  return(list(xy, adjMat))
}

## ---
## WORKHORSE LOOP START
## ---

## register cluster
cl <- makeCluster(4, setup_strategy = 'sequential'); registerDoParallel(cl)

## init
for(w in 1:length(web.list)){

  ## temp directions
  temp.dir <- web.list[w]
  rep.list <- list.dirs(path = temp.dir, recursive = F) %>% paste(., '/', sep = '') %>% gsub('//', '/', .)
  
  ## storage for sub matrices
  sub.master <- list()
  sub.sub.master <- list()
  
  ## look at this replicate
  for(q in 1:length(rep.list)){
  
    ## temp directions
    temp.temp.dir <- rep.list[q]
    
    ##
    ## file lists
    ##
    jacobian.files <- list.files(path = temp.temp.dir, pattern = 'JacobianD', recursive = T) %>% paste(temp.temp.dir, ., sep = '')
    parm.files <- list.files(path = temp.temp.dir, pattern = 'ParametersD', recursive = T) %>% paste(temp.temp.dir, ., sep = '')
    adj.files <- list.files(path = temp.temp.dir, pattern = 'AdjacencyD', recursive = T) %>% paste(temp.temp.dir, ., sep = '')
    
    ## ---
    ## (2) PREP FOR CALCULATIONS
    ## ---
    ind.webs <- list()
    lev <- c()
    
    ## (a) load i's adjacency matrix
    data.adj <- as.matrix(read.table(adj.files[1])) # only need [1] because they are identical
    
    ## who is my top predator?
    top.predator <- which.max(TrophInd(Tij = t(data.adj))$TL)
    
    ## who are producers?
    producers <- which(colSums(data.adj) == 0)
    
    ## blank list for branching parms
    branches <- list()
    
    ## inner workhorse
    for(i in 1:length(jacobian.files)){
    
      ## (b) read jacobian
      data.jac <- as.matrix(read.table(jacobian.files[i]))
      
      ## prep parms for species
      temp.parm <- read.table(parm.files[i], fill = T)
      
      ## (c) split species level parameters
      temp.sp <- temp.parm[1:(apply(temp.parm, 1, FUN = function(p) any(is.na(p))) %>% which() %>% min() - 1), ]
      colnames(temp.sp) <- c('sp', 'alpha', 'gamma', 'mu', 'phi', 'psi', 'rho', 'sigm')
      
      ## (d) split flow level parameters
      temp.fl <- temp.parm[-c(1:(apply(temp.parm, 1, FUN = function(p) any(is.na(p))) %>% which() %>% min() - 1)), 1:5]
      colnames(temp.fl) <- c('from', 'to', 'chi_to', 'beta_fr', 'lambda')
      
      ## store
      branches[[i]] <- matrix(nrow = 1, temp.fl$chi_to, byrow = T)
      ind.webs[[i]] <- data.jac
      lev[i] <- data.jac %>% eigen %$% values %>% Re() %>% max()
    
    }
    
    ## store alphas
    alp <- temp.sp$alpha
    
    ##
    ## get body masses 
    ##
    masses <- alp ^ -4
    
    ##
    ## paste jacobians on diagonal
    ##
    JJ <- bdiag(ind.webs) %>% as.matrix()
    
    ## (e) extract S (richness) and C (connectance) and B (mass scaling, = body mass ratio ^ -1) from each web
    Ss <- nrow(data.adj)
    Cs <- sum(data.adj) / (Ss ^ 2)
    
    ##---
    ## (4) make it spatial
    ##---
    ## make empty dispersal matrices
    data.disp <- as.matrix(0 * data.jac)
    
    ##
    ## dispersal parm sweep
    ##
    sub.sub.master <- foreach(x = 1 : nrow(spatials), .inorder = TRUE, .packages = c('tidyverse', 'igraph', 'NetIndices'), .verbose = F) %dopar% {
      
      ## make sure you have a fully connected spatial graph
      for(h in 1:1e4){
        gg <- geo.graph(10, 0.32)
        g <- graph_from_adjacency_matrix(gg[[2]], mode = 'undirected')
        layout.g <- gg[[1]]
        if(components(g)$no == 1) break}
      
      ## dispersal formula parameters
      diffusion.coef <- spatials$dif[x]
      z <- spatials$z[x]
      
      ## populate. this is the C matrix based on body mass
      pre.disp <- diffusion.coef * masses ^ z 
      
      ## no normalization
      diag(data.disp) <- pre.disp
      
      ## calculate (i) adjacency, (ii) node degrees and (iii) populate matrix diagonal
      L.adj <- g %>% get.adjacency() %>% as.matrix
      L.degree <- L.adj * 0
      diag(L.degree) <- g %>% degree()
      
      ## make Laplacian
      L <- (L.degree - L.adj)
      
      ## ---
      ## (5) merging local jacobian, dispersal matrix, and spatial laplacian starts here
      ## ---
      ## grab dispersal matrix
      C <- data.disp
      
      ## mash it all together
      J <- JJ - (L %x% C)
      
      ## eigenvalues
      et.space <- eigen(J)$values
      
      ##
      ## sort and store
      ##
      l1 <- lev               ## largest leading evs among patch set
      l2 <- et.space          ## metaweb evs
      max.l2 <- max(Re(l2))
      
      ## ---
      
      ## vector to fill in with results of permutations
      vec.of.permutations <- c()
      
      ## -
      ## PERMUTATIONS, shuffle the diagonal of the dispersal matrix
      ## -
      for(pp in 1:number.of.dispersal.permutations){
      
        ## create a permuted version of data.disp (linearization of the movement part of our ODEs)
        C.permuted <- data.disp
        diag(C.permuted) <- sample(diag(C.permuted))
        
        ## mash it all together
        J.permuted <- JJ - (L %x% C.permuted)
        
        ## eigenvalues
        et.space.permuted <- eigen(J.permuted)$values
        
        ## store leading ev
        vec.of.permutations[pp] <- max(Re(et.space.permuted))
      }
      
      ## extract median
      median.eigenvalue <- median(vec.of.permutations)
      mean.eigenvalue <- mean(vec.of.permutations)
      se.eigenvalue <- sd(vec.of.permutations) / sqrt(length(vec.of.permutations))
      frac.shuffles.with.larger.eigenvalue <- mean(vec.of.permutations > max.l2)
      
      ## ---
      
      ## ---
      ## (6) data output
      ## ---
      
      ## store
      resy <- data.frame('S' = Ss, 'C' = Cs, 
                         'ev.1' = l1[1],    ## leading eigenvalue of each patch
                         'ev.2' = l1[2],    ## leading eigenvalue of each patch
                         'ev.3' = l1[3],    ## leading eigenvalue of each patch
                         'ev.4' = l1[4],    ## leading eigenvalue of each patch
                         'ev.5' = l1[5],    ## leading eigenvalue of each patch
                         'ev.6' = l1[6],    ## leading eigenvalue of each patch
                         'ev.7' = l1[7],    ## leading eigenvalue of each patch
                         'ev.8' = l1[8],    ## leading eigenvalue of each patch
                         'ev.9' = l1[9],    ## leading eigenvalue of each patch
                         'ev.10' = l1[10],  ## leading eigenvalue of each patch
                         'space.ev' = max(Re(l2)), ## leading eigenvalue of metacommunity
                         'median.permuted.eigenvalue' = median.eigenvalue, ## median eigenvalues of N permutations of dispersal rates
                         'mean.permuted.eigenvalue' = mean.eigenvalue,     ## mean eigenvalues of N permutations of dispersal rates
                         'se.permuted.eigenvalue' = se.eigenvalue,         ## std. err. of eigenvalues of N permutations of dispersal rates
                         'frac.permutations.with.larger.eigenvalue' = frac.shuffles.with.larger.eigenvalue,
                         'z' = z, 
                         'diff.co' = diffusion.coef,
                         'sim' = paste('S', Ss, 'C', Cs, 'sim', q, sep = '.')
                        )
      
      ## store
      resy
    
    }
    
    ## bind
    res.first <- do.call('rbind', sub.sub.master)
      
    ## hold this
    sub.master[[ q ]] <- res.first
    
  }
  
  ## bind
  res <- do.call('rbind', sub.master)
  
  ## store
  master[[w]] <- res
  
  ## drop a line
  cat('Directory', web.list[w], 'complete (', round(100 * w / length(web.list), 2), '%)\n')
  
}

## stop clust
stopCluster(cl)


## ---
## bind and prep output
## ---

## glue the list together and calculate some things
res <- do.call('rbind', master)
res$max.local <- apply(res[, grepl('ev.', names(res), fixed = T)], 1, max)
res$prop.stab <- apply(res[, grepl('ev.', names(res), fixed = T)], 1, function(v) sum(v < 0) / length(v))

##
## summarize
##
res.sum <- res %>%
  group_by(prop.stab, diff.co, z) %>%
  summarise(prop.reg = mean(space.ev < 0),
            se.reg = sd(space.ev < 0) / sqrt(length(space.ev))) %>%
  as.data.frame()


## plot the fraction of webs that are more stable with dispersal allometry
ggplot(data = res.sum, aes(x = prop.stab, y = prop.reg, fill = as.factor(z), colour = as.factor(z), group = z)) +
  theme_classic() +
  xlab('Prop. stable webs') +
  ylab('Prop. stable metacommunities') +
  geom_errorbar(aes(ymax = prop.reg + 1*se.reg, 
                    ymin = prop.reg - 1*se.reg), 
                width = 0, size = 0.5, position = position_dodge(width = 0.045)) +
  geom_point(shape = 21, size = 2, alpha = 0.9, stroke = 0.75, colour = '#000000', 
             position = position_dodge(width = 0.045)) +
  scale_fill_manual(values = brewer.pal(length(zs), 'RdBu'), name = expression('Dispersal exponent, z'),
                       guide = guide_legend(label = TRUE,
                                            title.position = 'left',
                                            title.theme = element_text(size = 7, angle = 90),
                                            label.position = 'right',
                                            direction = 'vertical')) +
  scale_colour_manual(values = brewer.pal(length(zs), 'RdBu'), guide = F) +
  facet_wrap(~ diff.co, nrow = 1) +
  theme(legend.position = 'right',
        legend.title = element_text(size = 8), legend.text = element_text(size = 5.5),
        panel.border = element_blank(), panel.background = element_blank(),
        strip.background = element_rect(fill = '#d0d1e6', size = 0), 
        strip.text = element_text(size = 10, margin = margin(0.02, 0, 0.02, 0, "cm")),
        axis.title = element_text(size = 10))

