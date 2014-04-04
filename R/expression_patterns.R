#` return a named matrix of patterns from a list
pattern_matrix <- function(pattern_list, ncond) {
  patterns <- t(matrix(unlist(pattern_list), nrow=ncond))
  rownames(patterns) <- names(pattern_list)
  colnames(patterns) <- as.character(1:ncond)
  return(patterns)
}

#' Write out DE information for each DE gene expression pattern
output_pattern_sets <- function(de_data, conditions, 
                                named_patterns, prob_cutoff) {
  final <- de_data[['final']]
  mean_cols <- de_data[['mean_cols']]
  de_prob_cols <- de_data[['de_prob_cols']]
  ee_prob_col <- de_data[['ee_prob_col']]
  prob_cols <- c(de_prob_cols, ee_prob_col)
  patterns <- de_data[['results']]
  # add patterns
  n <- length(unique(conditions))
  if (n==2) {
    final$pattern <- apply(final[,mean_cols],
                         1,
                         function(x) paste(pattern(x), collapse='_'))
  } else {
    final$pattern <- multiway_directional_patterns(de_data, prob_cutoff)
  }
  num_patterns <- lapply(unique(final$pattern),function(x) {
    if (x == "no significant pattern") {
      return(paste(rep(0, length(unique(conditions))), collapse='_'))
    } else {
      x
    }
  })

  # genes between EE and DE cutoffs should have no pattern
  final$pattern[which(apply(final[,prob_cols], 1, 
    function(x){ 
      all(x < prob_cutoff) 
    }))] <- "no significant pattern"
  # same for genes with counts too low to infer DE
    final$pattern[which(apply(final[,prob_cols], 1, 
    function(x){ 
      any(is.na(x)) 
    }))] <- "no significant pattern"
  # genes above EE cutoff are always labelled as equally expressed
  final$pattern[which(final[,ee_prob_col] >= prob_cutoff)] <- "equal expression"
  # genes with flat pattern too
  flat_pattern <- paste(rep('1', n), collapse="_")
  flat_pattern_idx <- which(sapply(final$pattern, function(x) { x == flat_pattern}))
  final$pattern[flat_pattern_idx] <- "equal expression"
  # replace named patterns, saving the originals
  final$raw_pattern <- final$pattern
  final$pattern <- replace_patterns(final$pattern, named_patterns)
  # select probable DE/EE genes above cutoff
  sig <- final[which(apply(final[,prob_cols], 1, function(x) {any(x >= prob_cutoff)})),]
  if (!nrow(sig)) {
    stop("There are no rows with significantly differential or equal expression")
  } else {
    print(paste("There were", nrow(sig), "rows (out of", nrow(final), 
      "tested) with significantly differential or equal expression (PP >=", prob_cutoff, ")"))
    write.table(x=table(sig$pattern),
                file="significant_genes_per_pattern.csv",
                sep=",",
                row.names=F,
                col.names=T)
  }
  return(list(final=final, num_patterns=num_patterns))
}

multiway_directional_patterns <- function(de_data, prob_cutoff) {
  final <- de_data[['final']]
  mean_cols <- de_data[['mean_cols']]
  de_prob_cols <- de_data[['de_prob_cols']]
  ee_prob_col <- de_data[['ee_prob_col']]
  prob_cols <- c(de_prob_cols, ee_prob_col)
  patterns <- de_data[['results']]$AllParti
  final$pattern <- "no significant pattern"

  for (pattern in prob_cols) {
    pat_string <- paste(patterns[pattern,], collapse="_")
    print(paste("converting pattern", pat_string, "to directional"))
    pat_rows <- which(final[,pattern] >= prob_cutoff)
    if (length(pat_rows) == 0) next
    means <- final[pat_rows, mean_cols]
    basic_pattern <- patterns[pattern,]
    names(means) <- names(basic_pattern)
    all <- all_directional_patterns(basic_pattern)
    all.sorted <- lapply(all, sort)
    final$pattern[pat_rows] <- apply(means, 1, function(x) {
      patrows <- lapply(all.sorted, function(y) {
        all(names(y) == names(sort(x)))
      })
      pat <- all[which(unlist(patrows))]
      pat <- paste(unlist(pat), collapse="_")
      if (pat == "") {
        print("wtf")
      }
      return(pat)
    })
  }
  return(final$pattern)
}

#' Generate all valid directional patterns corresponding to
#' an undirected (difference) pattern.
all_directional_patterns <- function(x) {
  x <- unlist(x)
  uniq <- unique(x)
  if (length(uniq) == 1) {
    return(x)
  } else if (length(uniq) == length(x)) {
    return(permute_pattern(x))
  } else {
    # the pattern contains repeats
    orig.idx <- NULL
    decomp_reps <- decompose_repeats(x)
    reps <- decomp_reps[['reps']]
    noreps <- decomp_reps[['noreps']]
    if (length(noreps) != length(unique(noreps))) {
      # there are non-consecutive repeats, columns need reordering
      # before repeat name permutation, and then the column ordering
      # needs to be restored
      orig <- names(x)
      sorted <- sort(x)
      # capture the original ordering
      orig.idx <- sapply(names(sorted), function(y) {
        which(orig == y)
      })
      x <- sort(x)
      # and re-decompose the repeats
      decomp_reps <- decompose_repeats(x)
      reps <- decomp_reps[['reps']]
      noreps <- decomp_reps[['noreps']]
    }
    perms <- permute_pattern(noreps)
    recomp_perms <- list()
    for (perm in perms) {
      # permute the names for the repeats
      r <- recompose_repeats(perm, reps)
      names(r) <- names(x)
      name_perms <- permute_repeat_names(list(r), reps)
      # output one list per permutation
      for (p in name_perms) {
        r2 <- r
        names(r2) <- p
        recomp_perms <- c(recomp_perms, list(r2))
      }
    }
    print(recomp_perms)
    # restore column ordering if it was changed
    if (!is.null(orig.idx)) {
      recomp_perms <- 
        lapply(recomp_perms, function(perm) {
          newperm <- list()
          i <- 1
          for (m in names(perm)) {
            idx <- orig.idx[[i]]
            newperm[idx] <- perm[[m]]
            names(newperm)[idx] <- m
            i <- i + 1
          }
          return(unlist(newperm))
        })
    }
    return(recomp_perms)
  }
}

#' Return the mirror of the pattern of numbers.
#' Note: this is *not* the same as reversing the pattern.
#' Examples:
#' mirror_sequence(c(1,2,2,3))
#' c(3,2,2,1)
#' mirror_sequence(c(1,2,3,2))
#' c(3,2,1,2)
mirror_pattern <- function(x) {
  return((max(x)+1)-x)
}

#' Get all permutations of the pattern
permute_pattern <- function(x) {
  get_package('combinat')
  p <- permn(x)
  p <- lapply(p, function(y) {
    names(y) <- names(x)
    return(y)
  })
  return(unique(p))
}

#' Identify all repeats in the container, returning a vector of the original
#' sequence with repeats removed, and a data frame
#' with one row per repeat, the first column being the start position
#' (1-indexed) and the second column being the length
decompose_repeats <- function(x) {
  reps <- data.frame()
  idx <- 0
  replen <- 1
  last <- NULL
  i <- 1
  noreps <- c()
  for (y in x) {
    seqend <- (i == length(x))
    repend <- FALSE
    isrep <- FALSE
    if (i == 1) {
      # first entry
      last <- y
      idx <- i
      i = i + 1
      noreps <- append(noreps, y)
      next
    } else if (last == y) {
      # we're in a repeat
      isrep <- TRUE
      replen = replen + 1
      if (!seqend) {
        i = i + 1
        next
      }
    }
    if (replen > 1) {
      # we've just finished a repeat
      reps <- rbind(reps, c(idx, replen))
      repend <- TRUE
    }
    if (!(seqend && repend) || !isrep) {
      noreps <- append(noreps, y)
    }
    # just a regular entry
    replen <- 1
    idx <- i
    last <- y
    i = i + 1
  }
  if (ncol(reps) == 2)
    names(reps) <- c('start', 'length')
  return(list(noreps=noreps, reps=reps))
}

#' Insert repeats into a pattern vector.
recompose_repeats <- function(x, reps) {
  y <- x
  if (nrow(reps)) {
    for (i in 1:nrow(reps)) {
      start <- reps$start[i]
      len <- reps$length[i]
      seq <- rep(y[start], len-1)
      y <- append(y, seq, after=start)
    }
  }
  return(y)
}

#' Take a pattern with repeats and for each repeat,
#' output a list of new patterns with each possible permutation
#' of the names for the repeat elements.
permute_repeat_names <- function(x, reps) {
  x <- unlist(x, recursive=FALSE)
  all <- list()
  laststop <- 0
  if (nrow(reps)) {
    # collect all name components
    for (i in 1:nrow(reps)) {
      start <- reps$start[i]
      len <- reps$length[i]
      stop <- start + len - 1
      if (start > (laststop + 1)) {
        # insert preceding names that don't need permuting
        noperm_names <- names(x)[(laststop+1):(start-1)]
        noperm_names <- sapply(noperm_names, as.character)
        all <- append(all, noperm_names)
      } 
      # insert permuted names
      permute_names <- names(x)[start:stop]
      all[[length(all)+1]] <- permn(permute_names)
      laststop <- stop
    }
    if (laststop < length(x)) {
      all <- append(all, names(x)[(laststop+1):length(x)])
    }
  }
  # permute all combinations of the sets of repeat names
  all.grid <- expand.grid(all, stringsAsFactors=FALSE)
  all.split <- split(all.grid, rownames(all.grid))
  perms <- lapply(all.split, function(y) {
    pattern <- x
    names(pattern) <- unlist(y)
  })
  return(perms)
}

#' Reduce a sequence of numbers to its pattern of changes
pattern <- function(x, p=c(1), i=1, j=2) {
  if (x[i] == x[j]) {
    # no change
    n <- p[length(p)]
  } else if (x[i] > x[j]) {
    # decrease
    n <- p[length(p)] - 1
  } else {
    # increase
    n <- p[length(p)] + 1
  }
  # recurse along the pattern
  if (length(x) > j) {
    n <- pattern(x, n, i+1, j+1)
  }
  # add this to the pattern so far
  pattern <- c(p, n)
  # before the top-level return, adjust to so min is 1
  if (length(pattern) == length(x) && min(pattern) < 1) {
    pattern = pattern + abs(min(pattern)) + 1
  }
  return(pattern)
}

#' Replace patterns with the specified strings
replace_patterns <- function(patterns, named_patterns) {
  sapply(patterns,
    function(x) {
      if (x %in% names(named_patterns)) {
        return(named_patterns[[x]])
      } else {
        return(x)
      }
    })
}

#' Extract the pattern corresponding to equal expression
#' across all conditions. Returns a character vector of the pattern
#' name.
get_EE_pattern <- function(patterns) {
  flatrows <- apply(patterns, 1, function(x) { sum(x) == length(x) })
  return(rownames(patterns)[which(flatrows)]) 
}

#' sort a names list according to the specified name ordering
#' the 'order' container must contain all names in the list
sort_list <- function(x, order) {
  sorted <- list()
  for (n in order) {
    sorted <- c(sorted, x[[n]])
  }
  return(sorted)
}