library(R6)

#' @title SIMBLE
#'
#' @description
#' Simulation of Breakpoint Ligation Evolution.
#'
#'
#' @import R6
#' @import data.table
#' @import GenomicRanges
#' @import IRanges
#' @import gUtils
#' @import gTrack
#' @import JaBbA
#' @import gGnome
#' @importFrom gUtils %&%
NULL


#' simble
#'
#' Simulation object
#' Object has an "epoch" state
#' Each epoch is a state

simble = R6Class("simble",
                 public = list(
                     ## epoch = NULL,
                     ## tile = NULL,
                     ## junctions = NULL,
                     initialize = function(junctions = NULL, sort = TRUE, seed = NULL, sex = NULL) {
                         if (! is.null(junctions)) {
                             if (inherits(junctions, "character")) {
                                 if (grepl("\\.vcf$", junctions)) {
                                     junctions = gGnome::ra_breaks(junctions)
                                 } else {
                                     stop("Unrecognized file type. extension must be .vcf")
                                 }
                             }
                             if (inherits(junctions, "GRangesList")) {
                                 orig_junctions = junctions
                                 junctions = endoapply(sortSeqlevels(junctions), function(gr) sort(gr, ignore.strand = TRUE))
                                 remove_idx = unique(gr2dt(grl.unlist(junctions))[! seqnames %in% c(1:24, "X", "Y")]$grl.ix)
                                 if (length(remove_idx) > 0) {
                                     junctions = junctions[-remove_idx]
                                 }
                                 tmp = private$karyograph(junctions)$tile
                                 if (sort) {
                                     seqlevels(tmp) = seqlevels(tmp)[orderSeqlevels(seqlevels(tmp), TRUE)]
                                     tmp = sort(tmp)
                                     pos_tile = tmp[as.character(strand(tmp)) == "+"]
                                     new_tmp = c(pos_tile, pos_tile)
                                     names(new_tmp) = c(1:length(pos_tile), -(1:length(pos_tile)))
                                     strand(new_tmp[(length(new_tmp)/2+1):length(new_tmp)]) = "-"
                                     tmp = new_tmp
                                 }

                                 private$tile = tmp
                                 private$pos_tile = pos_tile
                                 private$junctions = junctions
                                 discarded_junc = orig_junctions[elementNROWS(setdiff(orig_junctions, private$junctions)) > 0]
                                 private$discarded_junc = discarded_junc

                                 message("Junctions supplied to simble: \n", length(orig_junctions))
                                 message("Junctions that will be used for simulation: \n", length(private$junctions))
                                 message("Junctions discarded (junctions found on mitochondrial or decoy chromosomes): \n", length(private$discarded_junc))
                                 message("Genome broken into GRanges tile of length: \n", length(private$tile))

                                         
                                 ########## initializing contigs
                                 grep_m = c("Male", "M", "male", "m")
                                 grep_f = c("Female", "F", "female", "f")
                                 
                                 if (is.null(sex)) {
                                     set.seed(8)
                                     sex = sample(c("Male", "Female"), 1)
                                 }
                                 
                                 junc = private$junctions
                                 u_junc = grl.unlist(junc)
                                 u_junc = GenomeInfoDb::sortSeqlevels(u_junc)
                                 ## seqlevels(u_junc) = seqlevels(u_junc)[orderSeqlevels(seqlevels(u_junc), TRUE)]
                                 
                                 sort_u_junc = sort(u_junc, ignore.strand = T)
                                 sort_u_junc$idx = 1:length(u_junc)


                                 ## segged = keepSeqlevels(private$pos_tile, c(1:22, "X", "Y"))
                                 segged = private$gr_keep_seq(pos_tile, as.character(c(1:22, "X", "Y")))
                                 segged = sort(GenomeInfoDb::sortSeqlevels(segged))
                                 ## seqlevels(segged) = seqlevels(segged)[orderSeqlevels(seqlevels(segged), TRUE)]
                                 ## segged = sort(segged)
                                 segged$id = seg_int = 1:length(segged)

                                 ref_contigs = split(segged, seqnames(segged), drop = F)
                                 ref_int_contigs = contigs_lst = as.list(split(seg_int, seqnames(segged)))
                                 tmp_tel_1 = unlist(lapply(ref_int_contigs, function(x) head(x,1)), use.names = F)
                                 tmp_tel_2 = unlist(lapply(ref_int_contigs, function(x) tail(x,1)), use.names = F)
                                 left_tel = c(tmp_tel_1, -tmp_tel_2)
                                 right_tel = c(-tmp_tel_1, tmp_tel_2)
                                 tel_seg =  unlist(lapply(ref_int_contigs, function(x) c(head(x,1), tail(x, 1))), use.names = F)
                                 rev_contigs  = lapply(contigs_lst, function(x) -rev(x))

                                 ## seg_idx = nearest(sort_u_junc, segged, ignore.strand = T)
                                 ## names(seg_idx) = as.character(strand(sort_u_junc))
                                 map = private$map_junc2seg(sort_u_junc, segged)
                                 values(sort_u_junc)[map$queryHits,"seg_idx"] = as.integer(names(segged)[map$subjectHits])
                                 


                                 ## sort_u_junc$seg_idx = seg_idx

                                 dt_junc = gr2dt(sort_u_junc[,c("grl.ix", "grl.iix", "seg_idx")])
                                 connections = dt_junc[,list(seqnames, seg_idx, strand), by = grl.ix]
                                 connections[, strand_val := ifelse(strand == "-", 1, 2)]

                                 lst = split(connections, connections$grl.ix)
                                 bint_lst = lapply(lst, function(x) {
                                     tmp = x[,seg_idx] * (x[,strand_val]*2-3)
                                     ## tmp = ifelse(tmp > 0, tmp + 1, tmp) ## no longer relevant
                                     names(tmp) = as.character(x[,strand])
                                     return(tmp)
                                 })

                                 private$j_segments = bint_lst

                                 junc_pool = 1:length(junc)
                                 
                                 iteration_list = NULL

                                 hold_out = data.table(contigs_lst = list(),
                                                       contig_nm = character(),
                                                       contig_idx = integer(),
                                                       iter_idx = integer(),
                                                       loose_left = logical(),
                                                       loose_right = logical())

                                 ## Utility functions ##
                                 
                                 rename_cont_lst = function(subcont, nms) {
                                     mapply(function(subcont, nms) { names(subcont)=nms; return(subcont) }, subcont, nms, SIMPLIFY = FALSE)
                                 }
                                 
                                 
                                 autosomal_hi = sapply(contigs_lst[! names(contigs_lst) %in% c("X", "Y")],
                                                       function(x) {
                                                           names(x) = paste0(rep("hi", length(x)), "_", x)
                                                           return(x)
                                                       })
                                 autosomal_lo = sapply(contigs_lst[! names(contigs_lst) %in% c("X", "Y")],
                                                       function(x) {
                                                           names(x) = paste0(rep("lo", length(x)), '_', x)
                                                           return(x)
                                                       })
                                 
                                 diploid_autosomal = c(autosomal_hi, autosomal_lo)
                                 current_contigs = data.table(contigs_lst = diploid_autosomal,
                                                              contig_nm = names(diploid_autosomal),
                                                              contig_idx = seq_along(diploid_autosomal),
                                                              iter_idx = 0,
                                                              loose_left = FALSE,
                                                              loose_right = FALSE) 

                                 tmp_seg_1 = contigs_lst$X
                                 names(tmp_seg_1) = paste0(rep("hi", length(tmp_seg_1)), "_", tmp_seg_1)
                                 tmp_seg_1 = list(X = tmp_seg_1)

                                 
                                 if (!"Y" %in% seqlevelsInUse(orig_junctions)) {
                                     ## if (rbinom(1, 1, prob = 0.5)) {
                                     if (sex %in% grep_f) {
                                         ## if (1 == 0) {
                                         tmp_seg_2 = unlist(tmp_seg_1)
                                         names(tmp_seg_2) = paste0(rep("lo", length(tmp_seg_2)), "_", tmp_seg_2)
                                         tmp_seg_2 = list(X = tmp_seg_2)
                                     } else if (sex %in% grep_m) {
                                         tmp_seg_2 = contigs_lst$Y
                                         names(tmp_seg_2) = paste0(rep("lo", length(tmp_seg_2)), "_", tmp_seg_2)
                                         tmp_seg_2 = list(Y = tmp_seg_2)
                                     }
                                 } else {
                                     if (sex %in% grep_f) {
                                         message("Junctions on chr Y detected \nOverriding female designation")
                                     }
                                     tmp_seg_2 = contigs_lst$Y
                                     names(tmp_seg_2) = paste0(rep("lo", length(tmp_seg_2)), '_', tmp_seg_2)
                                     tmp_seg_2 = list(Y = tmp_seg_2)
                                 }
                                 
                                 sex_chrom = c(tmp_seg_1, tmp_seg_2)
                                 sex_chrom_dt = data.table(contigs_lst = sex_chrom,
                                                           contig_nm = names(sex_chrom),
                                                           contig_idx = c(45, 46),
                                                           iter_idx = 0,
                                                           loose_left = FALSE,
                                                           loose_right = FALSE)

                                 current_contigs = rbind(current_contigs, sex_chrom_dt)
                                 private$epoch_0 = current_contigs
                                 
                                 message("\n\nEpoch 0 instantiated\n\n")

                                 ######### scrambling junctions

                                 message("Establishing scrambled index of junctions for simulation...\n\n")
                                 if (is.null(seed)) {
                                     message("No seed set. \nDefaulting to set.seed(1987)")
                                     set.seed(1987)
                                 } else {
                                     if (! inherits(seed, "numeric") & ! inherits(seed, "integer")) {
                                         stop("seed arg must be set to class numeric or class integer value")
                                     }
                                     message("seed set to: ", seed) 
                                     set.seed(seed)
                                 }


                                 scrambled_index = sample(1:length(private$junctions), length(private$junctions))
                                 private$scrambled_index = scrambled_index

                                 message("\n\nScrambled index of junctions assigned.")

                                 private$current_index = scrambled_index


                                 ########## number of total steps
                                 private$tot_steps = length(junctions)
                                 
                             }
                         } else {
                             stop("junctions must be supplied as .vcf file or a GRangesList object")
                         }
                     },
                     initialize_contigs = function(sex = NULL) {
                         grep_m = c("Male", "M", "male", "m")
                         grep_f = c("Female", "F", "female", "f")

                         if (is.null(sex)) {
                             set.seed(8)
                             sex = sample(c("Male", "Female"), 1)
                         }
                         
                         junc = private$junctions
                         u_junc = grl.unlist(junc)
                         seqlevels(u_junc) = seqlevels(u_junc)[orderSeqlevels(seqlevels(u_junc), TRUE)]
                         
                         sort_u_junc = sort(u_junc, ignore.strand = T)
                         sort_u_junc$idx = 1:length(u_junc)


                         segged = keepSeqlevels(private$pos_tile, c(1:22, "X", "Y"))
                         seqlevels(segged) = seqlevels(segged)[orderSeqlevels(seqlevels(segged), TRUE)]
                         segged = sort(segged)
                         segged$id = seg_int = 1:length(segged)

                         ref_contigs = split(segged, seqnames(segged), drop = F)
                         ref_int_contigs = contigs_lst = as.list(split(seg_int, seqnames(segged)))
                         tmp_tel_1 = unlist(lapply(ref_int_contigs, function(x) head(x,1)), use.names = F)
                         tmp_tel_2 = unlist(lapply(ref_int_contigs, function(x) tail(x,1)), use.names = F)
                         private$left_tel = left_tel = c(tmp_tel_1, -tmp_tel_2)
                         private$right_tel = c(-tmp_tel_1, tmp_tel_2)

                         
                         tel_seg =  unlist(lapply(ref_int_contigs, function(x) c(head(x,1), tail(x, 1))), use.names = F)
                         rev_contigs  = lapply(contigs_lst, function(x) -rev(x))

                         seg_idx = nearest(sort_u_junc, segged, ignore.strand = T)
                         names(seg_idx) = as.character(strand(sort_u_junc))

                         sort_u_junc$seg_idx = seg_idx

                         dt_junc = gr2dt(sort_u_junc[,c("grl.ix", "grl.iix", "seg_idx")])
                         connections = dt_junc[,list(seqnames, seg_idx, strand), by = grl.ix]
                         connections[, strand_val := ifelse(strand == "-", 1, 2)]

                         lst = split(connections, connections$grl.ix)
                         bint_lst = lapply(lst, function(x) {
                             tmp = x[,seg_idx] * (x[,strand_val]*2-3)
                             tmp = ifelse(tmp > 0, tmp + 1, tmp)
                             names(tmp) = as.character(x[,strand])
                             return(tmp)
                         })

                         private$j_segments = bint_lst

                         junc_pool = 1:length(junc)
                         
                         iteration_list = NULL

                         ## hold_out = data.table(contigs_lst = list(),
                         ##                       contig_nm = character(),
                         ##                       contig_idx = integer(),
                         ##                       iter_idx = integer(),
                         ##                       loose_left = logical(),
                         ##                       loose_right = logical())

                         ## Utility functions ##
                         
                         rename_cont_lst = function(subcont, nms) {
                             mapply(function(subcont, nms) { names(subcont)=nms; return(subcont) }, subcont, nms, SIMPLIFY = FALSE)
                         }
                         
                         
                         autosomal_hi = sapply(contigs_lst[! names(contigs_lst) %in% c("X", "Y")],
                                               function(x) {
                                                   names(x) = paste0(rep("hi", length(x)), "_", x)
                                                   return(x)
                                               })
                         autosomal_lo = sapply(contigs_lst[! names(contigs_lst) %in% c("X", "Y")],
                                               function(x) {
                                                   names(x) = paste0(rep("lo", length(x)), '_', x)
                                                   return(x)
                                               })
                         
                         diploid_autosomal = c(autosomal_hi, autosomal_lo)
                         current_contigs = data.table(contigs_lst = diploid_autosomal,
                                                      contig_nm = names(diploid_autosomal),
                                                      contig_idx = seq_along(diploid_autosomal),
                                                      iter_idx = 0,
                                                      loose_left = FALSE,
                                                      loose_right = FALSE) 

                         tmp_seg_1 = contigs_lst$X
                         names(tmp_seg_1) = paste0(rep("hi", length(tmp_seg_1)), "_", tmp_seg_1)
                         tmp_seg_1 = list(X = tmp_seg_1)

                         
                         if (!"Y" %in% seqlevelsInUse(orig_junctions)) {
                             ## if (rbinom(1, 1, prob = 0.5)) {
                             if (sex %in% grep_f) {
                                 ## if (1 == 0) {
                                 tmp_seg_2 = unlist(tmp_seg_1)
                                 names(tmp_seg_2) = paste0(rep("lo", length(tmp_seg_2)), "_", tmp_seg_2)
                                 tmp_seg_2 = list(X = tmp_seg_2)
                             } else if (sex %in% grep_m) {
                                 tmp_seg_2 = contigs_lst$Y
                                 names(tmp_seg_2) = paste0(rep("lo", length(tmp_seg_2)), "_", tmp_seg_2)
                                 tmp_seg_2 = list(Y = tmp_seg_2)
                             }
                         } else {
                             if (sex %in% grep_f) {
                                 message("Junctions on chr Y detected \nOverriding female designation")
                             }
                             tmp_seg_2 = contigs_lst$Y
                             names(tmp_seg_2) = paste0(rep("lo", length(tmp_seg_2)), '_', tmp_seg_2)
                             tmp_seg_2 = list(Y = tmp_seg_2)
                         }
                         
                         sex_chrom = c(tmp_seg_1, tmp_seg_2)
                         sex_chrom_dt = data.table(contigs_lst = sex_chrom,
                                                   contig_nm = names(sex_chrom),
                                                   contig_idx = c(45, 46),
                                                   iter_idx = 0,
                                                   loose_left = FALSE,
                                                   loose_right = FALSE)

                         current_contigs = rbind(current_contigs, sex_chrom_dt)
                         return(current_contigs)
                     },
                     iterate_sim = function(num_steps = length(private$current_index), padding = 3e4) {
                         these_indices = 1:num_steps

                         ## if the user supplies more steps than there are junctions
                         ## left to sample... then return message
                         if (length(these_indices) > length(private$current_index)) {
                             how_many_iterations_left = length(private$current_index)
                             message("num_steps supplied: ", num_steps)
                             message("number of iterations left: ", how_many_iterations_left)
                             these_indices = 1:how_many_iterations_left
                             message("num_steps defaulting to: ", how_many_iterations_left, '\n\n')
                         }
                         
                         these_scrambled_indices = private$current_index[these_indices]
                         curr_step = private$curr_step
                         private$past_index = c(private$past_index, these_scrambled_indices)
                         
                         for (i in these_scrambled_indices) {
                             if (curr_step == 1) {
                                 prev_contigs = private$epoch_0
                             } else {
                                 prev_step = curr_step - 1
                                 prev_contigs = private$epoch[[prev_step]]
                             }

                             prev_contigs[, n := as.character(1:.N)]
                             prev_contigs[, c('loose_left', 'loose_right') := list(unlist(loose_left), unlist(loose_right))]
                             
                             junc_match = private$j_segments[[i]]
                             junction_grl = private$junctions[[i]]
                             junction_grl_val = values(private$junctions[i])

                             tmp_match = private$contigs2dt(prev_contigs$contigs_lst, private$tile)
                             setnames(tmp_match, "unlisted", "seg_id")

                             tmp_match[, loose_left := FALSE]
                             tmp_match[, loose_right := FALSE]
                             setkey(prev_contigs, n)

                             tmp_match[, contig_idx := prev_contigs[as.character(.ix), contig_idx]]
                             ## tmp_match[, pool := prev_contigs[as.character(.ix), pool]]
                             tmp_match[, loose_left := ifelse(.iix == .iix[1],
                                                              prev_contigs[as.character(.ix), loose_left],
                                                              FALSE), by = .ix]
                             tmp_match[, loose_right := ifelse(.iix == .iix[.N],
                                                               prev_contigs[as.character(.ix), loose_right],
                                                               FALSE), by = .ix]
                             tmp_match[, loose_left_reference_face := ifelse(loose_left,
                                                                      ifelse(strand == "+", "left", "right"),
                                                                      NA)]
                             tmp_match[, loose_right_reference_face := ifelse(loose_right,
                                                                       ifelse(strand == "+", "right", "left"),
                                                                       NA)]

                             tmp_match[, is_terminal := .iix %in% c(.iix[1], .iix[.N]), by = .ix]

                             tmp_match[, seqnames := factor(seqnames, c(as.character(1:24), "X", "Y", "M"))] ## necessary for later to order by reference distance

                             
                             tmp_match[, junc_match1 := ifelse(abs(junc_match[1]) == abs(seg_id), junc_match[1], NA)]
                             tmp_match[, junc_match2 := ifelse(abs(junc_match[2]) == abs(seg_id), junc_match[2], NA)]

                             tmp_match[,junc_contig_face1 :=  ifelse(!is.na(junc_match1),
                                                              ifelse(sign(junc_match1) == sign(seg_id),
                                                                     "right",
                                                                     "left"), NA)]

                             tmp_match[,junc_contig_face2 :=  ifelse(!is.na(junc_match2),
                                                              ifelse(sign(junc_match2) == sign(seg_id),
                                                                     "right",
                                                                     "left"), NA)]




                             tmp_match[, cumdist_from_left := cumsum(width), by = .ix]
                             tmp_match[, cumdist_from_right := rev(cumsum(rev(width))), by = .ix]

                             tmp_match[, num_segs := .N, by = .ix]
                             tmp_match[, num_from_left := 1:.N, by = .ix]
                             tmp_match[, num_from_right := rev(1:.N), by = .ix]

                             tmp_match[, num_rearrangements := length(which(diff(seg_id) != 1)), by = .ix]
                             tmp_match[, num_rearrangements_from_left := cumsum(c(0, ifelse(diff(seg_id) == 1, 0, 1))), by = .ix]
                             tmp_match[, num_rearrangements_from_right := rev(cumsum(c(0, ifelse(diff(rev(seg_id)) == -1, 0, 1)))), by = .ix]
                             

                             

                             tmp_match = tmp_match[order(.ix, seqnames, start, end)] ## order by reference distance with respect to every contig (hence .ix first in ordering)


                             ## get reference distances of adjacent segments in contigs 
                             tmp_match[, ref_dist_to_next := c(as.numeric(ifelse(seqnames[-1] == seqnames[-.N],
                                                                                 c(start[-1]) - c(end[-.N]),
                                                                                 -1)), NA), by = .ix]

                             tmp_match[, ref_dist_to_prev := c(NA,as.numeric(ifelse(seqnames[-.N] == seqnames[-1],
                                                                                    c(start[-1]) - c(end[-.N]),
                                                                                    -1))), by = .ix]

                             tmp_match = tmp_match[order(seqnames, start, end)] ## now order with respect to the whole genome


                             gr_match = gUtils::dt2gr(tmp_match)
                             

                             all_matches = private$matching_algo2(gr_match, junction_grl, junc_match, padding = padding)

                             match_lst1 = private$select_match(all_matches, 1, junc_match)
                             match_lst2 = private$select_match(all_matches, 2, junc_match)
                             

                             is_on_same_contig = match_lst1$final_match$.ix == match_lst2$final_match$.ix

                             ra1 = private$process_matching_algo(match_lst1, tmp_match, 1, prev_contigs)
                             ra2 = private$process_matching_algo(match_lst2, tmp_match, 2, prev_contigs)

                             combined_contigs = list(ra1$single_ra$contigs_lst[[1]], ra2$single_ra$contigs_lst[[1]])
                             combined_indices = list(ra1$adjusted_indices, ra2$adjusted_indices)
                             combined_matches = c(ra1$single_ra$cont_match, ra2$single_ra$cont_match)
                             
                             ## in_current = ra1$in_current | ra2$in_current
                             ## in_hold_out = ra1$in_hold_out | ra2$in_hold_out

                             ## need to adjust the indices of the matches
                             ## only matters if the matches are on the same contig
                             
                             ## adjust ra1 by ra2 indices
                             if (is_on_same_contig) {
                                 

                                 
                                 tmp_ra = list(copy(ra1$single_ra), copy(ra2$single_ra))
                                 tmp_contig = combined_contigs[[1]][ra1$adjusted_indices]

                                 if (length(ra1$add_to_first) > 0 |
                                     length(ra2$add_to_first) > 0) {
                                     select = which.max(c(length(ra1$add_to_first),
                                                          length(ra2$add_to_first)))
                                 } else {
                                     select = NULL
                                 }
                                 will_fuse_loose = c(tmp_ra[[1]]$will_fuse_loose, tmp_ra[[2]]$will_fuse_loose)
                                 new_loose_index = ifelse(will_fuse_loose, ifelse(combined_matches == 1, "first", "last"), NA)
                                 
                                 other = (1:2)[-1*select]
                                 add_to_first = list(ra1$add_to_first, ra2$add_to_first)[select][1][[1]]
                                 ## new_cont_match = combined_indices[[select]][combined_matches[other]]
                                 new_cont_match = combined_indices[select][1][[1]][combined_matches[other]]
                                 combined_matches[other] = new_cont_match
                                 segs_add_first = combined_contigs[select][1][[1]][add_to_first]

                                 ## if (!isNone(other)) {
                                 ##     tmp_ra[[other]][, cont_match := new_cont_match]
                                 ## }
                                 if (length(ra1$add_to_last) > 0 |
                                     length(ra2$add_to_last) > 0) {
                                     select = which.max(c(length(ra1$add_to_last),
                                                          length(ra2$add_to_last)))
                                 } else {
                                     select = NULL
                                 }
                                 
                                 add_to_last = list(ra1$add_to_last, ra2$add_to_last)[select][1][[1]]
                                 segs_add_last = combined_contigs[select][1][[1]][add_to_last]
                                 tmp_contig = c(segs_add_first, tmp_contig, segs_add_last)

                                 new_contigs_len = length(tmp_contig)

                                 ensure_correct_contig_matches = ifelse(!is.na(new_loose_index),
                                                                 ifelse(new_loose_index == "first", 1, length(tmp_contig)),
                                                                 combined_matches)
                                 tmp_ra[[1]][, cont_match := ensure_correct_contig_matches[1]]
                                 tmp_ra[[2]][, cont_match := ensure_correct_contig_matches[2]]
                                 tmp_ra[[1]][, contigs_len := new_contigs_len]
                                 tmp_ra[[2]][, contigs_len := new_contigs_len]
                                 tmp_ra[[1]][, contigs_lst := list(list(tmp_contig))]
                                 tmp_ra[[2]][, contigs_lst := list(list(tmp_contig))]

                                 ra1$single_ra = tmp_ra[[1]]
                                 ra2$single_ra = tmp_ra[[2]]
                             }
                             ra_table = rbind(ra1[["single_ra"]], ra2[["single_ra"]])
                             
                             setkey(ra_table, junc_match_idx, contig_idx)
                             ra_table[, f_side := private$LR_tbl[which(junc_sign == sign(junc_match) & seg_sign == sign(seg_match)), side], by = junc_match][, u_side := f_side * -1]

                             fused = mapply(function(full, c_match, j_match, loose_fuse) {
                                 if (!loose_fuse) {
                                     if (j_match < 0) {
                                         if (full[c_match] > 0) m = 1:(c_match) else m = (c_match):length(full)
                                     } else {
                                         if (full[c_match] < 0) m = 1:(c_match) else m = (c_match):length(full)
                                     }
                                     return(full[m])
                                 } else {
                                     full
                                 }
                             }, ra_table$contigs_lst, ra_table$cont_match, ra_table$junc_match, ra_table$will_fuse_loose, SIMPLIFY = F)

                             unfused = mapply(function(full, c_match, j_match, loose_fuse) {
                                 if (!loose_fuse) {
                                     if (j_match < 0) {
                                         if (full[c_match] < 0) m = 1:(c_match-1) else m = (c_match+1):length(full)
                                     } else {
                                         if (full[c_match] > 0) m = 1:(c_match-1) else m = (c_match+1):length(full)
                                     }
                                     return(full[m])
                                 } else {
                                     return(integer(0))
                                 }
                             }, ra_table$contigs_lst, ra_table$cont_match, ra_table$junc_match, ra_table$will_fuse_loose, SIMPLIFY = F)


                             
                             ra_table = cbind(ra_table, data.table(fused))
                             ra_table = cbind(ra_table, data.table(unfused))


                             new_nm = ra_table[,list(list(unique(unlist(contig_nm))))][,V1]

                             if (identical(ra_table$contig_idx[[1]], ra_table$contig_idx[[2]]) &
                                 identical(ra_table$contigs_lst[[1]], ra_table$contigs_lst[[2]]) &
                                 (! all(ra_table$cont_match == 1) & ! all(ra_table$cont_match == ra_table$contigs_len))) { ## this last one is for if matches happen on exact same contig at ends... BFB!

                                 contig = ra_table$contigs_lst[[1]]
                                 c_match = ra_table$cont_match
                                 ord = order(c_match)
                                 sorted_c_match = sort(c_match)
                                 ord_junc_match = ra_table$junc_match[ord]
                                 ord_seg_match = ra_table$seg_match[ord]
                                 split_contigs = private$junc_break(contig, ord_junc_match, sorted_c_match, ord_seg_match)
                                 resected_segs = private$junc_join(split_contigs, ord_junc_match, ord_seg_match)
                                 unfused_on_same_contig = list(resected_segs$unfused)
                                 new_contigs_lst = list(resected_segs$ra_contig)
    
                             } else {        
                                 f_side1 = ra_table[1,f_side]
                                 f_side2 = ra_table[2,f_side]
                                 
                                 fused1 = ra_table[1,fused[[1]]]
                                 fused2 = ra_table[2,fused[[1]]]
                                 new_contigs_lst = list(private$link_fused(f_side1, fused1, f_side2, fused2))

                             }

                             ra_contig = data.table(contigs_lst = new_contigs_lst,
                                                    contig_nm = new_nm,
                                                    contig_idx = prev_contigs[,max(contig_idx) + 1],
                                                    iter_idx = ra_table[, list(list(sort(c(unique(unlist(iter_idx)), i))))][,V1])
                             ra_contig[,c("loose_left", "loose_right") := list(sapply(contigs_lst, function(x) !(head(x, 1) %in% private$left_tel)), sapply(contigs_lst, function(x) !((tail(x, 1)) %in% private$right_tel)))]
                             


                             ##[, contig_idx := 1:.N]

                             
                             
                             if (identical(ra_table$contig_idx[[1]], ra_table$contig_idx[[2]]) &
                                 identical(ra_table$contigs_lst[[1]], ra_table$contigs_lst[[2]]) &
                                 (! all(ra_table$cont_match == 1) & ! all(ra_table$cont_match == ra_table$contigs_len))) {
                                 unfused_ra = ra_table[,list(contigs_lst = unfused_on_same_contig,
                                                contig_nm = contig_nm[1],
                                                contig_idx = contig_idx[1],
                                                iter_idx = list(sort(c(unlist(iter_idx[1]), i))))][,c("loose_left", "loose_right") := list(sapply(contigs_lst, function(x) !((head(x, 1)) %in% private$left_tel)), sapply(contigs_lst, function(x) !((tail(x, 1)) %in% private$right_tel)))]
                             } else {
                                 unfused_ra = ra_table[,list(contigs_lst = unfused,
                                                contig_nm,
                                                contig_idx,
                                                iter_idx = lapply(iter_idx, function(y) sort(c(y,i))), loose_left = sapply(unfused, function(x) !((head(x, 1)) %in% private$left_tel)), loose_right = sapply(unfused, function(x) !((tail(x, 1)) %in% private$right_tel)))]
                             }
                             
                             updated_contigs = rbind(ra_contig, prev_contigs[, n := NULL][!contig_idx %in% ra_table[,contig_idx]], unfused_ra)

                             private$epoch[[curr_step]] = updated_contigs
                             private$junction_history[[curr_step]] = junction_grl
                             private$j_segments_history[[curr_step]] = junc_match
                             private$ra_history[[curr_step]] = ra_table
                             private$curr_step = curr_step + 1

                         }
                         
                         private$current_index = private$current_index[-these_indices]

                         
                     }),
                 private = list(
                     junctions = NULL,
                     j_segments = NULL,
                     epoch = list(),
                     epoch_0 = NULL,
                     tile = NULL,
                     pos_tile = NULL,
                     scrambled_index = NULL,
                     current_index = NULL,
                     past_index = NULL,
                     curr_step = 1,
                     tot_steps = NULL,
                     discarded_junc = NULL,
                     left_tel = NULL,
                     right_tel = NULL,
                     ra_history = list(),
                     contigs_history = list(),
                     junction_history = GRangesList(),
                     j_segments_history = list(),
                     ## iterations = length(epoch),
                     LR_tbl = data.table(junc_sign = c(-1, -1, 1, 1),
                                         seg_sign = c(-1, 1, -1, 1),
                                         side = c(1, -1, -1, 1)),
                     link_fused = function(side1, fused1, side2, fused2) {
                         tmp1 = fused1
                         tmp2 = fused2
                         if (side1 > 0) {
                             tmp1 = -rev(tmp1)
                         }
                         if (side2 < 0) {
                             tmp2 = -rev(tmp2)
                         }
                             c(tmp1, tmp2)
                     },
                     list2dt = function(lst) {
                         num_elements = elementNROWS(lst)
                         .ix = 1:length(num_elements)
                         unlst = unlist(lst)
                         unlst_nm = names(unlst)
                         dt = data.table(unlisted = unlst, .ix = rep(.ix[num_elements>0], times = num_elements[num_elements > 0]))
                         dt[, .iix := 1:.N , by = .ix]
                         if (is.null(unlst_nm)) {
                             dt[, nm := NA]
                         } else {
                             dt[, nm := unlst_nm]
                         }     
                         return(dt)
                     },
                     contigs2dt = function(lst, tile) {
                         tmp_dt = private$list2dt(lst)
                         tmp_gr = tile[match(as.character(tmp_dt$unlisted),names(tile))]
                         tmp_dt = cbind(tmp_dt, gr2dt(tmp_gr)[, list(seqnames, start, end, strand, width)])
                         return(tmp_dt)
                     },
                     na2false = function(v) {
                         ## v = ifelse(is.na(v), v, FALSE)
                         v[is.na(v)] = FALSE
                         as.logical(v)
                     },
                     gv = function(gr, field) {
                         values(gr)[, field]
                     },
                     gr_keep_seq = function(gr, seqlev) {
                         
                         tmp = gr[seqnames(gr) %in% seqlev]
                         si = as.data.frame(seqinfo(tmp))
                         new_si = as(si[seqlev,], "Seqinfo")
                         ir = ranges(tmp)
                         seq = seqnames(tmp)
                         nm = names(tmp)
                         st = strand(tmp)
                         gr_val = values(tmp)
                         new_gr = GRanges(seq, ir, st, seqinfo = new_si)
                     },
                     map_junc2seg = function(query, subject) {
                         ## query is junction
                         ## subject is tiled
                         if (inherits(query, "GRangesList")) {
                             is_grl = TRUE
                             query = grl.unlist(query)
                         } else {
                             is_grl = FALSE
                         }
                         w = width(query)

                         if (any(w>1)) {
                             query = gr.start(query, width = 1)
                         }
                         
                         
                         query[strand(query) == "+"] = shift((query[strand(query) == "+"] ), 1)

                         tmp=as.data.table(as.data.frame(GenomicRanges::findOverlaps(query, subject, maxgap = 0, minoverlap = 1, type = "within", ignore.strand = T)))
                         return(tmp)
                         
                     },
                     karyograph = function (junctions, tile = NULL, label.edges = FALSE) 
                     {
                         require(gplots)
                         require(igraph)
                         require(Matrix)
                         require(data.table)
                         require(RColorBrewer)
                         if (length(junctions) > 0) {
                             bp.p = grl.pivot(junctions)
                             bp1 = gr.end(gr.fix(bp.p[[1]]), 1, ignore.strand = F)
                             bp2 = gr.start(gr.fix(bp.p[[2]]), 1, ignore.strand = F)
                             if (any(as.logical(strand(bp1) == "*") | as.logical(strand(bp2) == 
                                                                                 "*"))) {
                                 stop("Error: bp1 and bp2 must be signed intervals (i.e. either + or -)")
                             }
                             if (length(bp1) != length(bp2)) {
                                 stop("Error: bp1 and bp2 inputs must have identical lengths")
                             }
                             values(bp1)$bp.id = 1:length(bp1)
                             values(bp2)$bp.id = 1:length(bp1) + length(bp1)
                             pgrid = sgn1 = c(`-` = -1, `+` = 1)[as.character(strand(bp1))]
                             sgn2 = c(`-` = -1, `+` = 1)[as.character(strand(bp2))]
                             tmp.sl = seqlengths(grbind(bp1, bp2))
                             tmp.sl.og = tmp.sl
                             tmp.sl = gr2dt(grbind(bp1, bp2))[, max(end + 1, na.rm = TRUE), 
                                                              keyby = seqnames][names(tmp.sl.og), structure(pmax(V1, 
                                                                                                                 tmp.sl.og, na.rm = TRUE), names = names(tmp.sl.og))]
                             bp1 = gr.fix(bp1, tmp.sl)
                             bp2 = gr.fix(bp2, tmp.sl)
                         }
                         else {
                             if (is.null(tile)) {
                                 tile = si2gr(junctions)
                                 if (length(tile) == 0) {
                                     warning("Empty input given, producing empty output")
                                     return(NULL)
                                 }
                                 A = sparseMatrix(1, 1, x = 0, dims = rep(length(tile), 
                                                                          2))
                                 return(list(tile = tile, adj = A, G = graph.adjacency(A), 
                                             ab.adj = A != 0, ab.edges = NULL, junctions = junctions))
                             }
                             junctions = GRangesList()
                             bp1 = bp2 = GRanges()
                         }
                         if (!is.null(tile)) {
                             tile = gr.fix(tile)
                             tile = gr.fix(tile, bp1)
                             bp1 = gr.fix(bp1, tile)
                             strand(tile) = "+"
                             tile = disjoin(tile)
                             tile = sort(c(tile, gaps(tile)))
                             if (!identical(sort(seqlevels(tile)), seqlevels(junctions))) {
                                 tile = gr.fix(tile, junctions)
                                 junctions = gr.fix(junctions, tile)
                             }
                             if (length(junctions) > 0) {
                                 tbp = setdiff(gr.stripstrand(gr.trim(tile, 1)), gr.stripstrand(grbind(bp1, 
                                                                                                       bp2)))
                                 bp1 = gr.fix(bp1, tbp)
                                 bp2 = gr.fix(bp2, tbp)
                                 tbp = gr.fix(tbp, bp1)
                             }
                             else {
                                 tbp = gr.stripstrand(gr.trim(tile, 1))
                             }
                             tbp = tbp[start(tbp) != 1]
                             if (length(tbp) > 0) {
                                 tbp$seg.bp = TRUE
                             }
                         }
                         else {
                             tbp = NULL
                         }
                         if (length(junctions) > 0) {
                             if (length(tbp) > 0) {
                                 g = gaps(gr.stripstrand(sort(c(bp1[, c()], bp2[, 
                                                                                c()], tbp[, c()]))))
                             }
                             else {
                                 g = gaps(gr.stripstrand(sort(c(bp1[, c()], bp2[, 
                                                                                c()]))))
                             }
                         }
                         else {
                             g = gaps(gr.stripstrand(sort(tbp)))
                         }
                         g = g[strand(g) == "*"]
                         strand(g) = "+"
                         values(g)$bp.id = NA
                         values(g)$seg.bp = NA
                         tile.og = tile
                         tile = grbind(bp1, bp2, g, tbp)
                         tile = disjoin(gr.stripstrand(tile[order(gr.stripstrand(tile))]))
                         strand(tile) = "+"
                         tile = gr.fix(tile)
                         tile$is.tel = start(tile) == 1 | end(tile) == seqlengths(tile)[as.character(seqnames(tile))]
                         values(tile)$tile.id = 1:length(tile)
                         junc.bp = grbind(bp1, bp2)
                         junc.bpix = numeric()
                         if (length(junc.bp) > 0) {
                             junc.bpix = which(paste(seqnames(tile), end(tile)) %in% 
                                               paste(seqnames(junc.bp), start(junc.bp)))
                         }
                         tile = gr.fix(tile, bp1)
                         tile = gr.fix(tile, bp2)
                         bp1 = gr.fix(bp1, tile)
                         bp2 = gr.fix(bp2, tile)
                         all.bp = grbind(bp1, bp2, tbp)
                         all.bpix = numeric()
                         if (length(all.bp) > 0) {
                             all.bpix = which(paste(seqnames(tile), end(tile)) %in% 
                                              paste(seqnames(all.bp), start(all.bp)))
                         }
                         if (length(all.bpix > 0)) {
                             to.fuse = all.bpix[which(all.bpix > 1 & !((all.bpix - 
                                                                        1) %in% all.bpix))]
                             end(tile)[to.fuse - 1] = end(tile)[to.fuse - 1] + 1
                             tile = tile[-to.fuse]
                         }
                         if (length(junc.bpix) > 0) {
                             ab.pairs = cbind(ifelse(as.logical(strand(bp1) == "+"), 
                                                     gr.match(GenomicRanges::shift(gr.start(bp1), 1), 
                                                              gr.start(tile)), gr.match(gr.start(bp1), gr.end(tile))), 
                                              ifelse(as.logical(strand(bp2) == "+"), gr.match(GenomicRanges::shift(gr.start(bp2), 
                                                                                                                   1), gr.start(tile)), gr.match(gr.start(bp2), 
                                                                                                                                                 gr.end(tile))))
                             ab.pairs.bpid = bp1$bp.id
                             pp = (sgn1 * sgn2) > 0 & sgn1 > 0
                             mm = (sgn1 * sgn2) > 0 & sgn1 < 0
                             mp = sgn1 > 0 & sgn2 < 0
                             ab.pairs[pp, 1] = -ab.pairs[pp, 1]
                             ab.pairs[mm, 2] = -ab.pairs[mm, 2]
                             ab.pairs[mp, ] = -ab.pairs[mp, ]
                             edge.id = rep(1:nrow(ab.pairs), 2)
                             ab.pairs = rbind(ab.pairs, cbind(-ab.pairs[, 2], -ab.pairs[, 
                                                                                        1]))
                             ab.pairs.bpid = c(ab.pairs.bpid, ab.pairs.bpid)
                             adj.ab = Matrix(0, nrow = 2 * length(tile), ncol = 2 * 
                                                                             length(tile), dimnames = rep(list(as.character(c(1:length(tile), 
                                                                                                                              -(1:length(tile))))), 2))
                             tmp.ix = cbind(match(as.character(ab.pairs[, 1]), rownames(adj.ab)), 
                                            match(as.character(ab.pairs[, 2]), colnames(adj.ab)))
                             adj.ab[tmp.ix[!duplicated(tmp.ix), , drop = F]] = ab.pairs.bpid[!duplicated(tmp.ix)]
                         }
                         else {
                             ab.pairs.bpid = edge.id = c()
                             ab.pairs = matrix(nrow = 0, ncol = 2)
                             adj.ab = Matrix(FALSE, nrow = 2 * length(tile), ncol = 2 * 
                                                                                 length(tile), dimnames = rep(list(as.character(c(1:length(tile), 
                                                                                                                                  -(1:length(tile))))), 2))
                         }
                         seg.ix = 1:length(tile)
                         ref.pairs = cbind(seg.ix[1:(length(seg.ix) - 1)], seg.ix[2:(length(seg.ix))])
                         ref.pairs = ref.pairs[which(as.character(seqnames(tile[ref.pairs[, 
                                                                                          1]])) == as.character(seqnames(tile[ref.pairs[, 2]]))), 
                                               ]
                         if (nrow(ref.pairs) > 0) {
                             edge.id = c(edge.id, max(edge.id) + rep(1:nrow(ref.pairs), 
                                                                     2))
                             ref.pairs = rbind(ref.pairs, cbind(-ref.pairs[, 2], -ref.pairs[, 
                                                                                            1]))
                             adj.ref = Matrix(0, nrow = 2 * length(tile), ncol = 2 * 
                                                                              length(tile), dimnames = rep(list(as.character(c(1:length(tile), 
                                                                                                                               -(1:length(tile))))), 2))
                             adj.ref[cbind(match(as.character(ref.pairs[, 1]), rownames(adj.ref)), 
                                           match(as.character(ref.pairs[, 2]), colnames(adj.ref)))] = nrow(ab.pairs) + 
                                 1:nrow(ref.pairs)
                         }
                         else {
                             adj.ref = Matrix(FALSE, nrow = 2 * length(tile), ncol = 2 * 
                                                                                  length(tile), dimnames = rep(list(as.character(c(1:length(tile), 
                                                                                                                                   -(1:length(tile))))), 2))
                         }
                         tmp.nm = as.character(c(1:length(tile), -(1:length(tile))))
                         tile = c(tile, gr.flipstrand(tile))
                         names(tile) = tmp.nm
                         adj.source = sign(adj.ref) + 2 * sign(adj.ab)
                         adj = sign(adj.ref) + sign(adj.ab)
                         edges = Matrix::which(adj != 0, arr.ind = TRUE)
                         adj[edges] = 1:nrow(edges)
                         rownames(adj) = colnames(adj) = 1:nrow(adj)
                         G = graph.adjacency(adj, weighted = "edge.ix")
                         node.ind = abs(as.numeric(V(G)$name))
                         V(G)$chrom = as.character(seqnames(tile))[node.ind]
                         V(G)$start = start(tile)[node.ind]
                         V(G)$end = end(tile)[node.ind]
                         V(G)$width = width(tile)[node.ind]
                         V(G)$strand = sign(as.numeric(V(G)$name))
                         V(G)$size = 5
                         V(G)$shape = c("rectangle", "crectangle")[1 + as.numeric(V(G)$strand < 
                                                                                      0)]
                         V(G)$border.width = c(1, 2)[1 + as.numeric(V(G)$strand == 
                                                                        "-")]
                         V(G)$label = paste(V(G)$chrom, ":", round(V(G)$start/1e+06, 
                                                                   0), "-", round(V(G)$end/1e+06, 0), sep = "")
                         V(G)$label[V(G)$strand < 0] = paste(V(G)$chrom, ":", round(V(G)$end/1e+06, 
                                                                                    0), "-", round(V(G)$start/1e+06, 0), sep = "")[V(G)$strand < 
                                                                                                                                       0]
                         col.map = structure(brewer.master(length(seqlevels(tile))), 
                                             names = seqlevels(tile))
                         V(G)$chrom.ord = levapply(as.numeric(V(G)$start), list(V(G)$chrom), 
                                                   "rank")
                         V(G)$y = V(G)$chrom.ord * 30
                         V(G)$x = chr2num(V(G)$chrom) * 300 + 100 * rep(c(0, 1), each = length(tile)/2)
                         V(G)$col = col.map[V(G)$chrom]
                         E(G)$weight = 1
                         E(G)$from = edges[E(G)$edge.ix, 1]
                         E(G)$to = edges[E(G)$edge.ix, 2]
                         E(G)$col = c(col2hex("gray20"), col2hex("red"))[adj.source[edges[E(G)$edge.ix, 
                                                                                          ]]]
                         E(G)$type = c("reference", "aberrant", "aberrant")[adj.source[edges[E(G)$edge.ix, 
                                                                                             ]]]
                         E(G)$line.style = "SEPARATE_ARROW"
                         E(G)$arrow.shape = "ARROW"
                         E(G)$width = 1
                         ab.ix = E(G)$type == "aberrant"
                         E(G)$bp.id = NA
                         if (length(ab.pairs.bpid) > 0) {
                             E(G)$bp.id[ab.ix] = ab.pairs.bpid[adj.ab[cbind(E(G)$from[ab.ix], 
                                                                            E(G)$to[ab.ix])]]
                         }
                         E(G)$eid = NA
                         E(G)$eid[ab.ix] = edge.id[adj.ab[cbind(E(G)$from[ab.ix], 
                                                                E(G)$to[ab.ix])]]
                         E(G)$eid[!ab.ix] = edge.id[adj.ref[cbind(E(G)$from[!ab.ix], 
                                                                  E(G)$to[!ab.ix])]]
                         values(tile) = values(tile)[, c("tile.id", "is.tel")]
                         tile$ab.source = 1:length(tile) %in% E(G)$from[ab.ix]
                         tile$ab.target = 1:length(tile) %in% E(G)$to[ab.ix]
                         ab.edges = array(NA, dim = c(length(junctions), 3, 2), dimnames = list(NULL, 
                                                                                                c("from", "to", "edge.ix"), c("+", "-")))
                         dupped = duplicated(ab.pairs.bpid)
                         ab.edges[, 1:2, 1] = cbind(match(ab.pairs[!dupped, 1], names(tile)), 
                                                    match(ab.pairs[!dupped, 2], names(tile)))
                         ab.edges[, 1:2, 2] = cbind(match(ab.pairs[dupped, 1], names(tile)), 
                                                    match(ab.pairs[dupped, 2], names(tile)))
                         ab.edges[, 3, 1] = match(paste(ab.edges[, 1, 1], "|", ab.edges[, 
                                                                                        2, 1]), paste(E(G)$from, "|", E(G)$to))
                         ab.edges[, 3, 2] = match(paste(ab.edges[, 1, 1], "|", ab.edges[, 
                                                                                        2, 1]), paste(E(G)$from, "|", E(G)$to))
                         if (label.edges & nrow(ab.edges) > 0) {
                             ix = c(ab.edges[, 1, 1], ab.edges[, 2, 1], ab.edges[, 
                                                                                 1, 2], ab.edges[, 2, 2])
                             tile$edges.out = tile$edges.in = ""
                             tile$edges.in[ix] = sapply(ix, function(x) {
                                 ix = which(adj[, x] != 0)
                                 paste(ix, "->", sep = "", collapse = ",")
                             })
                             tile$edges.out[ix] = sapply(ix, function(x) {
                                 ix = which(adj[x, ] != 0)
                                 paste("->", ix, sep = "", collapse = ",")
                             })
                         }
                         return(list(tile = tile, adj = adj, G = G, ab.adj = adj.ab != 
                                                                        0, ab.edges = ab.edges, junctions = junctions))
                     },
                     junc_break = function(contig, ord_junc_match, sorted_c_match, ord_seg_match) {

                         contig_side_to_fuse = ifelse(sign(ord_junc_match) == sign(ord_seg_match), "right", "left")

                         where_fuse_loose = ifelse(sorted_c_match == 1, "fuse_loose_first",
                                            ifelse(sorted_c_match == length(contig), "fuse_loose_last",
                                                   NA))
                         

                         ## ints = mapply(function(junc, c_match, seg_match) {
                         ##     if (junc != seg_match) {
                         ##         m = 1:c_match
                         ##         return(m)
                         ##     } else {
                         ##         m = c_match:(length(contig))
                         ##     }
                         ## }, ord_junc_match, sorted_c_match, ord_seg_match, SIMPLIFY = F)



                         ir_tmp = unlist(IRangesList(mapply(function(junc, c_match, seg_match) {
                             if (junc != seg_match) {
                                 m = IRanges(1, c_match)
                                 return(m)
                             } else {
                                 m = IRanges(c_match, length(contig))
                             }
                         }, ord_junc_match, sorted_c_match, ord_seg_match)))


                         ## ## ord = order(c_match)
                         ## ## c_match = sort(c_match)
                         ## ## junc_match = junc_match[ord]

                         ## if (! all(sorted_c_match[1]==sorted_c_match[-1])) {
                         ##     i = 1
                         ##     tmp1 = ifelse(contig[sorted_c_match][i] > 0,
                         ##            ifelse(ord_junc_match[i] < 0, sorted_c_match[1] + 1, sorted_c_match[i] - 1),
                         ##            ifelse(ord_junc_match[i] > 0, sorted_c_match[1] + 1, sorted_c_match[i] - 1))

                         ##     i = -1
                         ##     tmp2 = ifelse(contig[sorted_c_match][i] > 0,
                         ##            ifelse(ord_junc_match[i] < 0, sorted_c_match[i], sorted_c_match[i] - 1),
                         ##            ifelse(ord_junc_match[i] > 0, sorted_c_match[i], sorted_c_match[i] - 1))

                         ##     tmp = c(tmp1, tmp2)


                         ##     first_match = tmp[1]
                         ##     second_match = tmp[2]
                         ## } else {
                         ##     first_match = second_match = unique(sorted_c_match)
                         ## }
                         
                         ir = disjoin(c(IRanges(1, length(contig)), ir_tmp))

                         ## ir = disjoin(c(IRanges(1, length(contig)), IRanges(first_match, second_match)))

                         if (length(ir) == 1) {
                             if (identical(contig_side_to_fuse, c("right", "left"))) {
                                 begin = integer(0)
                                 middle = as.integer(ir[1])
                                 end = integer(0)
                             }
                         } else if (length(ir) == 2) { ## i.e. this only happens when both breakpoints occur on the same contig and one of the two matches to a loose end
                             if (all(contig_side_to_fuse == "right")) {
                                 begin = integer(0)
                                 middle = as.integer(ir[1])
                                 end = as.integer(ir[2])
                             } else if (all(contig_side_to_fuse == "left")) {
                                 begin = as.integer(ir[1])
                                 middle = as.integer(ir[2])
                                 end = integer(0)
                             } else if (identical(contig_side_to_fuse, c("right", "left"))) {
                                 if (all(where_fuse_loose %in% c(NA, "fuse_loose_last"))) {
                                     begin = as.integer(ir[1])
                                     middle = as.integer(ir[2])
                                     end = integer(0)
                                 } else if (all(where_fuse_loose %in% c("fuse_loose_first", NA))) {
                                     begin = as.integer(0)
                                     middle = as.integer(ir[1])
                                     end = as.integer(ir[2])
                                 }
                             }
                             
                         } else if (length(ir) == 3) {
                             begin = as.integer(ir[1])
                             middle = as.integer(ir[2])
                             end = as.integer(ir[3])
                         }
                         
                         return(list(begin = contig[begin], middle = contig[middle], end = contig[end]))
                     },
                     junc_join = function(c_parts, ord_junc_match, ord_seg_match) {

                         directionality = ifelse(ord_junc_match == ord_seg_match, "right", "left")
                         
                         if (directionality[1] == directionality[2]) {
                             if (all(directionality == "left")) {
                                 new_join = c(c_parts$begin, -rev(c_parts$middle))
                                 unjoin = c_parts$end
                             } else {
                                 new_join = c(-rev(c_parts$middle), c_parts$end)
                                 unjoin = c_parts$begin
                             }
                         } else if (directionality[1] == "left") {
                             new_join = c(c_parts$begin, c_parts$end)
                             unjoin = c_parts$middle
                             ## unjoin = integer(0) ## changing the behavior of a deletion to remove the segments forever
                         } else {
                             new_join = c(c_parts$begin, rep(c_parts$middle, 2), c_parts$end)
                             unjoin = integer(0)
                         }
                         return(list(ra_contig = new_join, unfused = unjoin))
                     },
                     matching_algo2 = function(gr_match, ## matches contig within window
                                               junction_grl,
                                               junc_match,
                                               padding = 3e4) {

                         possible_matches1 = private$get_matches(gr_match, junction_grl, junc_match,  1, padding = padding)
                         possible_matches2 = private$get_matches(gr_match, junction_grl, junc_match,  2, padding = padding)

                         dist1 = distance(junction_grl[1], possible_matches1, ignore.strand = T)
                         exact1 = which(dist1 == 0)
                         far1 = which(dist1 > 0)
                         closest1 = grbind(possible_matches1[exact1], possible_matches1[far1][which.min(dist1[far1])])
                         dist2 = distance(junction_grl[2], possible_matches2, ignore.strand = T)
                         exact2 = which(dist2 == 0)
                         far2 = which(dist2 > 0)
                         closest2 = grbind(possible_matches2[exact2], possible_matches2[far2][which.min(dist1[far2])])

                         same_contig = intersect(closest1$.ix, closest2$.ix)

                         junc_dist = distance(junction_grl[1], junction_grl[2], ignore.strand = T)
                         strand1 = as.character(strand(junction_grl[1]))
                         strand2 = as.character(strand(junction_grl[2]))

                         if (!is.na(junc_dist)) {
                             if (strand1 == "+" & strand2 == "+") {
                                 boundary_end = max(start(closest2), start(junction_grl[2])) + 1
                                 boundary_start = start(junction_grl[1]) + 1
                                 boundary = GRanges(seqnames(junction_grl[1]), IRanges(boundary_start, boundary_end))
                             } else if (strand1 == "-" & strand2 == "-") {
                                 boundary_end = start(junction_grl[2])
                                 boundary_start = min(start(closest1), start(junction_grl[1]))
                                 boundary = GRanges(seqnames(junction_grl[1]), IRanges(boundary_start, boundary_end))
                             } else if (strand1 == "+" & strand2 == "-") {
                                 boundary_end = start(junction_grl[2])
                                 boundary_start = start(junction_grl[1]) + 1
                                 boundary = GRanges(seqnames(junction_grl[1]), IRanges(boundary_start, boundary_end))
                             } else if (strand1 == "-" & strand2 == "+") {
                                 boundary = GRanges(seqnames(junction_grl[1]), IRanges(-1e9, 1e9))
                             }
                         } else {
                             boundary = GRanges(seqlevels(junction_grl), IRanges(-1e9, 1e9))
                         }
                         
                         ## `%&%` = gUtils::`%&%` ## temporary fix
                         ## gUtils::`%&%`(closest2, boundary) ## consistent, but annoying fix
                         ## eligible1 = closest1 %&% boundary
                         ## eligible2 = closest2 %&% boundary

                         eligible1 = gUtils::`%&%`(closest1, boundary)
                         eligible2 = gUtils::`%&%`(closest2, boundary)

                         return(list(j1_matches = eligible1, j2_matches = eligible2))

                     },
                     select_match = function(possible_matches, which_junc, junc_match) {

                         junction_reference_face = ifelse(junc_match[which_junc] < 0, "left", "right")

                         possible_matches = possible_matches[[which_junc]]
                         
                         junc_match_field = paste0("junc_match", which_junc)
                         will_fuse_loose = private$gv(possible_matches, "will_fuse_loose")
                         ## below enforces fusion of loose end always: using loose end that fuses the most number of segments
                         num_exact_matches = length(which(abs(possible_matches$seg_id) == abs(private$gv(possible_matches, junc_match_field))))
                         
                         ## 
                         if (any(will_fuse_loose)) {        
                             if (length(which(will_fuse_loose)) > 1) {

                                 ## max_fused = max(possible_matches[possible_matches$will_fuse_loose]$number_of_segments_to_be_fused)
                                 ## max_fused = max(possible_matches[possible_matches$will_fuse_loose]$contig_size_to_be_fused) ## changing heuristic to maximum length and also not excluding non loose fusion
                                 ## tmp = possible_matches[possible_matches$number_of_segments_to_be_fused == max_fused]
                                 ## tmp = sample(tmp[tmp$contig_size_to_be_fused == max(tmp$contig_size_to_be_fused)], 1)
                                 
                                 selection = (possible_matches$will_fuse_loose & possible_matches$num_rearrangements_to_be_fused > 0)

                                 if (any(selection)) {
                                     tmp = possible_matches[selection]
                                 } else {
                                     max_fused = max(possible_matches$num_rearrangements_to_be_fused) ## changing heuristic to maximum length and also not excluding non loose fusion
                                     tmp = possible_matches[possible_matches$num_rearrangements_to_be_fused == max_fused]
                                 }

                                 tmp = sample(tmp[tmp$num_rearrangements_to_be_fused == max(tmp$num_rearrangements_to_be_fused)], 1)

                                 final_match = tmp
                                 sign_seg = sign(private$gv(final_match, "seg_id"))
                             } else {
                                 ## max_fused = max(possible_matches$contig_size_to_be_fused) ## added this in

                                 ## tmp = possible_matches[possible_matches$will_fuse_loose]
                                 ## tmp = possible_matches[possible_matches$contig_size_to_be_fused == max_fused] ## added this in

                                 ## tmp = sample(tmp[tmp$contig_size_to_be_fused == max(tmp$contig_size_to_be_fused)], 1) ## added this in


                                 selection = (possible_matches$will_fuse_loose & possible_matches$num_rearrangements_to_be_fused > 0)

                                 if (any(selection)) {
                                     tmp = possible_matches[selection]
                                 } else {
                                     max_fused = max(possible_matches$num_rearrangements_to_be_fused) ## changing heuristic to maximum length and also not excluding non loose fusion
                                     tmp = possible_matches[possible_matches$num_rearrangements_to_be_fused == max_fused]
                                 }

                                 tmp = sample(tmp[tmp$num_rearrangements_to_be_fused == max(tmp$num_rearrangements_to_be_fused)], 1)

                                 ## max_fused = max(possible_matches$num_rearrangements_to_be_fused) ## added this in
                                 ## tmp = possible_matches[possible_matches$num_rearrangements_to_be_fused == max_fused] ## added this in
                                 ## tmp = sample(tmp[tmp$num_rearrangements_to_be_fused == max(tmp$num_rearrangements_to_be_fused)], 1) ## added this in
                                 final_match = tmp
                                 sign_seg = sign(private$gv(final_match, "seg_id"))
                             }
                             
                             if (junction_reference_face == "left") {
                                 ## if (private$gv(final_match, junc_contig_face_field) == "right") {
                                 if (abs(private$gv(final_match, "seg_id")) < abs(private$gv(final_match, junc_match_field))) {
                                     fill_in_segs = seq(abs(private$gv(final_match, "seg_id")) + 1, abs(private$gv( final_match, junc_match_field))) * sign_seg
                                     segs_to_trim = NULL
                                     where_to_fill = private$gv(final_match, ".iix")
                                     where_to_trim = NULL
                                 } else if (abs(private$gv(final_match, "seg_id")) > abs(private$gv(final_match, junc_match_field))) {
                                     fill_in_segs = NULL
                                     segs_to_trim = seq(abs(private$gv( final_match, junc_match_field)) + 1, abs(private$gv(final_match, "seg_id"))) * sign_seg
                                     where_to_fill = NULL
                                     segment_diff = abs(private$gv(final_match, "seg_id")) - abs(private$gv(final_match, junc_match_field))
                                     where_to_trim = (private$gv(final_match, ".iix")-segment_diff+1):private$gv(final_match, ".iix")
                                 } else {
                                     fill_in_segs = NULL
                                     segs_to_trim = NULL
                                     where_to_fill = NULL
                                     where_to_trim = NULL
                                 }
                             }
                             else if (junction_reference_face == "right") {
                                 ## else if ((private$gv(final_match, junc_contig_face_field)) == "left") {
                                 if (abs(private$gv(final_match, "seg_id")) < abs(private$gv(final_match, junc_match_field))) {
                                     segs_to_trim = seq(abs(private$gv(final_match, "seg_id")), abs(private$gv( final_match, junc_match_field)) - 1) * sign_seg
                                     fill_in_segs = NULL
                                     segment_diff = abs(private$gv(final_match, junc_match_field)) - abs(private$gv(final_match, "seg_id"))
                                     where_to_trim = 1:segment_diff
                                     where_to_fill = NULL
                                 } else if (abs(private$gv(final_match, "seg_id")) > abs(private$gv(final_match, junc_match_field))){
                                     fill_in_segs = seq(abs(private$gv(final_match, junc_match_field)), abs(private$gv( final_match, "seg_id")) - 1) * sign_seg
                                     segs_to_trim = NULL
                                     where_to_fill = private$gv(final_match, ".iix")
                                     where_to_trim = NULL
                                 } else {
                                     fill_in_segs = NULL
                                     segs_to_trim = NULL
                                     where_to_fill = NULL
                                     where_to_trim = NULL
                                 }
                             }
                         } else { 
                             if (num_exact_matches > 1) {
                                 ## max_fused = max(possible_matches$number_of_segments_to_be_fused)
                                 ## tmp = possible_matches[possible_matches$number_of_segments_to_be_fused == max_fused]
                                 max_fused = max(possible_matches$contig_size_to_be_fused) ## added this in
                                 tmp = possible_matches[possible_matches$contig_size_to_be_fused == max_fused] ## added this in
                                 tmp = sample(tmp[tmp$contig_size_to_be_fused == max(tmp$contig_size_to_be_fused)], 1)
                                 final_match = tmp
                             }
                             else if (length(possible_matches) == 1) {
                                 tmp = possible_matches
                                 final_match = tmp
                             }
                             fill_in_segs = NULL
                             segs_to_trim = NULL
                             where_to_fill = NULL
                             where_to_trim = NULL
                         }
                         ## where_to_add = ifelse(sign(private$gv(final_match, "seg_id")) == sign(private$gv(final_match, junc_match_field)),
                         ##                       private$gv(final_match, "num_from_right"),
                         ##                       private$gv(final_match, "num_from_left"))
                         where_to_add = private$gv(final_match, "num_from_left")
                         values(final_match)[, "where_to_add"] = where_to_add


                         return(list(final_match = final_match,
                                     fill_in_segs = fill_in_segs,
                                     segs_to_trim = segs_to_trim,
                                     where_to_fill = where_to_fill,
                                     where_to_trim = where_to_trim))

                     },
                     get_matches = function(gr_match,
                                            junction_grl,
                                            junc_match,
                                            which_breakpoint_ix,
                                            padding = 3e4)
                     {
                         this_id = which_breakpoint_ix
                         junction_reference_face = ifelse(junc_match[which_breakpoint_ix] < 0, "left", "right")
                         junc_match_field = paste0("junc_match", this_id)
                         junc_contig_face_field = paste0("junc_contig_face", this_id)
                         this_junc = junction_grl[this_id]
                         this_junc_match = junc_match[this_id]
                         ## junc_match < 0 means on REFERENCE, junction is left facing, 
                         ## this is so any match does not have its reference "end" (right side if oriented on reference)
                         ## beyond the window of query (30kb for now?)
                         this_win = padding
                         if (this_win > 1e9) {
                             this_win = pmin(this_win, 1e9)
                         }
                         junction_window = (this_junc+this_win)
                         if (junction_reference_face == "left") {
                             ## no_further = end(gr_match) <= end(junction_window)
                             no_further = end(gr_match) <= end(this_junc)
                         } else {
                             ## no_further = start(gr_match) >= start(junction_window)
                             no_further = start(gr_match) >= start(this_junc)
                         }
                         
                         possible_matches = gr_match[gr_match %^% junction_window & no_further]
                         values(possible_matches)[, "junc_match_idx"] = which_breakpoint_ix


                         

                         ## grabbing only matches that are loose reference facing to the right or that actually match to junction
                         if (junc_match[this_id] < 0)  {## if junction is left facing ON REFERENCE

                             sel = private$na2false(possible_matches$loose_left_reference_face == "right") |
                                 private$na2false(possible_matches$loose_right_reference_face == "right") |
                                 !is.na(values(possible_matches[, junc_match_field])[,1])
                                        
                             possible_matches = possible_matches[sel]
                         } else if (junc_match[this_id] > 0) {

                              sel = private$na2false(possible_matches$loose_left_reference_face == "left") |
                                 private$na2false(possible_matches$loose_right_reference_face == "left") |
                                 !is.na(values(possible_matches[, junc_match_field])[,1])
                                        
                             possible_matches = possible_matches[sel]
                         }
                         

                         ## values(possible_matches)[, junc_match_field] = na.omit(private$gv(possible_matches,junc_match_field ))
                         values(possible_matches)[, junc_match_field] = this_junc_match

                         ## junction contig face also follows same sign convention as if the junction actually matches to the contig
                         values(possible_matches)[, junc_contig_face_field] = ifelse(sign(private$gv(possible_matches, junc_match_field)) == sign(possible_matches$seg_id), "right", "left")
                         

                         ## now how to select?
                         ## if there is loose end, always grab that one?
                         ## or do selection based on # segs and/or contig distance?

                         ## collecting stats for possible matches

                         contig_side_to_fuse = ifelse(sign(private$gv(possible_matches, "seg_id")) == sign(private$gv(possible_matches, junc_match_field)),
                                                      "right",
                                                      "left")
                         
                         values(possible_matches)[, "contig_side_to_fuse"] = contig_side_to_fuse

                         num_exact_matches = length(which(abs(possible_matches$seg_id) == abs(private$gv(possible_matches, junc_match_field))))

                         ## will_fuse_loose = !is.na(possible_matches$loose_right_reference_face) | !is.na(possible_matches$loose_left_reference_face)
                         if (this_junc_match < 0) {
                             will_fuse_loose = private$na2false(possible_matches$loose_right_reference_face == "right" & this_junc_match < 0 | possible_matches$loose_left_reference_face == "right" & this_junc_match < 0)
                         } else if (this_junc_match > 0) {
                             will_fuse_loose = private$na2false(possible_matches$loose_right_reference_face == "left" & this_junc_match > 0 | possible_matches$loose_left_reference_face == "left" & this_junc_match > 0)
                         }
                         
                         
                         possible_matches$will_fuse_loose = will_fuse_loose
                         number_of_segments_to_be_fused = ifelse(private$gv(possible_matches, junc_contig_face_field) == "left",
                                                                 possible_matches$num_from_left,
                                                                 possible_matches$num_from_right)
                         values(possible_matches)[,"number_of_segments_to_be_fused"] = number_of_segments_to_be_fused

                         contig_size_to_be_fused = ifelse(private$gv(possible_matches, junc_contig_face_field) == "left",
                                                          possible_matches$cumdist_from_left,
                                                          possible_matches$cumdist_from_right)
                         possible_matches$contig_size_to_be_fused = contig_size_to_be_fused

                         num_rearrangements_to_be_fused = ifelse(private$gv(possible_matches, junc_contig_face_field) == "left",
                                                                 possible_matches$num_rearrangements_from_left,
                                                                 possible_matches$num_rearrangements_from_right)
                         possible_matches$num_rearrangements_to_be_fused = num_rearrangements_to_be_fused

                         num_rearrangements_on_unfused = ifelse(private$gv(possible_matches, junc_contig_face_field) == "left",
                                                                possible_matches$num_rearrangements_from_right,
                                                                possible_matches$num_rearrangements_from_left)
                         possible_matches$num_rearrangements_on_unfused = num_rearrangements_on_unfused
                         return(possible_matches)
                     },
                     process_matching_algo = function(match_lst,
                                                      contigs_dt,
                                                      which_breakpoint_ix,
                                                      aggregated_contigs_tbl) {

                         
                         final_match = match_lst[['final_match']]
                         ## which_pool = private$gv(final_match, "pool")
                         matched_ix = private$gv(final_match, ".ix")
                         junc_match_idx = private$gv(final_match, "junc_match_idx")
                         contigs_dt = contigs_dt[order(.ix, .iix)] ## back to somatic ordering    
                         matched = contigs_dt[.ix == matched_ix]
                         junc_match_field = paste0("junc_match", which_breakpoint_ix)
                         ## junc_match = na.omit(c(private$gv(final_match, "junc_match1"), private$gv(final_match, "junc_match2")))
                         junc_match = private$gv(final_match, junc_match_field)
                         where_to_add = private$gv(final_match, "where_to_add")
                         select_from_agg = private$gv(final_match, "contig_idx")
                         ## cont_match = grab_contig_match(match_lst)
                         fill_in_segs = match_lst[["fill_in_segs"]]
                         which_side = private$gv(final_match, "contig_side_to_fuse")
                         will_fuse_loose = private$gv(final_match, "will_fuse_loose")
                         contigs_len = private$gv(final_match, "num_segs")
                         segs_to_trim = match_lst[["segs_to_trim"]]
                         where_to_trim = match_lst[["where_to_trim"]]
                         where_to_fill = match_lst[["where_to_fill"]]

                         ## in_current = private$gv(final_match, "pool") == "current"
                         ## in_hold_out = private$gv(final_match, "pool") == "hold_out"

                         

                         adjust_cont_match = ifelse(is.null(fill_in_segs), 0, length(fill_in_segs))

                         if (will_fuse_loose) {
                             ## cont_match = ifelse(which_side == "left", where_to_add + adjust_cont_match, 1)
                             cont_match = ifelse(which_side == "left", contigs_len + adjust_cont_match, 1)
                             contigs_len = contigs_len +  adjust_cont_match
                         } else {
                             ## cont_match = where_to_add
                             cont_match = private$gv(final_match, "num_from_left")
                         }

                         contigs_dt = contigs_dt[order(.ix, .iix)] ## back to somatic ordering    
                         matched = contigs_dt[.ix == matched_ix]

                         allelic_nm = strsplit(matched[.iix == where_to_add, nm], "_")[[1]][1]

                         if (!is.null(fill_in_segs)) fill_in_segs = setNames(fill_in_segs, paste(allelic_nm, sub("\\-", "", fill_in_segs), sep = "_"))

                         matched_contig = original_contig = setNames(matched[, seg_id], matched[, nm])
                         original_indices = 1:length(original_contig)
                         

                         ## if (private$gv(final_match, "contig_side_to_fuse") == "left") {
                         if (junc_match > 0) {
                             if (private$gv(final_match, "contig_side_to_fuse") == "left") {
                                 if (!is.null(fill_in_segs)) matched_contig = c(matched_contig, rev(fill_in_segs))
                                 new_indices = original_indices
                                 add_to_last = seq((original_indices[length(original_indices)] + 1), length.out = length(fill_in_segs))
                                 add_to_first = NULL
                             } else {
                                 if (!is.null(fill_in_segs)) matched_contig = c(fill_in_segs, matched_contig)
                                 new_indices = original_indices + length(fill_in_segs)
                                 add_to_first = seq(1, length.out = length(fill_in_segs))
                                 add_to_last = NULL
                             }
                         } else if (junc_match < 0) {
                             if (private$gv(final_match, "contig_side_to_fuse") == "left") {
                                 if (!is.null(fill_in_segs)) matched_contig = c(matched_contig, fill_in_segs)
                                 new_indices = original_indices
                                 add_to_last = seq((original_indices[length(original_indices)] + 1), length.out = length(fill_in_segs))
                                 add_to_first = NULL
                             } else {
                                 if (!is.null(fill_in_segs)) matched_contig = c(rev(fill_in_segs), matched_contig)
                                 new_indices = original_indices + length(fill_in_segs)
                                 add_to_first = seq(1, length.out = length(fill_in_segs))
                                 add_to_last = NULL
                             }
                         }


                         ## ra =  aggregated_contigs_tbl[contig_idx == select_from_agg & pool == which_pool]
                         ra =  aggregated_contigs_tbl[contig_idx == select_from_agg]
                         ## ra[, pool := NULL][, n := NULL][, c("in_contig",
                         ##                                     "junc_match_idx",
                         ##                                     "junc_match",
                         ##                                     "cont_match")
                         ##                                 := list(TRUE,junc_match_idx, junc_match, cont_match)]
                         ra[, n := NULL][, c("in_contig",
                                             "junc_match_idx",
                                             "junc_match",
                                             "cont_match")
                                         := list(TRUE,junc_match_idx, junc_match, cont_match)]
                         ra[, contigs_lst := list(list(matched_contig))]
                         ra[, seg_match := contigs_lst[[1]][cont_match], by = list(contig_idx,junc_match)]
                         ra[, contigs_len := contigs_len]
                         ra[, will_fuse_loose := will_fuse_loose]
                         return(list(single_ra  = ra,
                                     ## in_current = in_current,
                                     ## in_hold_out = in_hold_out,
                                     adjusted_indices = new_indices,
                                     add_to_first = add_to_first,
                                     add_to_last = add_to_last))
                         
                     }),
                 active = list(
                     ret_junctions = function(i = NULL) {
                         if (is.null(i)) {
                             return(private$junctions)
                         } else {
                             return(private$junctions[i])
                         }
                     },
                     junction_segments = function() {
                         return(private$j_segments)
                     },
                     ret_epoch = function(i = NULL) {
                         if (is.null(i)) {
                             return(private$epoch)
                         } else {
                             return(private$epoch[i])
                         }
                     },
                     ret_epoch0 = function() {
                         return(private$epoch_0)
                     },
                     ## ret_iterations = function(i = NULL) {
                     ##     if (is.null(i)) {
                     ##         return(private$epoch)
                     ##     } else {
                     ##         return(private$epoch[i])
                     ##     }
                     ## },
                     discarded = function() {
                         return(private$discarded_junc)
                     },
                     scrambled_ix = function() {
                         return(private$scrambled_index)
                     },
                     current_ix = function() {
                         return(private$current_index)
                     },
                     past_ix = function() {
                         return(private$past_index)
                     },
                     total_steps = function() {
                         return(private$tot_steps)
                     },
                     current_step = function() {
                         message("Total steps: ", private$tot_steps)
                         return(private$curr_step)
                     }
                 )
                 )

                     
