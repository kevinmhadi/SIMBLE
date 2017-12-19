sim_jab = function(junc, padding = 3e4, seed = NULL) {

    ## browser()
    
    kag = karyograph(junc)

    tile =  karyograph(junc)$tile
    seqlevels(tile) = seqlevels(tile)[orderSeqlevels(seqlevels(tile), TRUE)]
    tile = sort(tile)
    pos_tile = tile %Q% (strand == "+")
    new_tile = c(pos_tile, pos_tile)
    names(new_tile) = c(1:length(pos_tile), -(1:length(pos_tile)))
    strand(new_tile[(length(new_tile)/2+1):length(new_tile)]) = "-"

    tile = new_tile
    pos_tile = tile %Q% (strand == "+")
   
    u_junc = grl.unlist(junc)
    seqlevels(u_junc) = seqlevels(u_junc)[orderSeqlevels(seqlevels(u_junc), TRUE)]
    ## pos = as.logical(strand(u_junc) == "+")
    ##u_junc[pos] = shift(u_junc[pos], 1)
    
    
    sort_u_junc = sort(u_junc, ignore.strand = T)
    sort_u_junc$idx = 1:length(u_junc)

    segged = dropSeqlevels(pos_tile, setdiff(seqlevelsInUse(pos_tile), seqlevelsInUse(sort_u_junc)))
    seqlevels(segged) = seqlevels(segged)[orderSeqlevels(seqlevels(segged), TRUE)]
    segged = sort(segged)
    segged$id = seg_int = 1:length(segged)
    ## names(seg_int) = cbind(seg_int, as.data.frame("+"))[,2]
    ##browser()

    ref_contigs = split(segged, seqnames(segged), drop = F)
    ref_int_contigs = contigs_lst = as.list(split(seg_int, seqnames(segged)))
    tmp_tel_1 = unlist(lapply(ref_int_contigs, function(x) head(x,1)), use.names = F)
    tmp_tel_2 = unlist(lapply(ref_int_contigs, function(x) tail(x,1)), use.names = F)
    left_tel = c(tmp_tel_1, -tmp_tel_2)
    right_tel = c(-tmp_tel_1, tmp_tel_2)
    tel_seg =  unlist(lapply(ref_int_contigs, function(x) c(head(x,1), tail(x, 1))), use.names = F)
    rev_contigs  = lapply(contigs_lst, function(x) -rev(x))

    seg_idx = nearest(sort_u_junc, segged, ignore.strand = T)
    names(seg_idx) = as.character(strand(sort_u_junc))

    sort_u_junc$seg_idx = seg_idx

    dt_junc = gr2dt(sort_u_junc[,c("grl.ix", "grl.iix", "seg_idx")])
    connections = dt_junc[,list(seqnames, seg_idx, strand), by = grl.ix]

    lst = split(connections, connections$grl.ix)
    bint_lst = lapply(lst, function(x) {
        tmp = x[,seg_idx] * (as.numeric(x[,strand])*2-3)
        names(tmp) = as.character(x[,strand])
        return(tmp)
    })

    LR_tbl = data.table(junc_sign = c(-1, -1, 1, 1), seg_sign = c(-1, 1, -1, 1), side = c(1, -1, -1, 1))

    junc_pool = 1:length(junc)
    
    iteration_list = NULL

    hold_out = data.table(contigs_lst = list(), contig_nm = character(), contig_idx = integer(), iter_idx = integer(), loose_left = logical(), loose_right = logical())

    ## Utility functions ##
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
    }

    rename_cont_lst = function(subcont, nms) {
        mapply(function(subcont, nms) { names(subcont)=nms; return(subcont) }, subcont, nms, SIMPLIFY = FALSE)
    }
    
    
    autosomal_hi = sapply(contigs_lst[! names(contigs_lst) %in% c("X", "Y")], function(x) {names(x) = paste0(rep("hi", length(x)), "_", x); return(x)})
    autosomal_lo = sapply(contigs_lst[! names(contigs_lst) %in% c("X", "Y")], function(x) {names(x) = paste0(rep("lo", length(x)), '_', x); return(x)})
    
    diploid_autosomal = c(autosomal_hi, autosomal_lo)
    current_contigs = data.table(contigs_lst = diploid_autosomal, contig_nm = names(diploid_autosomal), contig_idx = seq_along(diploid_autosomal), iter_idx = 0, loose_left = FALSE, loose_right = FALSE) 

    
    ## current_contigs = data.table(contigs_lst, contig_nm = names(contigs_lst), contig_idx = seq_along(contigs_lst), iter_idx = 0, loose_left = FALSE, loose_right = FALSE) ## for haploid

    tmp_seg_1 = contigs_lst$X
    names(tmp_seg_1) = paste0(rep("hi", length(tmp_seg_1)), "_", tmp_seg_1); tmp_seg_1 = list(X = tmp_seg_1)
    
    if (!"Y" %in% seqlevelsInUse(sort_u_junc)) {
        if (rbinom(1, 1, prob = 0.5)) {
            tmp_seg_2 = unlist(tmp_seg_1); names(tmp_seg_2) = paste0(rep("lo", length(tmp_seg_2)), "_", tmp_seg_2); tmp_seg_2 = list(X = tmp_seg_2)
        } else {
            tmp_seg_2 = contigs_lst$Y; names(tmp_seg_2) = paste0(rep("lo", length(tmp_seg_2)), "_", tmp_seg_2); tmp_seg_2 = list(Y = tmp_seg_2)
        }
    } else {
        tmp_seg_2 = contigs_lst$Y; names(tmp_seg_2) = paste0(rep("lo", length(tmp_seg_2)), '_', tmp_seg_2); tmp_seg_2 = list(Y = tmp_seg_2)
    }
    
    sex_chrom = c(tmp_seg_1, tmp_seg_2)
    sex_chrom_dt = data.table(contigs_lst = sex_chrom, contig_nm = names(sex_chrom), contig_idx = c(45, 46), iter_idx = 0, loose_left = FALSE, loose_right = FALSE)

    current_contigs = rbind(current_contigs, sex_chrom_dt)
    ## browser()
        
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

    for (i in 1:length(junc)) {
        message("commencing iteration: ", i)

        
        rand_id = ifelse(length(junc_pool)>1, sample(c(junc_pool), 1), junc_pool)
        sample_junc = bint_lst[[rand_id]]
        junc_match = ifelse(sample_junc>0, sample_junc+1, sample_junc)
        junction_grl = junc[[rand_id]]





        ## new searching using GRanges to get window of contigs!!!
        ## this is to start the search for contigs

        agg = rbind(copy(current_contigs)[, pool := "current"], copy(hold_out)[, pool := "hold_out"]); agg = agg[elementNROWS(agg$contigs_lst) > 0,]; agg[, n := as.character(1:.N)]; agg[, c('loose_left', 'loose_right') := list(unlist(loose_left), unlist(loose_right))]

        tmp_match = contigs2dt(agg$contigs_lst, tile); setnames(tmp_match, "unlisted", "seg_id")

        tmp_match[, loose_left := FALSE]; tmp_match[, loose_right := FALSE];setkey(agg, n)

        tmp_match[, contig_idx := agg[as.character(.ix), contig_idx]]
        tmp_match[, pool := agg[as.character(.ix), pool]]
        tmp_match[, loose_left := ifelse(.iix == .iix[1], agg[as.character(.ix), loose_left], FALSE), by = .ix]
        tmp_match[, loose_right := ifelse(.iix == .iix[.N], agg[as.character(.ix), loose_right], FALSE), by = .ix]
        tmp_match[, loose_left_reference_face := ifelse(loose_left, ifelse(strand == "+", "left", "right"), NA)]
        tmp_match[, loose_right_reference_face := ifelse(loose_right, ifelse(strand == "+", "right", "left"), NA)]

        tmp_match[, is_terminal := .iix %in% c(.iix[1], .iix[.N]), by = .ix]

        tmp_match[, seqnames := factor(seqnames, c(as.character(1:24), "X", "Y", "M"))] ## necessary for later to order by reference distance

        tmp_match[, junc_match1 := ifelse(abs(junc_match[1]) == abs(seg_id), junc_match[1], NA)];tmp_match[, junc_match2 := ifelse(abs(junc_match[2]) == abs(seg_id), junc_match[2], NA)]

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


        gr_match = dt2gr(tmp_match)

        
        ## now get match within a certain window for ra1
        ## browser(expr = (i == 93))
        
        all_matches = matching_algo2(gr_match, junction_grl, junc_match, padding = padding)
        ## browser(expr = (i >= 198))
        
        ## browser(expr = (any(as.character(seqnames(junction_grl)) == "X") & i >= 200))
        
        match_lst1 = select_match(all_matches, 1, junc_match)
        match_lst2 = select_match(all_matches, 2, junc_match)

        is_on_same_contig = match_lst1$final_match$.ix == match_lst2$final_match$.ix



        
        

        ra1 = process_matching_algo(match_lst1, tmp_match, 1, agg)
        ra2 = process_matching_algo(match_lst2, tmp_match, 2, agg)

        combined_contigs = list(ra1$single_ra$contigs_lst[[1]], ra2$single_ra$contigs_lst[[1]])
        combined_indices = list(ra1$adjusted_indices, ra2$adjusted_indices)
        combined_matches = c(ra1$single_ra$cont_match, ra2$single_ra$cont_match)
        
        in_current = ra1$in_current | ra2$in_current
        in_hold_out = ra1$in_hold_out | ra2$in_hold_out



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

            ensure_correct_contig_matches = ifelse(!is.na(new_loose_index), ifelse(new_loose_index == "first", 1, length(tmp_contig)), combined_matches)
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
            
            

            
            
            
        
        ## adjust ra2 by ra1 indices
        


        setkey(ra_table, junc_match_idx, contig_idx)

        
        ra_table[, f_side := LR_tbl[which(junc_sign == sign(junc_match) & seg_sign == sign(seg_match)), side], by = junc_match][, u_side := f_side * -1]


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


        


        
        ## ra_contig = data.table(contigs_lst = list(link_fused(f_side1, fused1, f_side2, fused2)), contig_nm = paste0('junc_', rand_id), contig_idx = current_contigs[,max(contig_idx) + 1])

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

            split_contigs = junc_break(contig, ord_junc_match, sorted_c_match, ord_seg_match)

            resected_segs = junc_join(split_contigs, ord_junc_match, ord_seg_match)

            unfused_on_same_contig = list(resected_segs$unfused)

            new_contigs_lst = list(resected_segs$ra_contig)
            
        } else {
            f_side1 = ra_table[1,f_side]
            f_side2 = ra_table[2,f_side]
            
            fused1 = ra_table[1,fused[[1]]]
            fused2 = ra_table[2,fused[[1]]]
            new_contigs_lst = list(link_fused(f_side1, fused1, f_side2, fused2))
        }

        ra_contig = data.table(contigs_lst = new_contigs_lst,
                               contig_nm = new_nm,
                               contig_idx = current_contigs[,max(contig_idx) + 1],
                               iter_idx = ra_table[, list(list(sort(c(unique(unlist(iter_idx)), i))))][,V1])
        ra_contig[,c("loose_left", "loose_right") := list(sapply(contigs_lst, function(x) !(head(x, 1) %in% left_tel)), sapply(contigs_lst, function(x) !((tail(x, 1)) %in% right_tel)))]
            


        updated_contigs = rbind(ra_contig, current_contigs[!contig_idx %in% ra_table[,contig_idx]])##[, contig_idx := 1:.N]

                

        if (in_hold_out) {
            updated_hold_out = hold_out[!contig_idx %in% ra_table[,contig_idx]]
        } else {
            updated_hold_out = hold_out
        }
        

        

        contigs_lst = updated_contigs$contigs
        ##names(contigs_lst) = updated_contigs$contig_nm

        not_to_throw_out = ra_contig$contigs_lst[[1]]
        
        if (ra_table[,identical(contigs_lst[[1]],contigs_lst[[2]])]) {
            updated_hold_out = rbind(updated_hold_out, ra_table[,list(contigs_lst = unfused_on_same_contig, contig_nm = contig_nm[1], contig_idx = contig_idx[1], iter_idx = list(sort(c(unlist(iter_idx[1]), i))))][,c("loose_left", "loose_right") := list(sapply(contigs_lst, function(x) !((head(x, 1)) %in% left_tel)), sapply(contigs_lst, function(x) !((tail(x, 1)) %in% right_tel)))])
        } else {
            updated_hold_out = rbind(updated_hold_out, ra_table[,list(contigs_lst = unfused, contig_nm, contig_idx, iter_idx = lapply(iter_idx, function(y) sort(c(y,i))), loose_left = sapply(unfused, function(x) !((head(x, 1)) %in% left_tel)), loose_right = sapply(unfused, function(x) !((tail(x, 1)) %in% right_tel)))])
        }



        ## contigs in hold_out need marked as having a rearrangement on them which sticks even after further perturbations


        iteration_list = c(iteration_list, list(list(junc_idx = rand_id,
                                                     junc_match = junc_match,
                                                     input_contigs = current_contigs,
                                                     output_contigs = updated_contigs,
                                                     ra_table = ra_table,
                                                     input_hold_out = hold_out,
                                                     output_hold_out = updated_hold_out,
                                                     in_current = in_current,
                                                     in_hold_out = in_hold_out,
                                                     junction_grl = junc[rand_id])))
        current_contigs = updated_contigs

        hold_out = updated_hold_out

        junc_pool = junc_pool[-match(rand_id, junc_pool)]





        message("completed iteration ", i, "\n")
    }
    
    return(iteration_list)
}
