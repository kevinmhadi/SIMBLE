## contig = ra_table$contigs_lst[[1]]
## c_match = ra_table$cont_match

errnot = function(...) {
    return(tryCatch(..., error = function(e) NULL))
}

unlapply = function(...) {
    return(unlist(lapply(...)))
}

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
    
        
        

    

    



    ## print(contig[begin])
    ## print(contig[middle])
    ## print(contig[end])

    return(list(begin = contig[begin], middle = contig[middle], end = contig[end]))
    

}

junc_join = function(c_parts, ord_junc_match, ord_seg_match) {
    ## browser()

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
}




            
        
write.fasta = function(grl, genome, path, append = T, compress = F, overwrite = FALSE) {
    if (file.exists(path) & !overwrite) {
        stop("file already exists... \n remove file first manually or set overwrite = TRUE")
    } else if (file.exists(path) & overwrite) {
        message("overwrite set to TRUE... removing file designated by path... giving you a moment...")
        Sys.sleep(7)
        cmd = paste0("rm ", path)
        system(cmd)
    }
    message("Total number of contigs: ", length(grl))
    message("Total number of bases: ", sum(as.numeric(width(grl.unlist(grl)))))
    for (i in 1:length(grl)) {
        message("Writing grl[[i]] where i=", i, " to FASTA")
        biostr = DNAStringSet(do.call(xscat, gr2seq(grl[[i]], genome)))
        names(biostr) = names(grl[i])
        message("Length of contig: ", width(biostr))
        ## writeXStringSet(biostr, path, append = append, compress = compress, format = 'fasta', width = 60)
        message("Chromosome/Contig ", seqnames(grl[[i]]), " with Allele ", grl[[i]]$allele, ", named ", names(grl[i]), " written to FASTA.")
    }
    message("Finished!")
    message("FASTA written to: ", path)
}

    
## NA12878 = read_vcf("~/DB/Platinum_Genomes/NA12878.vcf.gz")
## NA12877 = read_vcf("~/DB/Platinum_Genomes/NA12877.vcf.gz")

## NA12877 = renameSeqlevels(NA12877, sub("^chr", "", seqlevels(NA12877)))
## NA12878 = renameSeqlevels(NA12878, sub("^chr", "", seqlevels(NA12878)))

## alt_hi = CharacterList(NA12878$ALT)
## alt_lo = CharacterList(NA12877$ALT)

## ## hi_var = DNAStringSetList(mclapply(NA12878$ALT, function(x) x[1], mc.cores = 10)) ## Slow as HELL


## ##trick to get SINGLE haplotype for each "parent" faster
## test_altsplit = unstrsplit(alt_hi, sep = ",")
## this = test_altsplit[1430:1432]
## ## for 1st haplotype
## sub("^([A-Z]*)(\\,)([A-Z]*)", "\\1", this)
## ## for 2nd haplotype
## sub("^([A-Z]*)(\\,)([A-Z]*)", "\\3", this)


## alt_hi_split = unstrsplit(alt_hi, sep = ",")
## alt_lo_split = unstrsplit(alt_lo, sep = ",") 

## ## for 1st haplotype
## alt_hi_1 = sub("^([A-Z]*)(\\,)([A-Z]*)", "\\1", alt_hi_split)
## alt_lo_1 = sub("^([A-Z]*)(\\,)([A-Z]*)", "\\1", alt_lo_split)

## ## for 2nd haplotype
## alt_hi_2 = sub("^([A-Z]*)(\\,)([A-Z]*)", "\\3", alt_hi_split)
## alt_lo_2 = sub("^([A-Z]*)(\\,)([A-Z]*)", "\\3", alt_lo_split)



gr2seq = function(gr, genome) {

    NA12878 = read_vcf("~/DB/Platinum_Genomes/NA12878.vcf.gz")
    NA12877 = read_vcf("~/DB/Platinum_Genomes/NA12877.vcf.gz")

    NA12877 = renameSeqlevels(NA12877, sub("^chr", "", seqlevels(NA12877)))
    NA12878 = renameSeqlevels(NA12878, sub("^chr", "", seqlevels(NA12878)))

    alt_hi = CharacterList(NA12878$ALT)
    alt_lo = CharacterList(NA12877$ALT)

    alt_hi_split = unstrsplit(alt_hi, sep = ",")
    alt_lo_split = unstrsplit(alt_lo, sep = ",") 

    ## for 1st haplotype
    alt_hi_1 = sub("^([A-Z]*)(\\,)([A-Z]*)", "\\1", alt_hi_split)
    alt_lo_1 = sub("^([A-Z]*)(\\,)([A-Z]*)", "\\1", alt_lo_split)

    ## for 2nd haplotype
    alt_hi_2 = sub("^([A-Z]*)(\\,)([A-Z]*)", "\\3", alt_hi_split)
    alt_lo_2 = sub("^([A-Z]*)(\\,)([A-Z]*)", "\\3", alt_lo_split)


    hi_ix = which(gr$allele == "hi")
    lo_ix = which(gr$allele == "lo")
    hi_iix = findOverlaps(gr[gr$allele == "hi"], NA12878)
    lo_iix = findOverlaps(gr[gr$allele == "lo"], NA12877)
    hi_start = start(NA12878[to(hi_iix)]) - start(gr[hi_ix][from(hi_iix)]) +1
    hi_end = end(NA12878[to(hi_iix)]) - start(gr[hi_ix][from(hi_iix)]) +1
    hi_ranges = IRanges(hi_start, hi_end)
    hi_RList = split(hi_ranges, factor(hi_ix, 1:length(gr))[from(hi_iix)])
    hi_VList = CharacterList(split(alt_hi_1[to(hi_iix)], factor(hi_ix, 1:length(gr))[from(hi_iix)]))
    lo_start = start(NA12877[to(lo_iix)]) - start(gr[lo_ix][from(lo_iix)]) +1
    lo_end = end(NA12877[to(lo_iix)]) - start(gr[lo_ix][from(lo_iix)]) +1
    lo_ranges = IRanges(lo_start, lo_end)
    lo_RList = split(lo_ranges, factor(lo_ix, 1:length(gr))[from(lo_iix)])
    lo_VList = CharacterList(split(alt_lo_1[to(lo_iix)], factor(lo_ix, 1:length(gr))[from(lo_iix)]))

    ## SHIT -- need to trim the indel coordinates that hang off the breakpoint ends
    ## trim from beginning

    test_trim = function(RList, ends = NULL) {
        if (is.null(ends)) {
            which(unlist(mapply(function(x) {
                any(start(x) < 1) }, RList, SIMPLIFY = FALSE)))
        } else {
            which(unlist(mapply(function(x, y) {
                any(end(x) > y) }, RList, ends, SIMPLIFY = FALSE)))
        } }
        
        

    trim_these_hi_starts = test_trim(hi_RList)
    trim_these_hi_ends = test_trim(hi_RList, width(gr))
    trim_these_lo_starts = test_trim(lo_RList)
    trim_these_lo_ends = test_trim(lo_RList, width(gr))

    ## browser(expr = any(length(c(trim_these_hi_starts, trim_these_hi_ends, trim_these_lo_starts, trim_these_lo_ends)) > 0))


    trim_f = function(RList, VList, ends = NULL, trim_start = TRUE, trim_end = !trim_start) {
        if (trim_start) {
            new_start = mapply(function(RList_i, VList_i) {
                these_s = which(start(RList_i) < 1)
                repl = substr(VList_i[these_s], start = 1 + (1 - start(RList_i)[these_s]), stop = nchar(VList_i[these_s]))
                VList_i[these_s] = repl
                start(RList_i[these_s]) = pmax(start(RList_i[these_s]), 1)
                return(list(RList_i, VList_i))
                }, RList, VList, SIMPLIFY = F)
            output_list = list(IRangesList(lapply(new_start, function(x) x[[1]])), CharacterList(lapply(new_start, function(x) x[[2]])))
            return(output_list)
        } else if (trim_end) {
            new_end = mapply(function(RList_i, VList_i, ends_i) {                
                these_e = which(end(RList_i) > ends_i)
                repl = substr(VList_i[these_e], start = 1, stop = width(RList_i[these_e]) - (end(RList_i[these_e]) - ends_i))
                VList_i[these_e] = repl
                end(RList_i[these_e]) = pmin(end(RList_i[these_e]), ends_i)
                return(list(RList_i, VList_i))
                }, RList, VList, ends, SIMPLIFY = F)
            output_list = list(IRangesList(lapply(new_end, function(x) x[[1]])), CharacterList(lapply(new_end, function(x) x[[2]])))
            return(output_list)
        }
    }


    if (length(trim_these_hi_starts) > 0) {
        update_lst = trim_f(hi_RList[trim_these_hi_starts], hi_VList[trim_these_hi_starts])
        hi_RList[trim_these_hi_starts] = update_lst[[1]]
        hi_VList[trim_these_hi_starts] = update_lst[[2]]
    }
    if (length(trim_these_lo_starts) > 0) {
        update_lst = trim_f(lo_RList[trim_these_lo_starts], lo_VList[trim_these_lo_starts])
        lo_RList[trim_these_lo_starts] = update_lst[[1]]
        lo_VList[trim_these_lo_starts] = update_lst[[2]]
    }
    if (length(trim_these_hi_ends) > 0) {
        update_lst = trim_f(hi_RList[trim_these_hi_ends], hi_VList[trim_these_hi_ends], width(gr)[trim_these_hi_ends], trim_start = FALSE, trim_end = TRUE)
        hi_RList[trim_these_hi_ends] = update_lst[[1]]
        hi_VList[trim_these_hi_ends] = update_lst[[2]]
    }
    if (length(trim_these_lo_ends) > 0) {
        update_lst = trim_f(lo_RList[trim_these_lo_ends], lo_VList[trim_these_lo_ends], width(gr)[trim_these_lo_ends], trim_start = FALSE)
        lo_RList[trim_these_lo_ends] = update_lst[[1]]
        lo_VList[trim_these_lo_ends] = update_lst[[2]]
    }
    
    to_revComp = unlist(which(strand(gr)=="-"))
    this_seq = getSeq(genome, gr)

    conc = function(List1, List2) {
        mapply(function(x, y) c(x, y), List1, List2, SIMPLIFY = FALSE)
    }

    all_RList = IRangesList(conc(hi_RList, lo_RList))
    all_VList = CharacterList(conc(hi_VList, lo_VList))

    ## var_seq = replaceAt(this_seq, hi_RList, hi_VList)
    ## var_seq = replaceAt(this_seq, lo_RList, lo_VList)

    var_seq = replaceAt(this_seq, all_RList, all_VList)
    ## may need to do another apply function with c() to concatenate all the rangeslists and variantlists to do replaceAt() in one step
    ## currently this works because hi vars are only applied to one set of segments and lo vars to others
    ## if the 2nd set of vars were applied to a segment after changes were made, the coordinate system changes for a given segment
    var_seq[to_revComp] = reverseComplement(var_seq[to_revComp])
    ## this_seq[to_revComp] = reverseComplement(this_seq[to_revComp])
    return(var_seq)
}



find_int_junc = function(int_contigs) {
    ind = which(diff(int_contigs) != 1)
    tmp_1 = data.table(c = ind)
    tmp_2 = as.vector(rbind(tmp_1[,c], copy(tmp_1)[, c := c+1][,c]))
    split(tmp_2, 1:2)
}

## search_jmatch = function(junc_match, contigs) {
##     lapply(1:length(junc_match), function(x) {
##         rearr = lapply(1:length(contigs), function(y) {
##             ## browser(expr = (x == 3 & y == 101))
##             browser(expr = (x == 11 & y == 78))
##             this = match(abs(contigs[[y]]), abs(junc_match[[x]]))
##             tmp1 = na.omit(diff(this))
##             idx_1 = (which(diff(attributes(tmp1)$na.action) != 1)) + 1
##             idx_2 = idx_1 + 1
##             tmp2 = contigs[[y]][c(idx_1, idx_2)]
##             tmp3 = ifelse(tmp1 %in% c(1, 0, -1),
##                    ifelse(diff(tmp2)==sum(junc_match[[x]]),
##                           y,
##                           logical(0)),
##                    logical(0))
##             ## tmp2 = ifelse(length(tmp2) == 2, y, logical(0))
##             return(tmp3)
                    
##         })
##         rearr = unlist(rearr)
##         return(na.omit(rearr))
##     })
## }



## a rewrite... search_j match is bad

isNested <- function(x) {
    if (class(x) != "list") {
        stop("Expecting 'x' to be a list")
    }
    out <- any(sapply(x, is.list))
    return(out)
}

intercalate = function(...) {
    args = list(...)
    if (isNested(args)) {
        args = unlist(args, recursive = F)
    }
    s_along = lapply(args, seq_along)
    ord = order(do.call(c, s_along))
    conc = do.call(c, args)
    return(conc[ord])
}

intercalate_lst = function(...) {
    args = list(...)
    s_along = lapply(args, seq_along)
    ord = order(do.call(c, s_along))
    conc = do.call(c, args)
    return(conc[ord])
}


search_jmatch = function(junc_match, contigs, return_contig_idx = TRUE) {
    suppressWarnings( {
        lapply(1:length(junc_match), function(x) {
            rearr = lapply(1:length(contigs), function(y) {
                matches = which(! is.na(match(abs(contigs[[y]]), abs(junc_match[[x]]))))
                m_idx = which(diff(matches) == 1)
                j_m_idx = lapply(m_idx, function(x) intercalate(x, x + 1))
                if (return_contig_idx) {
                    c_idx = unlist(lapply(j_m_idx, function(z) {
                        if (diff(contigs[[y]][matches[z]]) == sum(junc_match[[x]])) y
                        else integer(0)
                    }))
                    return(c_idx)
                }
            })
            rearr = unlist(rearr)
            return(na.omit(rearr))
        })
    })
}

find_match = function(junc_match, contigs) {
    matches = which(! is.na(match(abs(contigs), abs(junc_match))))
    m_idx = which(diff(matches) == 1)
    j_m_idx = unlapply(m_idx, function(x) intercalate(x, x + 1))
    return(matches[j_m_idx])
}

grab_all_contigs = function(sim_object, step) {
    final_simulation_output = sim_object[step]
    cont = final_simulation_output[[1]]$output_contigs
    ho = final_simulation_output[[1]]$output_hold_out
    ho = ho[sapply(ho$contigs_lst, function(x) length(x)>0)]
    ho[,loose_left := unlist(loose_left)]
    ho[,loose_right := unlist(loose_right)]

    ## has_rearr = unlapply(ho$contigs_lst, function(x) any(diff(x) != 1))
    ## both_loose = ho[,loose_left+loose_right == 2]

    agg = rbind(cont, ho)
}


make_final_contigs = function(sim_object, clean_up_loose_ends = TRUE) {
    ## final_simulation_output = sim_object[[length(sim_object)]]
    
    ## cont = final_simulation_output$output_contigs
    ## ho = final_simulation_output$output_hold_out
    ## ho = ho[sapply(ho$contigs_lst, function(x) length(x)>0)]
    ## ho[,loose_left := unlist(loose_left)]
    ## ho[,loose_right := unlist(loose_right)]

    ## has_rearr = unlapply(ho$contigs_lst, function(x) any(diff(x) != 1))
    ## both_loose = ho[,loose_left+loose_right == 2]

    ## agg = rbind(cont, ho)
    ## browser()





    if (clean_up_loose_ends) {
        junction_grl = GRangesList(lapply(sim_object, function(x) x$junction_grl[[1]]))
        tile = karyograph(junctions = junction_grl)$tile
        agg = loose_end_match(sim_object, tile)
    } else {
        agg = grab_all_contigs(sim_object, length(sim_object))
    }
    


    keep_lg = agg[,(loose_left + loose_right == 0) |  unlapply(agg$contigs_lst, function(x) any(diff(x) !=1))]

    ## keep_lg = agg[,unlapply(agg$contigs_lst, function(x) any(diff(x) !=1))]
    keep = agg[keep_lg,]

    c_names = lapply(keep$contig_nm, function(x) paste0(x, collapse = "_"))
    c_dt = data.table(cont_nms = unlist(c_names), contig_idx = 1:length(c_names))
    c_dt = c_dt[, list(contig_idx, deriv = paste0(cont_nms, "_der_", .I)), by = cont_nms]
    setkey(c_dt, contig_idx)

    ref_cont = sim_object[[1]]$input_contigs
    ref_list_f = ref_cont$contigs_lst
    ref_list_r = lapply(ref_list_f, function(x) rev(-x))

    tmp_tel_1 = unlist(lapply(ref_cont$contigs_lst, function(x) head(x,1)), use.names = F)
    tmp_tel_2 = unlist(lapply(ref_cont$contigs_lst, function(x) tail(x,1)), use.names = F)
    left_tel = c(tmp_tel_1, -tmp_tel_2)
    right_tel = c(-tmp_tel_1, tmp_tel_2)
    
    left_fill = keep[, list(lapply(contigs_lst, function(x) {
        if (loose_left) {
              tmp = head(x, 1)
              if (tmp < 0) {
                  s = ref_list_r
                  idx = unlapply(s, function(y) match(names(tmp), names(y)))
              } else {
                  s = ref_list_f
                  idx = unlapply(s, function(y) match(names(tmp), names(y)))
              }
              chr_idx = which(!is.na(idx))
              this_match = na.omit(idx)
              ## browser()
              return(s[[chr_idx]][1:(this_match-1)])}})),
        by = contig_idx]
    ## browser()

    right_fill = keep[, list(lapply(contigs_lst, function(x) {
        ## print(x)
        if (loose_right) {
              tmp = tail(x, 1)
              if (tmp < 0) {
                  s = ref_list_r
                  idx = unlapply(s, function(y) match(names(tmp), names(y)))
              } else {
                  s = ref_list_f
                  idx = unlapply(s, function(y) match(names(tmp), names(y)))
              }
              chr_idx = which(!is.na(idx))
              this_match = na.omit(idx)
              return(s[[chr_idx]][(this_match+1): length(s[[chr_idx]])])}})),
        by = contig_idx]

    final_cont = mapply(function(left,mid,right) {
        c(left, mid, right) },
        left_fill$V1, keep$contigs_lst, right_fill$V1)

    tl =  karyograph(junc)$tile
    seqlevels(tl) = seqlevels(tl)[orderSeqlevels(seqlevels(tl), TRUE)]
    tl = sort(tl)
    pos_tile = tl %Q% (strand == "+")
    new_tl = c(pos_tile, pos_tile)
    names(new_tl) = c(1:length(pos_tile), -(1:length(pos_tile)))
    strand(new_tl[(length(new_tl)/2+1):length(new_tl)]) = "-"
    tl = new_tl

    grl = GRangesList(lapply(final_cont, function(cont) {
        allele = ifelse(grepl("hi", names(cont)), "hi", "lo")
        gr_cont = tl[match(cont, names(tl))]
        gr_cont$allele = allele
        return(gr_cont)
    }))

    names(grl) = c_dt$deriv
    return(grl)
    
}


find_jx_instantiations = function(sim_object, contigs) {
    juncs = lapply(sim_object, '[[', 2)
    ind_of_contigs = search_jmatch(juncs, contigs)
    inst_ind = which(elementNROWS(ind_of_contigs) > 0)
    inst_junc = juncs[inst_ind]
}



get_true_tile_cn = function(grl) {
    true_cn = as(coverage(unlist(grl)), "GRanges")
    pos_true_cn = setNames(true_cn, as.character(1:length(true_cn))); strand(pos_true_cn) = "+"; 
    neg_true_cn = setNames(true_cn, as.character((length(true_cn)+1):(2*length(true_cn)))); strand(neg_true_cn) = "-"
    true_tile = c(pos_true_cn, neg_true_cn)
}


## test with last_chance contigs and g_junc (incorporated snowman junctions)
true_edges = function(final_contigs, grl_junctions, check_sort = TRUE) {
    ## browser()

    if (check_sort) {
        final_contigs = sortSeqlevels(final_contigs)
        grl_junctions = sort(sortSeqlevels(grl_junctions))
    }
    
    true_cn = as(coverage(unlist(final_contigs)), "GRanges")
    pos_true_cn = setNames(true_cn, as.character(1:length(true_cn))); strand(pos_true_cn) = "+"; 
    neg_true_cn = setNames(true_cn, as.character((length(true_cn)+1):(2*length(true_cn)))); strand(neg_true_cn) = "-"
    true_tile = c(pos_true_cn, neg_true_cn)

    if (any(elementNROWS(grl_junctions) != 2)) {
        stop("All elements specified in grl_junctions must be of length 2")
    }
    
    if (all(unlist(width(grl_junctions)) == 1)) {
        message("All junctions within GRangesList are width 1,\nConverting to width 2")
        grl_junctions = endoapply(grl_junctions, function(x) {end(x) = end(x) + 1; x })
    } else if (all(unlist(width(grl_junctions)) == 2)) {
        message("All junctions within GRangesList are width 2,\nContinuing to copy profile")
    } else {
        stop("Input junctions within GRangesList are not of consistent width,\nWidth must be 1 or 2 for valid junctions")
    }

    junc_unlisted = grl.unlist(grl_junctions)

    junc_unlisted_dt = gr2dt(junc_unlisted)

    cn_olaps = findOverlaps(true_tile, junc_unlisted)
    dt_cn_olaps = as.data.table(cn_olaps)
    dt_cn_olaps[, ":="(grl.ix = junc_unlisted_dt[subjectHits, grl.ix], grl.iix = junc_unlisted_dt[subjectHits, grl.iix]) ]
    dt_cn_olaps = copy(dt_cn_olaps)[,junc_copy := abs(diff(true_tile[queryHits]$score)), keyby = .(grl.ix, subjectHits)]
    junc_copy = dt_cn_olaps[,abs(diff(true_tile[queryHits]$score)), keyby = .(grl.ix, subjectHits)][, unique(V1), by = grl.ix][, list(grl.ix, cn = V1)]

    cn_olaps2 = findOverlaps(gr.stripstrand(true_tile), gr.stripstrand(junc_unlisted))
    dt_cn_olaps2 = as.data.table(cn_olaps2)
    dt_cn_olaps2[, ":="(grl.ix = junc_unlisted_dt[subjectHits, grl.ix], grl.iix = junc_unlisted_dt[subjectHits, grl.iix]) ]
    dt_cn_olaps2[, seg_strand := as.character(strand(true_tile)[queryHits])]
    dt_cn_olaps2 = copy(dt_cn_olaps2)[,junc_copy := abs(diff(true_tile[queryHits]$score)), keyby = .(grl.ix, subjectHits, seg_strand)]

    junc_unlisted = grl.unlist(grl_junctions)
    junc_unlisted$ALT = unlist(junc_unlisted$ALT)
    start(junc_unlisted[strand(junc_unlisted) == "+"]) = start(junc_unlisted[strand(junc_unlisted) == "+"]) + 1
    end(junc_unlisted[strand(junc_unlisted) == "-"]) = end(junc_unlisted[strand(junc_unlisted) == "-"]) - 1
    olaps = findOverlaps(gr.stripstrand(true_tile), gr.stripstrand(junc_unlisted))
    dt_olaps = as.data.table(olaps)[order(subjectHits)]

    dt_olaps[, ":="(grl.ix = junc_unlisted_dt[subjectHits, grl.ix], grl.iix = junc_unlisted_dt[subjectHits, grl.iix]) ]

    dt_olaps[, seg_strand := as.character(strand(true_tile[queryHits]))]


    dt_olaps[, junc_strand := junc_unlisted_dt[subjectHits]$strand]
    dt_olaps[, tl_nm := names(true_tile)[queryHits]]

    junc_copy = dt_cn_olaps[,abs(diff(true_tile[queryHits]$score)), keyby = .(grl.ix, subjectHits)][, unique(V1), by = grl.ix][, list(grl.ix, cn = V1)]

    dt_olaps[, seg_lab := c("A", "B", "C", "D"), by = grl.ix]


    edge_dt = dt_olaps[, if (junc_strand[grl.iix == 1][1] == "+" & junc_strand[grl.iix == 2][1] == "+")
                             list(from = c(queryHits[seg_lab == "D"],
                                           queryHits[seg_lab == "B"]),
                                  to = c(queryHits[seg_lab == "A"],
                                         queryHits[seg_lab == "C"]))
                         else if (junc_strand[grl.iix == 1][1] == "-" & junc_strand[grl.iix == 2][1] == "-")
                             list(from = c(queryHits[seg_lab == "A"],
                                           queryHits[seg_lab == "C"]),
                                  to = c(queryHits[seg_lab == "D"],
                                         queryHits[seg_lab == "B"]))
                         else if (junc_strand[grl.iix == 1][1] == "+" & junc_strand[grl.iix == 2][1] == "-")
                             list(from = c(queryHits[seg_lab == "C"],
                                           queryHits[seg_lab == "B"]),
                                  to = c(queryHits[seg_lab == "A"],
                                         queryHits[seg_lab == "D"]))
                         else if (junc_strand[grl.iix == 1][1] == "-" & junc_strand[grl.iix == 2][1] == "+")
                             list(from = c(queryHits[seg_lab == "A"],
                                           queryHits[seg_lab == "D"]),
                                  to = c(queryHits[seg_lab == "C"],
                                         queryHits[seg_lab == "B"]))
                     , by = .(grl.ix) ]

    edge_copy = merge(edge_dt, junc_copy, by = "grl.ix", all.x = T, all.y = T)
    from_cn = edge_copy[,list(cn), keyby = .(key = as.character(from))]
    to_cn = edge_copy[,list(cn), keyby = .(key = as.character(to))]

    ref_es = true_tile

    ref_es_lst = split(ref_es, strand(ref_es), drop = T)
    pos_ref = ref_es_lst[[1]]
    neg_ref = ref_es_lst[[2]]

    pos_ref = split(pos_ref, seqnames(pos_ref))
    neg_ref = split(neg_ref, seqnames(neg_ref))

    ref_edges = rbindlist(lapply(c(pos_ref, neg_ref), function(tile) {
        ## browser()
        if (length(tile) > 1) {
            tl = names(tile)
            it = intercalate(tl[-length(tl)], tl[-1])
            ## it = it[-length(it)]
            dt_int = data.table(from = it[seq(1, length(it)-1, 2)],
                                to = it[seq(2, length(it), 2)], stringsAsFactors = FALSE)
            dt_int[, total_edges := max(tile[from]$score, tile[to]$score), by = .(from, to)]
            dt_int[, ref_edges_1st := as.numeric(total_edges - ifelse(!is.na(from_cn[as.character(from)]$cn), from_cn[as.character(from)]$cn, 0)), by = from]
            dt_int[, ref_edges_2nd := as.numeric(ref_edges_1st - ifelse(!is.na(to_cn[as.character(to)]$cn), to_cn[as.character(to)]$cn, 0)), by = to]
        } else {
            ## browser()
            return(NULL)
        }
    }))

    es = ref_edges[, list(from, to, cn = ref_edges_2nd, type = "reference")]
    es = rbindlist(list(es, edge_copy[, list(from, to, cn, type = "aberrant")]))
    es[, from := as.integer(from)]
    es[, to := as.integer(to)]
    es[, cn := as.integer(cn)]

    es = es[order(to)]


    es[, ":="(lwd = ifelse(type=="aberrant", log2(0.2*cn+2)+1, 1),
              lty = ifelse(type=='loose', 3, 1),
              col = ifelse(type=="aberrant",
                    ifelse(cn>0,
                           alpha("red", 0.4),
                           alpha("purple", 0.3)),
                    ifelse(type=="loose",
                           alpha("blue",0.6),
                           alpha("grey",0.2))),
              cex.arrow = 0,
              not.flat = type=="aberrant",
              v = ifelse(type=="aberrant", 2, 1),
              h = ifelse(type=="aberrant", 2, 1),
              dangle.w = 0.5)]
    es = es[!is.na(from) & !is.na(to)]

    return(es)
}


search_epoch = function(contig_idx, iteration_list, allele = NULL) {
    contigs_pools = lapply(iteration_list, function(lst) {
        ## browser()
        contigs_tbl = lst[[3]] ## input contigs table
        ho_tbl = lst[[6]] ## input hold out table
        contigs_lst1 = unlist(contigs_tbl[[1]]) ## contigs list
        contigs_lst2 = unlist(ho_tbl[[1]]) ## contigs list
        c(contigs_lst1, contigs_lst2)

    })
    ## browser()
    if (!is.null(allele)) {
        contig_idx = paste0(allele, "_", contig_idx)
        lg = lapply(contigs_pools, function(ctg) {
            contig_idx %in% names(ctg)
        })
        return(unlist(lg))
    } else {
        lg = lapply(contigs_pools, function(ctg) {
            abs(contig_idx) %in% abs(ctg)
        })
        return(unlist(lg))
    }
}


plot_contigs = function(single_iteration_list, tile, filename = NULL, windows = NULL, remove_unrearranged = FALSE, remove_double_loose = FALSE, keep_unloose = TRUE, id_contigs = FALSE, format = c("pdf", "png")) {
    
    these_contigs = rbind(single_iteration_list$output_contigs, single_iteration_list$output_hold_out)
    these_contigs = these_contigs[elementNROWS(these_contigs$contigs_lst) > 0]

    if (keep_unloose & remove_unrearranged) {
        lg = these_contigs[, unlist(loose_left) + unlist(loose_right) == 0]
        tmp_contigs = these_contigs[lg,]
    }
    
    if (remove_unrearranged) {
        lg  = these_contigs[,unlapply(contigs_lst, function(x) any(diff(x) !=1))]
        these_contigs = these_contigs[lg,]
        if (keep_unloose) {
            these_contigs = rbind(these_contigs, tmp_contigs)
        }
    } else if (remove_double_loose) {
        lg = these_contigs[, unlist(loose_left) + unlist(loose_right) == 2]
        these_contigs = these_contigs[!lg,]
    }
    
    grl = GRangesList(lapply(these_contigs$contigs_lst, function(cont) {
        allele = ifelse(grepl("hi", names(cont)), "hi", "lo")
        gr_cont = tile[match(cont, names(tile))]
        gr_cont$allele = allele
        return(gr_cont)
    }))



    loose_sink = (GRanges(seqnames = "1", IRanges(-1,-1)))
    values(loose_sink) = data.frame(tile.id = NA,is.tel = FALSE,ab.source = NA,ab.target = NA,allele = "lo")
    loose_sink = GRangesList(loose_sink)
    grl2 = c(grl, loose_sink)
    ## grl2 = Filter(function(x) length(x) > 0, grl2)

    


    col_contigs = endoapply(grl2, function(gr) {
        values(gr)[,"col"] = ifelse(values(gr)[, "allele"] == "hi",
                                    alpha("red", 0.1),
                                    alpha("blue", 0.1))
        return(gr)
    })

    ## select circular contigs to add bp to for circular paths
    ## here

    index_of_loose = nrow(these_contigs) + 1
    ## left_loose = these_contigs[, which(unlist(loose_left))]
    left_loose = unlapply(these_contigs$loose_left, function(x) if (! isTRUE(x) | length(x) == 0) FALSE else TRUE)
    ## right_loose = these_contigs[, which(unlist(loose_right))]
    right_loose = unlapply(these_contigs$loose_right, function(x) if (! isTRUE(x) | length(x) == 0) FALSE else TRUE)
    left_dt = data.table(left_loose = left_loose, to_nowhere = FALSE)
    right_dt = data.table(right_loose = right_loose, to_nowhere = FALSE)

    if (is.null(windows)) {
        windows = si2gr(seqinfo(col_contigs))
        n_paths = pmax(elementNROWS(col_contigs) - 1, 0)
        sum_paths = sum(n_paths)
        path_id_to_change = find_match(single_iteration_list$junc_match, these_contigs$contigs_lst[[1]])[1]
        path.col = rep("black", times = sum_paths)
        path.col[path_id_to_change] = "red"
        path.cex.arrow = rep(1, times = sum_paths)
        path.cex.arrow[path_id_to_change] = 3

    } else {
        if (inherits(windows, "character")) {
            windows = gstring(windows)
        }
        if (!inherits(windows, "GRanges")) {
            stop("windows arg does not specify proper coordinates")
        }
        
        olaps = lapply(col_contigs, function(gr) {
            lg = gr %over% windows
        })

        
        is_right_seg_in_window = unlapply(col_contigs, function(gr) if (length(gr) > 0) gr[length(gr)] %over% windows else FALSE)
        is_right_seg_in_window = is_right_seg_in_window[-length(is_right_seg_in_window)]
        is_left_seg_in_window = unlapply(col_contigs, function(gr) if (length(gr) > 0) gr[1] %over% windows else FALSE)
        is_left_seg_in_window = is_left_seg_in_window[-length(is_left_seg_in_window)]
        left_dt[!is_left_seg_in_window, c("left_loose", "to_nowhere") := list(FALSE, TRUE)]
        right_dt[!is_right_seg_in_window, c("right_loose", "to_nowhere") := list(FALSE,TRUE)]
        
        col_contigs = GRangesList(Map(function(gr, lg) {
            gr[lg]
        }, col_contigs, olaps))
        
        n_paths = pmax(elementNROWS(col_contigs) - 1, 0)
        sum_paths = sum(n_paths)
        these_matches_orig_coordinates = find_match(single_iteration_list$junc_match, these_contigs$contigs_lst[[1]])
        these_matches_windowed_contigs = find_match(single_iteration_list$junc_match, lapply(col_contigs, function(x) as.integer(names(x)))[[1]])
        is_junction_in_window = all(olaps[[1]][these_matches_orig_coordinates])
        if (n_paths[1] > 0 & is_junction_in_window) {
            path_id_to_change = these_matches_windowed_contigs[1]
            names(olaps[[1]]) = 1:length(olaps[[1]])
            tmp = olaps[[1]][olaps[[1]]]
            change_seg_color =  na.omit((which(tmp)[match(these_matches_orig_coordinates, names(tmp))]))
            tmp_vals = values(col_contigs[[1]])
            tmp_vals$col[change_seg_color] = ifelse(as.character(tmp_vals$allele[change_seg_color]) == "hi", alpha("red", 0.9), alpha("blue", 0.9))
            values(col_contigs[[1]]) = tmp_vals
        } else {
            path_id_to_change = integer(0)
        }
        path.col = rep("black", times = sum_paths)
        path.col[path_id_to_change] = "red"
        path.cex.arrow = rep(1, times = sum_paths)
        path.cex.arrow[path_id_to_change] = 3
        ## these_contigs_should_have_mark_path = which(unlapply(olaps, function(o) if (length(o) >0) length(which(diff(o) != 0)) > 1 else FALSE))

        mark_path_dt = rbindlist(lapply(1:length(n_paths), function(i) if (n_paths[i] > 0) {
                                                                            return(data.table(n = 1:n_paths[i], i))
                                                                        } else {
                                                                            return(NULL)
                                                                        }))
        mark_path_dt = mark_path_dt[order(as.character(i))] ## OH MY GOD... I COULDN'T INDEX THE RIGHT PATH BEFORE BECAUSE THE CONTIG INDEXES ARE SORTED AS CHARACTERS AND NOT INTEGERS (I.E. 1, 10, 11, 12,... NOT 1,2,3,4) this is the fix for now
        mark_path_dt[, path.ix := 1:.N]
        ## mark_path_dt = mark_path_dt[i %in% these_contigs_should_have_mark_path, list(n, i, path.ix, mark_color = diff(which(olaps[[i]])) != 1), by = i]
        mark_path_dt = mark_path_dt[, list(n, path.ix, mark_color = diff(which(olaps[[i]])) != 1), by = i]
        
        path.col[mark_path_dt[(mark_color), path.ix]] = "#A020F0CC" ## these are to mark junctions that have an intervening segment whose contig contains segments outside of the window but are flanked by those inside the window (i.e. the junction should go outside of the window to another contig, but another edge should come back inside)
        path.cex.arrow[mark_path_dt[(mark_color), path.ix]] = 5
    }

    left_loose_df = errnot(data.frame(from = index_of_loose, to = which(left_dt[, left_loose]), cn = 1, type = "loose", col    = "#00000080", h = 10, lwd = 5, lty = 1, cex.arrow = 0, v = 1, not.flat = TRUE, dangle.w = 1))
    right_loose_df = errnot(data.frame(from = which(right_dt[, right_loose]), to = index_of_loose, cn = 1, type = "loose", col    = "#00000080", h = 10, lwd = 5, lty = 1, cex.arrow = 0, v = 1, not.flat = TRUE, dangle.w = 1))
    left_nowhere_df = errnot(data.frame(from = index_of_loose, to = which(left_dt[, to_nowhere]), cn = 1, type = "loose", col    = "#00FF00CC", h = 10, lwd = 2, lty = 1, cex.arrow = 1, v = 1, not.flat = TRUE, dangle.w = 1))
    right_nowhere_df = errnot(data.frame(from = which(right_dt[, to_nowhere]), to = index_of_loose, cn = 1, type = "loose", col    = "#00FF00CC", h = 10, lwd = 2, lty = 1, cex.arrow = 1, v = 1, not.flat = TRUE, dangle.w = 1))
    es = rbind(left_loose_df, right_loose_df, left_nowhere_df, right_nowhere_df)

    if (id_contigs) {
        col_contigs = setNames(col_contigs, as.character(1:length(col_contigs)))
    }

    if (identical(format, c("pdf", "png")) |
        identical(format, c("png", "pdf"))) {
            message("selecting pdf format by default")
            format = "pdf"
    }
    
    if (is.null(filename)) {
        filename = "plot"
        if (format == "pdf") {
            filename = paste0(filename, ".pdf")
        } else if (format == "png") {
            filename  = paste0(filename, ".png")
        }
    }
    
    if (format == "pdf") {
        ppdf(plot(gTrack(col_contigs, edges = es, draw.paths = TRUE, name = "Karyotype", lwd.border = 1.2), links = single_iteration_list$junction_grl, windows = windows, path.col = path.col, path.cex.arrow = path.cex.arrow), filename = filename, cex = 1)
    } else if (format == "png") {
        ppng(plot(gTrack(col_contigs, edges = es, draw.paths = TRUE, name = "Karyotype", lwd.border = 1.2), links = single_iteration_list$junction_grl, windows = windows, path.col = path.col, path.cex.arrow = path.cex.arrow), filename = filename, cex = 1)
    }
}

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
}

                                
contigs2dt = function(lst, tile) {
    tmp_dt = list2dt(lst)
    tmp_gr = tile[match(as.character(tmp_dt$unlisted),names(tile))]
    tmp_dt = cbind(tmp_dt, gr2dt(tmp_gr)[, list(seqnames, start, end, strand, width)])
    return(tmp_dt)
}

cdt2grl = function(lst, tile) {
    tmp_gr = dt2gr(contigs2dt(lst, tile)[, width0 := width][, width := NULL])
    names(tmp_gr) = tmp_gr$unlisted
    return(setNames(split(tmp_gr, tmp_gr$.ix), NULL))
}

na2false = function(v) {
    ## v = ifelse(is.na(v), v, FALSE)
    v[is.na(v)] = FALSE
    as.logical(v)
}

null2false = function(v) {
    !length(v) == 0
}


gv = function(gr, field) {
    values(gr)[, field]
}

setAllNames = function(v, nm) {
    allnames = rep(nm, length.out = length(v))
    setNames(v, allnames)
}



plot_all_contigs = function(iteration_list, directory = NULL, windows = NULL, mc.cores = 1, id_contigs = FALSE, format = "pdf") {
    grabbed_junctions = Reduce(c, (lapply(iteration_list, "[[", 10)))
    tile = sorted_tile(grabbed_junctions)
    if (! is.null(directory)) {
        directory = paste0(directory, "/")
        system(paste('mkdir -p', directory))
    }
    width = ceiling(log10(abs(length(new_sim))))
    mclapply(1:length(iteration_list), function(i) {
        filename = paste0(directory, formatC(i, width = width, flag = 0), "_", "junc_", iteration_list[[i]]$junc_idx, ".", format)
        suppressWarnings(plot_contigs(iteration_list[[i]], tile = tile, filename = filename, windows = windows, id_contigs = id_contigs, format = format))
        message("plotted: ", filename)
    }, mc.cores = mc.cores)
}



sorted_tile = function(junc) {
    kag = karyograph(junc)
    tile =  karyograph(junc)$tile
    seqlevels(tile) = seqlevels(tile)[orderSeqlevels(seqlevels(tile), TRUE)]
    tile = sort(tile)
    pos_tile = tile %Q% (strand == "+")
    new_tile = c(pos_tile, pos_tile)
    names(new_tile) = c(1:length(pos_tile), -(1:length(pos_tile)))
    strand(new_tile[(length(new_tile)/2+1):length(new_tile)]) = "-"
    tile = new_tile
}


    
## system("mkdir -p ~/public_html/Simulation/sim_run_seed_1987")
## paste0("Simulation/sim_run_seed_1987/", formatC(i, width = 3, flag = 0), "_", "junc_", rand_id, ".pdf")) ## the walk must be grangeslist object))





































matching_algo = function(gr_match,
                         junction_grl,
                         junc_match,
                         which_breakpoint_ix, ## either 1 or 2
                         padding = 3e4) {


    
    na2false = function(v) {
        ## v = ifelse(is.na(v), v, FALSE)
        v[is.na(v)] = FALSE
        as.logical(v)
    }


    gv = function(gr, field) {
        values(gr)[, field]
    }
 
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
        
        possible_matches = possible_matches %Q% (na2false(loose_left_reference_face == "right") |
                                                 na2false(loose_right_reference_face == "right") |
                                                 !is.na(eval((parse(text = junc_match_field)))))
    } else if (junc_match[this_id] > 0) {
        possible_matches = possible_matches %Q% (na2false(loose_left_reference_face == "left") |
                                                 na2false(loose_right_reference_face == "left") |
                                                 !is.na(eval((parse(text = junc_match_field)))))
    }
    

    ## values(possible_matches)[, junc_match_field] = na.omit(gv(possible_matches,junc_match_field ))
    values(possible_matches)[, junc_match_field] = this_junc_match

    ## junction contig face also follows same sign convention as if the junction actually matches to the contig
    values(possible_matches)[, junc_contig_face_field] = ifelse(sign(gv(possible_matches, junc_match_field)) == sign(possible_matches$seg_id), "right", "left")
    

    ## now how to select?
    ## if there is loose end, always grab that one?
    ## or do selection based on # segs and/or contig distance?

    ## collecting stats for possible matches

    contig_side_to_fuse = ifelse(sign(gv(possible_matches, "seg_id")) == sign(gv(possible_matches, junc_match_field)),
                                 "right",
                                 "left")
    
    values(possible_matches)[, "contig_side_to_fuse"] = contig_side_to_fuse

    num_exact_matches = length(which(abs(possible_matches$seg_id) == abs(gv(possible_matches, junc_match_field))))

    ## will_fuse_loose = !is.na(possible_matches$loose_right_reference_face) | !is.na(possible_matches$loose_left_reference_face)
    if (this_junc_match < 0) {
        will_fuse_loose = na2false(possible_matches$loose_right_reference_face == "right" & this_junc_match < 0 | possible_matches$loose_left_reference_face == "right" & this_junc_match < 0)
    } else if (this_junc_match > 0) {
        will_fuse_loose = na2false(possible_matches$loose_right_reference_face == "left" & this_junc_match > 0 | possible_matches$loose_left_reference_face == "left" & this_junc_match > 0)
    }
    
    
    possible_matches$will_fuse_loose = will_fuse_loose
    number_of_segments_to_be_fused = ifelse(gv(possible_matches, junc_contig_face_field) == "left",
                                            possible_matches$num_from_left,
                                            possible_matches$num_from_right)
    values(possible_matches)[,"number_of_segments_to_be_fused"] = number_of_segments_to_be_fused

    contig_size_to_be_fused = ifelse(gv(possible_matches, junc_contig_face_field) == "left",
                                     possible_matches$cumdist_from_left,
                                     possible_matches$cumdist_from_right)
    possible_matches$contig_size_to_be_fused = contig_size_to_be_fused
    

    ## below enforces fusion of loose end always: using loose end that fuses the most number of segments

    ## 
    if (any(will_fuse_loose)) {        
        if (length(which(will_fuse_loose)) > 1) {
            ## browser()
            ## max_fused = max(possible_matches[possible_matches$will_fuse_loose]$number_of_segments_to_be_fused)
            max_fused = max(possible_matches[possible_matches$will_fuse_loose]$contig_size_to_be_fused) ## changing heuristic to maximum length and also not excluding non loose fusion
            ## tmp = possible_matches[possible_matches$number_of_segments_to_be_fused == max_fused]
            tmp = possible_matches[possible_matches$will_fuse_loose & possible_matches$contig_size_to_be_fused == max_fused]
            tmp = sample(tmp[tmp$contig_size_to_be_fused == max(tmp$contig_size_to_be_fused)], 1)
            final_match = tmp
            sign_seg = sign(gv(final_match, "seg_id"))
        } else {
            ## max_fused = max(possible_matches$contig_size_to_be_fused) ## added this in
            tmp = possible_matches[possible_matches$will_fuse_loose]
            ## tmp = possible_matches[possible_matches$contig_size_to_be_fused == max_fused] ## added this in
            ## tmp = sample(tmp[tmp$contig_size_to_be_fused == max(tmp$contig_size_to_be_fused)], 1) ## added this in
            final_match = tmp
            sign_seg = sign(gv(final_match, "seg_id"))
        }
        
        if (junction_reference_face == "left") {
        ## if (gv(final_match, junc_contig_face_field) == "right") {
            if (abs(gv(final_match, "seg_id")) < abs(gv(final_match, junc_match_field))) {
                fill_in_segs = seq(abs(gv(final_match, "seg_id")) + 1, abs(gv( final_match, junc_match_field))) * sign_seg
                segs_to_trim = NULL
                where_to_fill = gv(final_match, ".iix")
                where_to_trim = NULL
            } else if (abs(gv(final_match, "seg_id")) > abs(gv(final_match, junc_match_field))) {
                fill_in_segs = NULL
                segs_to_trim = seq(abs(gv( final_match, junc_match_field)) + 1, abs(gv(final_match, "seg_id"))) * sign_seg
                where_to_fill = NULL
                segment_diff = abs(gv(final_match, "seg_id")) - abs(gv(final_match, junc_match_field))
                where_to_trim = (gv(final_match, ".iix")-segment_diff+1):gv(final_match, ".iix")
            } else {
                fill_in_segs = NULL
                segs_to_trim = NULL
                where_to_fill = NULL
                where_to_trim = NULL
            }
        }
        else if (junction_reference_face == "right") {
        ## else if ((gv(final_match, junc_contig_face_field)) == "left") {
            if (abs(gv(final_match, "seg_id")) < abs(gv(final_match, junc_match_field))) {
                segs_to_trim = seq(abs(gv(final_match, "seg_id")), abs(gv( final_match, junc_match_field)) - 1) * sign_seg
                fill_in_segs = NULL
                segment_diff = abs(gv(final_match, junc_match_field)) - abs(gv(final_match, "seg_id"))
                where_to_trim = 1:segment_diff
                where_to_fill = NULL
            } else if (abs(gv(final_match, "seg_id")) > abs(gv(final_match, junc_match_field))){
                fill_in_segs = seq(abs(gv(final_match, junc_match_field)), abs(gv( final_match, "seg_id")) - 1) * sign_seg
                segs_to_trim = NULL
                where_to_fill = gv(final_match, ".iix")
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
    ## where_to_add = ifelse(sign(gv(final_match, "seg_id")) == sign(gv(final_match, junc_match_field)),
    ##                       gv(final_match, "num_from_right"),
    ##                       gv(final_match, "num_from_left"))
    where_to_add = gv(final_match, "num_from_left")
    values(final_match)[, "where_to_add"] = where_to_add


    return(list(final_match = final_match,
                fill_in_segs = fill_in_segs,
                segs_to_trim = segs_to_trim,
                where_to_fill = where_to_fill,
                where_to_trim = where_to_trim))

}





process_matching_algo = function(
    match_lst,
    contigs_dt,
    which_breakpoint_ix,
    aggregated_contigs_tbl) {

    gv = function(gr, field) {
        values(gr)[, field]
    }

    
    final_match = match_lst[['final_match']]
    which_pool = gv(final_match, "pool")
    matched_ix = gv(final_match, ".ix")
    junc_match_idx = gv(final_match, "junc_match_idx")
    contigs_dt = contigs_dt[order(.ix, .iix)] ## back to somatic ordering    
    matched = contigs_dt[.ix == matched_ix]
    junc_match_field = paste0("junc_match", which_breakpoint_ix)
    ## junc_match = na.omit(c(gv(final_match, "junc_match1"), gv(final_match, "junc_match2")))
    junc_match = gv(final_match, junc_match_field)
    where_to_add = gv(final_match, "where_to_add")
    select_from_agg = gv(final_match, "contig_idx")
    ## cont_match = grab_contig_match(match_lst)
    fill_in_segs = match_lst[["fill_in_segs"]]
    which_side = gv(final_match, "contig_side_to_fuse")
    will_fuse_loose = gv(final_match, "will_fuse_loose")
    contigs_len = gv(final_match, "num_segs")
    segs_to_trim = match_lst[["segs_to_trim"]]
    where_to_trim = match_lst[["where_to_trim"]]
    where_to_fill = match_lst[["where_to_fill"]]

    in_current = gv(final_match, "pool") == "current"
    in_hold_out = gv(final_match, "pool") == "hold_out"

    

    adjust_cont_match = ifelse(is.null(fill_in_segs), 0, length(fill_in_segs))

    if (will_fuse_loose) {
        ## cont_match = ifelse(which_side == "left", where_to_add + adjust_cont_match, 1)
        cont_match = ifelse(which_side == "left", contigs_len + adjust_cont_match, 1)
        contigs_len = contigs_len +  adjust_cont_match
    } else {
        ## cont_match = where_to_add
        cont_match = gv(final_match, "num_from_left")
    }

    contigs_dt = contigs_dt[order(.ix, .iix)] ## back to somatic ordering    
    matched = contigs_dt[.ix == matched_ix]

    allelic_nm = strsplit(matched[.iix == where_to_add, nm], "_")[[1]][1]

    if (!is.null(fill_in_segs)) fill_in_segs = setNames(fill_in_segs, paste(allelic_nm, sub("\\-", "", fill_in_segs), sep = "_"))

    matched_contig = original_contig = setNames(matched[, seg_id], matched[, nm])
    original_indices = 1:length(original_contig)
    

    ## if (gv(final_match, "contig_side_to_fuse") == "left") {
    if (junc_match > 0) {
        if (gv(final_match, "contig_side_to_fuse") == "left") {
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
        if (gv(final_match, "contig_side_to_fuse") == "left") {
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


    ra =  aggregated_contigs_tbl[contig_idx == select_from_agg & pool == which_pool]
    ra[, pool := NULL][, n := NULL][, c("in_contig",
                                        "junc_match_idx",
                                        "junc_match",
                                        "cont_match")
                                    := list(TRUE,junc_match_idx, junc_match, cont_match)]
    ra[, contigs_lst := list(list(matched_contig))]
    ra[, seg_match := contigs_lst[[1]][cont_match], by = list(contig_idx,junc_match)]
    ra[, contigs_len := contigs_len]
    ra[, will_fuse_loose := will_fuse_loose]
    return(list(single_ra  = ra, in_current = in_current, in_hold_out = in_hold_out, adjusted_indices = new_indices, add_to_first = add_to_first, add_to_last = add_to_last))
    
}

isNone = function(v) {
    length(v) == 0
}


special_match = function(v1, v2) {
    d1 = data.table(v1 = v1)
    d2 = data.table(v2 = v2)
    d1[, count := 1:.N, by = v1]
    d2[, count := 1:.N, by = v2]
    browser()
    paste(d1$v1, d1$count)
    return(NULL)
}




matching_algo2 = function(gr_match, ## matches contig within window
                         junction_grl,
                         junc_match,
                         padding = 3e4) {

    na2false = function(v) {
        ## v = ifelse(is.na(v), v, FALSE)
        v[is.na(v)] = FALSE
        as.logical(v)
    }


    gv = function(gr, field) {
        values(gr)[, field]
    }

    

    possible_matches1 = get_matches(gr_match, junction_grl, junc_match,  1, padding = padding)
    possible_matches2 = get_matches(gr_match, junction_grl, junc_match,  2, padding = padding)

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
    
    
    eligible1 = closest1 %&% boundary
    eligible2 = closest2 %&% boundary

    return(list(j1_matches = eligible1, j2_matches = eligible2))

}

select_match = function(possible_matches, which_junc, junc_match) {

    junction_reference_face = ifelse(junc_match[which_junc] < 0, "left", "right")

    possible_matches = possible_matches[[which_junc]]
    
    junc_match_field = paste0("junc_match", which_junc)
    will_fuse_loose = gv(possible_matches, "will_fuse_loose")
    ## below enforces fusion of loose end always: using loose end that fuses the most number of segments
    num_exact_matches = length(which(abs(possible_matches$seg_id) == abs(gv(possible_matches, junc_match_field))))
    
    ## 
    if (any(will_fuse_loose)) {        
        if (length(which(will_fuse_loose)) > 1) {
            ## browser()
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
            sign_seg = sign(gv(final_match, "seg_id"))
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
            sign_seg = sign(gv(final_match, "seg_id"))
        }
        
        if (junction_reference_face == "left") {
        ## if (gv(final_match, junc_contig_face_field) == "right") {
            if (abs(gv(final_match, "seg_id")) < abs(gv(final_match, junc_match_field))) {
                fill_in_segs = seq(abs(gv(final_match, "seg_id")) + 1, abs(gv( final_match, junc_match_field))) * sign_seg
                segs_to_trim = NULL
                where_to_fill = gv(final_match, ".iix")
                where_to_trim = NULL
            } else if (abs(gv(final_match, "seg_id")) > abs(gv(final_match, junc_match_field))) {
                fill_in_segs = NULL
                segs_to_trim = seq(abs(gv( final_match, junc_match_field)) + 1, abs(gv(final_match, "seg_id"))) * sign_seg
                where_to_fill = NULL
                segment_diff = abs(gv(final_match, "seg_id")) - abs(gv(final_match, junc_match_field))
                where_to_trim = (gv(final_match, ".iix")-segment_diff+1):gv(final_match, ".iix")
            } else {
                fill_in_segs = NULL
                segs_to_trim = NULL
                where_to_fill = NULL
                where_to_trim = NULL
            }
        }
        else if (junction_reference_face == "right") {
        ## else if ((gv(final_match, junc_contig_face_field)) == "left") {
            if (abs(gv(final_match, "seg_id")) < abs(gv(final_match, junc_match_field))) {
                segs_to_trim = seq(abs(gv(final_match, "seg_id")), abs(gv( final_match, junc_match_field)) - 1) * sign_seg
                fill_in_segs = NULL
                segment_diff = abs(gv(final_match, junc_match_field)) - abs(gv(final_match, "seg_id"))
                where_to_trim = 1:segment_diff
                where_to_fill = NULL
            } else if (abs(gv(final_match, "seg_id")) > abs(gv(final_match, junc_match_field))){
                fill_in_segs = seq(abs(gv(final_match, junc_match_field)), abs(gv( final_match, "seg_id")) - 1) * sign_seg
                segs_to_trim = NULL
                where_to_fill = gv(final_match, ".iix")
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
    ## where_to_add = ifelse(sign(gv(final_match, "seg_id")) == sign(gv(final_match, junc_match_field)),
    ##                       gv(final_match, "num_from_right"),
    ##                       gv(final_match, "num_from_left"))
    where_to_add = gv(final_match, "num_from_left")
    values(final_match)[, "where_to_add"] = where_to_add


    return(list(final_match = final_match,
                fill_in_segs = fill_in_segs,
                segs_to_trim = segs_to_trim,
                where_to_fill = where_to_fill,
                where_to_trim = where_to_trim))

}
