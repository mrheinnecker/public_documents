######################
## tertiary methods ##
######################
extract_promoters_from_gencode <- function(gencode){
  apply(gencode, 1, function(x){
    tibble(probeID=str_split(x[["relevant_cpgs"]], ";") %>% unlist(),
           gene=x[["gene"]]) %>%
      return()
  }) %>% bind_rows()%>%
    return()
}
store_table <- function(opt, var_name, path, file_name){
  if(isTRUE(!is.null(opt[[var_name]]))){
    if(nrow(opt[[var_name]])!=0){
      write_tsv(opt[[var_name]] %>%
                  mutate(ID=opt$id), 
                file=file.path(path, file_name))
    }
  }
}
prepare_germline_output_message <- function(opt){
  if(is.null(opt$df_germ_topart_relevant)){
    opt$germline_mes <- 'no TOP-ART relevant mutations found -> 0 pts'
  } else {
    opt$germline_mes <- opt$df_germ_unfiltered %>% 
      filter(gene %in% opt$ref_gen_germline_topart_relevant) %>% 
      mutate_if(is.numeric, round, digits=2) %>%
      select(mutated_gene=gene, ACMG_class, reliability=rel) %>%
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column('names') %>% 
      apply(.,1,correct_space_in_printed_tables) %>% 
      apply(.,1,function(x){paste(x,collapse='\t')}) 
  }  
  return(opt)
}
prepare_somatic_output_message <- function(opt){
  if(!is.null(opt$pre_scoring)){
    if(!nrow(opt$pre_scoring_topart_relevant)==0){
      pre_scoring_mes <- opt$pre_scoring_topart_relevant %>% 
        mutate_at(vars(c("tcn", "maf", "aff_cp", "wt_cp")),
                  funs(as.numeric)) %>%
        mutate_if(is.numeric, round, digits=2) %>%
        select(-ref, -alt, -chr, -pos, -purity) %>%
        mutate(class=class %>% 
                 str_replace_all('nonsynonymous ', 'ns-') %>% 
                 str_replace_all('frameshift', 'fs') %>%
                 str_replace_all('insertion', 'ins') %>%
                 str_replace_all('deletion', 'del') %>%
                 str_replace_all('germline', 'germ'),
        ) %>%
        t() %>% as.data.frame() %>% rownames_to_column('names') 
      pre_scoring_mes[is.na(pre_scoring_mes)] <- ''
      opt$pre_scoring_mes <- pre_scoring_mes %>%
        apply(.,1,correct_space_in_printed_tables) %>% 
        apply(.,1,function(x){paste(x,collapse='\t')}) %>%
        c('\n  Pre-evaluation of every variant:',.)
      if(is.null(opt$final_scoring_topart_relevant)){
        opt$final_scoring_mes <- ""
      } else {
        opt$final_scoring_mes <- 
          sapply(c(1:nrow(opt$final_scoring_topart_relevant)),function(n){
          if(!is_null(opt$phasing_info_topart_relevant[[n]])){
            pi <- opt$phasing_info_topart_relevant[[n]]
            if(nrow(pi)==1){
              pi <- pi %>% select(-info)
            }
            pi[is.na(pi)] <- ''
            c(paste(opt$final_scoring_topart_relevant[n,], collapse='\t'), 
              pi %>% select(-score) %>%
                mutate_if(is.numeric, round, digits=2) %>%
                t() %>% as.data.frame() %>% rownames_to_column('names') %>% 
                apply(.,1,correct_space_in_printed_tables) %>% 
                apply(.,1,function(x){paste(x,collapse='\t')}) %>% 
                paste('  ',.) %>% yellow()) %>%
              return()
          } else {
            paste(opt$final_scoring_topart_relevant[n,], collapse='\t') %>% 
              return()
          }
        }) %>% 
          unlist() %>% 
          c(paste0('\n  Final estimation of somatic score for each gene: ',
                   "\n    gene\tn_mut\tscore"),.)
      }
    } else {
      opt$pre_scoring_mes <- 'no TOP-ART relevant mutations found -> 0 pts'
    }
  } else {
    opt$pre_scoring_mes <- 'no TOP-ART relevant mutations found -> 0 pts'
  }  
  return(opt)
}
print_full_output_message <- function(opt){
  mes <- c(
    PID=paste0('PID:\t\t', opt$pid, '\n', 
               'seq-method:\t', opt$seq_method, '\n',
               'sample:\t\t', opt$sample, '\n',
               'purity:\t\t', opt$purity, '\n',
               'TOP-ART-score:  ', opt$TOP_ART_score),
    P1_fin = c(underline(paste('P1 (signature) criterion:', 
                               opt$signature_score)),
               opt$signature_scoring,
               opt$signature_add) %>% stri_remove_empty() %>%
      paste(collapse = '\n    '),
    P2_fin = c(underline(paste('P2 (rearrangement) criterion:' , 
                               opt$rearrangement_score)),
               opt$rearrangement_scoring
    ) %>% 
      paste(collapse = '\n    '),
    G2_fin = c(underline(paste('G2 (germline) criterion:', 
                               opt$germline_score)),
               opt$germline_mes
    ) %>% stri_remove_empty() %>%
      paste(collapse = '\n    '),
    G1_fin = c(underline(paste('G1 (somatic) criterion:', 
                               opt$somatic_score)),
               opt$pre_scoring_mes, 
               opt$final_scoring_mes
    ) %>% 
      paste(collapse = '\n    ')
  ) %>% paste(collapse = '\n \n') %>%
    str_replace_all('unclear', red('unclear')) %>%
    cat(., '\n')
}
# function to calculate maf (minor allele frequency) for snvs
calc_maf_snv <- function(BASE){
  # it works for vcf files and their INFO column... each columns contains a
  # string which itself has many information about the snv. 
  # all that information is separated by ';'
  # first split by ';'
  x <- str_split(BASE, ';') %>% unlist()
  # extract the part of the vector which has the 'DP4'
  erg <- x[grepl('DP4', x)] %>%
    # split by = and take the part after it
    str_split('=') %>%
    map_chr(.,2) %>%
    # split by ',' and take all entries
    str_split(',') %>% unlist() %>%
    as.numeric()
  # calculate the maf and return result
  return(sum(erg[3], erg[4])/sum(erg))
}
# function to calculate MAF from the TUMOR field of vcf files for indels
calc_maf_indel <- function(BASE){
  x <- str_split(BASE, ':') %>% 
    unlist() %>% 
    .[c(length(.)-1,length(.))] %>% 
    as.numeric()
  x[2]/x[1] %>%
    return()
}
aff_germ_copies <- function(maf_tumor, tcn, purity, chr, sex, Cc){
  if(sex=="male"&chr=="X"){
    (maf_tumor*(purity*tcn+1-purity)-Cc*(1-purity))/purity %>%
      return()
  } else {
    (maf_tumor*(purity*tcn+2*(1-purity))-Cc*(1-purity))/purity %>%
      return()
  }  
}
aff_som_copies <- function(maf, tcn, purity, chr, sex){
  if(sex=="male"&chr=="X"){
    maf*(tcn+1/purity-1) %>%
      return()
  } else {
    maf*(tcn+2/purity-2) %>%
      return()
  }
}
correct_space_in_printed_tables <- function(x){
  max_char <- sapply(x, nchar) %>% 
    max()
  filled <- sapply(x, function(y){
    if(nchar(y)<max_char){
      erg <- paste(rep(' ', max_char-nchar(y)),collapse = '') %>% 
        paste0(y,.)
    } else {
      erg <- y
    }
    return(erg)
  })
  return(filled)
}
prepare_raw_bam_file <- function(file, chr1, chr2, pos1, pos2){
  ref_pos1 <- pos1 %>% as.numeric()
  ref_pos2 <- pos2 %>% as.numeric()
  ref_chr1 <- chr1 
  ref_chr2 <- chr2 
  if(pos1>pos2){
    ref_pos1 <- pos2 %>% as.numeric()
    ref_pos2 <- pos1 %>% as.numeric()
    ref_chr1 <- chr2 
    ref_chr2 <- chr1 
  }
  c2 <- GenomicAlignments::readGAlignmentPairs(
    file,
    param=Rsamtools::ScanBamParam(
      which=GenomicRanges::GRanges(seqnames = ref_chr1, 
                                   ranges = ref_pos1),
      what = c('seq')
      )) %>% 
    as.data.frame()
  if(nrow(c2)==0){
    return(tibble())
  } else {
    c2 %>% mutate(
      pos1_in_first = apply(.,1, function(read){
        ifelse(ref_pos1 %in% c(read[['start.first']]:read[['end.first']]),
               T,F)
      }),
      pos1_in_last = apply(.,1, function(read){
        ifelse(ref_pos1 %in% c(read[['start.last']]:read[['end.last']]),
               T,F)
      }),
      pos2_in_first = apply(.,1, function(read){
        ifelse(ref_pos2 %in% c(read[['start.first']]:read[['end.first']]),
               T,F)
      }),
      pos2_in_last = apply(.,1, function(read){
        ifelse(ref_pos2 %in% c(read[['start.last']]:read[['end.last']]),
               T,F)
      })
    ) %>% mutate(
      keep=apply(.,1, function(x){
        ifelse(TRUE %in% c(x[['pos1_in_first']], 
                           x[['pos1_in_last']])&
                 (TRUE %in% c(x[['pos2_in_first']], 
                              x[['pos2_in_last']])), 
               T, F)})
    ) %>% 
      filter(keep==T) %>% 
      return()
  }
}
correct_read_order <- function(x){
  if(x[['start.last']]<x[['start.first']]){
    x1 <- x %>% as.data.frame() %>% select(ends_with('first'))
    n1 <- names(x1) %>% str_replace_all('.first', '.last')
    x2 <- x %>% as.data.frame() %>%  select(ends_with('last'))
    n2 <- names(x2) %>% str_replace_all('.last', '.first')
    y <- cbind(x2, x1) %>% set_names(nm=c(n2,n1)) 
  } else {
    y <- x
  }
  return(y)
}
verify_reads <- function(bam, comb_line){
  # function to verify reads which were filtered from bam file and contain the 
  # positions of both mutations  
  ############################################################################## 
  # at first we load all reference positions/mutated bases and reference bases  
  ref_pos1 = comb_line[['pos1']] %>% as.numeric()
  ref_pos2 = comb_line[['pos2']] %>% as.numeric()
  ref_alt1 = comb_line[['alt1']]
  ref_alt2 = comb_line[['alt2']]
  ref_ref1 = comb_line[['ref1']]
  ref_ref2 = comb_line[['ref2']]
  ##############################################################################
  # now a df is created containing all selected reads from the bam file as rows  
  df <- bam %>%
    mutate(
      result = apply(.,1, function(input_read){
        ########################################################################
        # only main chromosomes will be supported
        if(!input_read[['seqnames.first']] %in% c(1:23, 'X', 'Y')){
          stop <- F
          # if the first part is not in a valid chromosome the complete 
          # sequence and the info about it is taken from the second read-end         
          cigar_complete <- input_read[['cigar.last']]
          seq_complete <- input_read[['seq.last']]
          # and the data which is required for processing 
          # is stored in a new variable "read"         
          read <- tibble(x=input_read %>% names(), 
                         y=input_read  %>% as.character()) %>%
            column_to_rownames(var='x') %>% 
            t() %>% 
            as.data.frame() %>%
            select(ends_with('last')) %>%
            set_names(nm=names(.)%>%
                        str_replace('last', 'first'))
        }
        ######################################################################## 
        # ... for the second part of the paired end read applies the same        
        if(!input_read[['seqnames.last']] %in% c(1:23, 'X', 'Y')){
          stop <- F
          cigar_complete <- input_read[['cigar.first']]
          seq_complete <- input_read[['seq.first']]
          read <- tibble(x=input_read %>% names(), 
                         y=input_read  %>% as.character()) %>%
            column_to_rownames(var='x') %>% 
            t() %>% 
            as.data.frame() %>%
            select(ends_with('first'))
        }
        ########################################################################  
        # if both ends are located in valid chromosomes... 
        if(input_read[['seqnames.last']] %in% c(1:23, 'X', 'Y')&
           input_read[['seqnames.first']] %in% c(1:23, 'X', 'Y')){
          # ... we will paste the two paired reads together to one read, 
          # therefore the first read must be the one that starts first...
          # therefore a function was created:
          read_pre <- tibble(x=input_read %>% names(), 
                             y=input_read  %>% as.character()) %>%
            column_to_rownames(var='x') %>% t() %>% as.data.frame()
          read_pre$start.first <- as.numeric(read_pre$start.first)
          read_pre$start.last <- as.numeric(read_pre$start.last)
          read_pre$end.first <- as.numeric(read_pre$end.first)
          read_pre$end.last <- as.numeric(read_pre$end.last)
          read <- correct_read_order(read_pre)
          ######################################################################
          # here sometimes a weird error occurs: paired RNA reads which include 
          # themselves im not sure what that means alignment failure??          
          if(isTRUE(read[['end.last']]<read[['end.first']])){
            # anyway they must be filtered out
            stop <- T
            erg <- 'read in read'
          } else {
            stop <- F
            ####################################################################
            # now the two paired reads are pasted to one read: this is 
            # done over the cigar columns:
            # we have two possibilities:
            #   1. the read is shorter than the sum of both sequenced ends,
            #      which means that the part in the middle is covered by both
            #      ends. Here bases must be discarded in order to get the 
            #      required features of the reads (cigar string and seq-string)
            #   2. or the overall length is long enough which means the part
            #      in between our sequenced sites is not sequenced and needs to
            #      be filled up with "N"
            if(isTRUE(read[['start.last']]<read[['end.first']])){
              # ... if bases must be removed
              calc_new_cigar1 <- str_match_all(read[['cigar.last']], 
                                               '\\d*[:upper:]') %>% 
                unlist() %>%
                tibble(element=.) %>%
                mutate(
                  width=str_sub(element, start = 1, end=nchar(element)-1) %>% 
                    as.numeric(),
                  counter=c(1:nrow(.))
                )
              diff <- as.numeric(read[['end.first']])-
                as.numeric(read[['start.last']])
              calc_new_cigar2 <- calc_new_cigar1 %>%
                mutate(
                  keep=sapply(counter, function(c){
                    width_sum <- calc_new_cigar1 %>%
                      filter(counter<=c) %>%
                      .$width %>% 
                      sum()
                    ifelse(width_sum > diff, T,F) %>% 
                      return()
                  })
                ) 
              open_diff <- diff-calc_new_cigar2 %>% 
                filter(keep==F) %>% 
                .$width %>% 
                sum()
              cigar_complete <- calc_new_cigar2 %>%  
                filter(keep==T) %>%
                select(-counter)%>%
                mutate(
                  counter=c(1:nrow(.))
                ) %>%
                mutate(
                  new_width = apply(.,1, function(row){
                    if(row[['counter']]==1){
                      nw <- as.numeric(row[['width']])-open_diff-1 
                      erg <- paste0(nw, 'M')
                    } else {
                      erg <- row[['element']]
                    }
                    return(erg)
                  })
                ) %>% 
                .$new_width %>%
                Reduce(function(x,y)paste0(x,y),.) %>%
                paste0(read[['cigar.first']],.)
              seq_complete <- paste0(read[['seq.first']], 
                                     str_sub(read[['seq.last']], 
                                             start = diff+2, 
                                             end=nchar(read[['seq.last']])))             
            } else {
              # ... if the read needs to be filled up
              cigar_complete <- paste0(as.numeric(read[['start.last']])-
                                         as.numeric(read[['end.first']])-1,
                                       'N') %>%
                paste0(read[['cigar.first']], 
                       ., 
                       read[['cigar.last']])
              seq_complete <- paste0(read[['seq.first']], read[['seq.last']])
            }
          } 
        }
        ########################################################################
        ########### here starts the main evaluation part
        ########################################################################
        if(stop==T){
        } else {         
          alignment_structure <- str_match_all(cigar_complete, 
                                               '\\d*[:upper:]') %>% 
            unlist() %>%
            tibble(element=.) %>%
            mutate(
              width=str_sub(element, start = 1, end=nchar(element)-1) %>% 
                as.numeric(),
              cat=str_sub(element, start = nchar(element), end=nchar(element)),
              counter=c(1:nrow(.))
            )
          alignment_structure_processed <- alignment_structure %>%
            mutate(
              start = sapply(counter, function(e){
                if(e==1){
                  s <-read[['start.first']]
                } else {
                  s <-as.numeric(read[['start.first']])+
                    sum(alignment_structure %>%
                          filter(!cat %in% c('S', 'I'),
                                 counter < e) %>%
                          .$width)
                }
                return(as.numeric(s))
              }),
              end = as.numeric(start)+as.numeric(width)-1,
              start_in_element = sapply(counter, function(e){
                if(e==1){
                  s <-read[['start.first']]
                } else {
                  s <-as.numeric(read[['start.first']])+
                    sum(alignment_structure %>%
                          filter(!cat %in% c('N','D'),
                                 counter < e) %>%
                          .$width)
                }
                return(as.numeric(s)-as.numeric(read[['start.first']])+1)
              }),
              end_in_element = as.numeric(start_in_element)+as.numeric(width)-1
            ) %>% 
            mutate(
              seq= str_sub(seq_complete, 
                           start=as.numeric(start_in_element), 
                           end=as.numeric(end_in_element))
            )        
          position <- alignment_structure_processed %>% 
            mutate(
              base1 = apply(.,1,function(x){
                extract_base(x, 
                             ref_pos1, 
                             alignment_structure_processed)})
            ) %>%
            mutate(
              ref1 = apply(.,1, function(x){
                extract_ref(x, 
                            ref_pos1, 
                            alignment_structure_processed, 
                            'base1')})
            ) %>%
            mutate(
              base2 = apply(.,1,function(x){
                extract_base(x, 
                             ref_pos2, 
                             alignment_structure_processed)})
            ) %>%  
            mutate(
              ref2 = apply(.,1, function(x){
                extract_ref(x, 
                            ref_pos2, 
                            alignment_structure_processed, 
                            'base2')})
            )
          base1 <- position %>% filter(base1!='X') %>% .$base1
          base2 <- position %>% filter(base2!='X') %>% .$base2     
          ref1 <- position %>% filter(ref1!='X') %>% .$ref1
          ref2 <- position %>% filter(ref2!='X') %>% .$ref2
          erg <- compare_final_base(base1, base2,
                                    ref1, ref2,
                                    ref_alt1, ref_alt2,
                                    ref_ref1, ref_ref2)
        }
        return(erg)
      })
    )
  return(df)
}
extract_ref <- function(row, ref_pos, data, fallback){
  if(row[['cat']]%in%c('M', 'D')&as.numeric(row[['start']])<=ref_pos&
     as.numeric(row[['end']])>=ref_pos){
    ret <- row[[fallback]]
    # if not yet at the last row of the alignment structure... 
    #check the next category
    if(row[['counter']]!=nrow(data)){
      nextcat <- data%>%filter(counter==as.numeric(row[['counter']])+1) %>% 
        .$cat
      if(nextcat=='D'){
        ret <- rep('D', 
                   as.numeric(
                     data %>%
                       filter(counter==as.numeric(row[['counter']])+1) %>% 
                       .$width)+1) %>% 
          Reduce(function(x,y)paste0(x,y),.)
      }
      if(nextcat=='I'){
        ret <- 'I'
      }
    } 
  } else {
    ret <- 'X'
  }
  return(ret)
}
extract_base <- function(row, ref_pos, data){
  if(row[['cat']]%in%c('M','D')&
     as.numeric(row[['start']])<=ref_pos&as.numeric(row[['end']])>=ref_pos){
    ret <- str_sub(row[['seq']],
                   start = ref_pos-as.numeric(row[['start']])+1,
                   end=ref_pos-as.numeric(row[['start']])+1)
    if(row[['counter']]!=nrow(data)){
      nextcat <- data %>%
        filter(counter==as.numeric(row[['counter']])+1) %>% 
        .$cat
      if(nextcat=='I'){
        insertion_start <- ret
        ret <- paste0(insertion_start, 
                      data %>%
                        filter(counter==as.numeric(row[['counter']])+1) %>% 
                        .$seq)
      }
    }
  } else {
    ret <- 'X'
  }
  return(ret)
}
compare_final_base <- function(b1, b2, r1, r2, RB1, RB2, RR1, RR2){
  if(length(b1)>0&length(b2)>0){
    # at first check if one of the mutations is an indel
    class1 <- ifelse(nchar(RB1)==1&nchar(RR1)==1,
                     'snv',
                     ifelse(nchar(RB1)!=1&nchar(RR1)==1,
                            'ins',
                            'del'))
    class2 <- ifelse(nchar(RB2)==1&nchar(RR2)==1,
                     'snv',
                     ifelse(nchar(RB2)!=1&nchar(RR2)==1,
                            'ins',
                            'del'))
    mut1_detected <- ifelse(class1=='snv',
                            ifelse(b1==RB1,
                                   T,
                                   F),
                            ifelse(nchar(b1)==nchar(RB1)&nchar(r1)==nchar(RR1),
                                   T,
                                   F))
    mut2_detected <- ifelse(class2=='snv',
                            ifelse(b2==RB2,
                                   T,
                                   F),
                            ifelse(nchar(b2)==nchar(RB2)&nchar(r2)==nchar(RR2),
                                   T,
                                   F))
    erg <- ifelse(mut1_detected,
                  ifelse(mut2_detected,
                         'both',
                         'mut1'),
                  ifelse(mut2_detected,
                         'mut2',
                         'none'))
  } else {    
    erg <- 'spanned_out'
  }
  return(erg)
}
calc_probability <- function(df, purity, tcn1, tcn2, aff_copies1, aff_copies2, 
                             class1, class2, nDNA_reads, nRNA_reads, dist, 
                             comb){
  number <- tibble(result = c('both', 'mut1', 'mut2', 'none', 
                              'read in read', 'spanned_out')) %>%
    mutate(n=sapply(result, function(x){df %>% filter(result==x) %>% 
        nrow()})) %>% 
    column_to_rownames(var='result') %>%t() %>% 
    as.data.frame()
  both <- number[['both']]
  mut1 <- number[['mut1']]
  mut2 <- number[['mut2']]
  none <- number[['none']]
  trav <- sum(as.numeric(number[['read in read']]), 
              as.numeric(number[['spanned_out']]))
  number[['none']] <- number[['none']]-sum(number[1,])*(1-as.numeric(purity))
  # now the purity is calculated out... actually all "none" detetions  count 
  # for mono-allelic inactivation.... but to get a more reliable result it 
  # would be necessary to detect reads which have both mutations
  # therefore we check if any "both" detections were found and if not, we 
  # assume that the "none" detections are artefacts or non tumor cell tissue    
  #adj_none <- number[['none']]
  number[['none']] <- ifelse(number[['none']]<0,
                             0,
                             number[['none']])
  clear_sum <- sum(none, both, mut1, mut2)
  prob_bi_allelic_occur <- ifelse(clear_sum>0,
                                  round((number[['mut1']]+
                                           number[['mut2']])/clear_sum,2),
                                  0)
  prob_mono_allelic_occur <- ifelse(clear_sum>0&number[['both']]>0,
                                    round((number[['both']]+
                                             number[['none']])/clear_sum,2),
                                    0)
  # check if all copies are hit
  # if RNA... the complete copies cannot be defined according to the reads...
  #diff_tcn <- as.numeric(tcn1)-as.numeric(tcn2)
  mtcn <- min(as.numeric(tcn1), as.numeric(tcn2))
  if(sum(mut1, mut2, both)==0){
    status <- 'null'
    info <- "phasing"
    left_wt_copies <- mtcn-max(as.numeric(aff_copies1), as.numeric(aff_copies2))
    pre_score <- 1
  } else if(prob_mono_allelic_occur>=prob_bi_allelic_occur){
    status <- 'at same'
    info <- paste("phasing",'muts on same read', sep=' -> ')
    left_wt_copies <- mtcn-max(as.numeric(aff_copies1), as.numeric(aff_copies2))
    pre_score <- 1
  } else {
    status <- 'at diff'
    left_wt_copies <- mtcn-sum(as.numeric(aff_copies1), as.numeric(aff_copies2))
    pre_score <- ifelse(left_wt_copies<0.5, 2,1)
    info <- paste("phasing",'muts on diff reads',paste(round(left_wt_copies, 2), 
                                                       'wt-copies left'), 
                  sep=' -> ')
  } 
  if(clear_sum==0){
    rel <- 'traversed reads only' %>% paste(paste('unclear:',info), ., 
                                            paste(pre_score, 'pts'), sep=' -> ')
  } else if(none==clear_sum){
    rel <- 'muts not found at any overlapping read' %>% 
      paste(paste('unclear:',info), ., paste(pre_score, 'pts'),sep=' -> ')
  } else if(sum(number[['mut1']], number[['mut2']])==1&status=='at diff'){
    rel <- 'no both & only one with one' %>% 
      paste(paste('unclear:',info), ., paste(pre_score, 'pts'),sep=' -> ')
  } else if(abs(prob_bi_allelic_occur-prob_mono_allelic_occur)<0.2){
    rel <- 'similar probability' %>% 
      paste(paste('unclear:',info), ., paste(pre_score, 'pts'),sep=' -> ')
  } else if(sum(number[['mut1']], number[['mut2']])!=0&number[['both']]!=0){
    rel <- 'found reads with both and reads with only one' %>% 
      paste(paste('unclear:',info), ., paste(pre_score, 'pts'),sep=' -> ')
  } else if(status=='at diff'&pre_score==1&round(mtcn)<=2){
    rel <- 'muts at different copies and tcn<=2, but left wt-copies >= 0.5' %>% 
      paste(paste('unclear:',info), ., paste(pre_score, 'pts'),sep=' -> ')
  } else {
    rel <- paste(paste(green('clear'),info,sep=': '), 
                 paste(pre_score, 'pts'),sep=' -> ')
  }
  return(tibble(
    comb=comb,
    dist=dist,
    status=status,
    score=pre_score,
    wt_cp=round(left_wt_copies,2),
    tcn=round(mtcn,2),
    DNA_rds=nDNA_reads,
    RNA_rds=nRNA_reads,
    both=both,
    none=none,
    mut1=mut1,  
    mut2=mut2,
    travrsd=trav,
    info=rel
  )
  )
}
phase <- function(opt, df_gene){
  # now a dataframe with all necesarry combination of input mutations is created
  all_combinations <- outer(df_gene$pos, df_gene$pos, `paste`) %>%
    .[which(upper.tri(.))] %>%
    as.data.frame() %>% 
    set_names(nm='raw') %>%
    mutate(
      pos1=str_split(raw, ' ') %>%map_chr(.,1) %>% as.numeric(),
      pos2=str_split(raw, ' ') %>%map_chr(.,2) %>% as.numeric(),
      distance = abs(pos1-pos2)
    ) %>%
    # add all information about every mutation to the data frame
    merge(df_gene %>% set_names(nm=names(df_gene) %>% paste0(.,'1')), 
          by='pos1', all.x=T) %>%
    merge(df_gene %>% set_names(nm=names(df_gene) %>% paste0(.,'2')), 
          by='pos2', all.x=T) %>%
    .[order(.$distance),] %>%
    mutate(comb_id=paste(mut_id1, mut_id2, sep='-')) 
  final_combinations <-  apply(all_combinations, 1, function(line){
    dist <- abs(as.numeric(line[['pos1']])-as.numeric(line[['pos2']]))
    if(line[['chr1']]==line[['chr2']]){
      dna_bam <- prepare_raw_bam_file(opt$bam_file, line[['chr1']], 
                                       line[['chr2']], line[['pos1']], 
                                       line[['pos2']])
      if(!is.null(opt$rna_bam)){
         rna_bam <- prepare_raw_bam_file(opt$rna_bam, line[['chr1']], 
                                          line[['chr2']], line[['pos1']], 
                                          line[['pos2']])
      } else {
        rna_bam <- tibble()
      }
      bam <- bind_rows(dna_bam, rna_bam)
    } else {
      bam <- tibble()
      info <- 'muts at diff chrom -> no overlapping DNA-reads'
    }
    if(nrow(bam)==0){
      RESULT <- tibble(
        comb=line[['comb_id']],
        dist=dist,
        status="null",
        score=1,
        wt_cp=min(as.numeric(line[['tcn1']]), as.numeric(line[['tcn2']])) %>% 
          round(2),
        tcn=min(as.numeric(line[['tcn1']]), as.numeric(line[['tcn2']])) %>% 
          round(2),
        DNA_rds=0,
        RNA_rds=0,
        both=0,
        none=0,
        mut1=0,  
        mut2=0,
        travrsd=0,
        info=ifelse(!is.null(opt$rna_bam),
                    'unclear: no ovrlp. reads -> 1 pts',
                    'unclear: no ovrlp. reads & no RNA file -> 1 pts'),
      )
    } else {
      pre_RESULT <- verify_reads(bam, line)
      RESULT <- calc_probability(pre_RESULT, 
                                 line[['purity1']], 
                                 line[['tcn1']],
                                 line[['tcn2']],
                                 line[['aff_cp1']],
                                 line[['aff_cp2']],
                                 line[['class1']],
                                 line[['class2']],
                                 nrow(dna_bam),
                                 nrow(rna_bam),
                                 dist,
                                 comb=line[['comb_id']])
    }
    return(RESULT %>% mutate(comb=line[['comb_id']]))
  }) %>%
    Reduce(function(x,y)rbind(x,y),.) %>%
    return()
}
#######################
## secondary methods ##
#######################
get_classification <- function(data){
  case_when(
    data$class=="homdel" ~ "homdel",
    nchar(data$alt)==nchar(data$ref) ~ "snv",
    nchar(data$alt)>nchar(data$ref) ~ "ins",
    nchar(data$alt)<nchar(data$ref) ~ "del",
  ) %>%
    return()
}
get_mutational_signatures <- function(yapsa_file, seq_method){
  if(is.na(yapsa_file)){
    return(NA)
  } else {
    # select the right normalization pattern according to the sequencing method
    # (X10 == WGS)
    filter_patt <- ifelse(seq_method=='WGS', '_abs', '_norm')
    df <- yapsa_file %>%
      # read yapsa file
      read_delim("\t", escape_double = FALSE, col_names = TRUE,
                 trim_ws = TRUE) %>%
      as.data.frame() %>%
      # filter the rows according to the normalisation method
      filter(str_detect(sample, filter_patt),
             # and the signatures (AC, SBSB and ID)
             # this step must be excluded as soon as all samples had ayapsa 
             #rerun at the new version by marc...
             # right now I need to remove all PCAWG artifact signatures.... 
             #because otherwise older yapsa files will cause errors
             !str_detect(sample, "PCAWG_artif")
      ) %>%
      mutate(
        # check if the confidence interval contains zero ... TRUE == low conf,
        # FALSE == high conf
        confidence = ifelse(lower < 0 & upper > 0, 'low',
                            'high')
      ) %>% 
      # select required columns
      select(sig, sample, exposure, norm_exposure, confidence) %>%
      # round exposure
      mutate(
        #norm_exposure = round(norm_exposure, 2),
        exposure = round(exposure, 2),
      )
    n_tot <- df %>% filter(sig=='total') %>% pull(exposure) %>% first()
    # for AC signatures also the value with artifact inclusion are required
    # therefore the AC sigts are filtered from the main df and the 'Valid'
    # is stored in another df
    # as well as the 'Artif'
    df_AC_A <- df %>% filter(str_detect(sig, 'AC'),
                             sample == paste0('Artif', filter_patt)) %>%
      # the artificial exposure gets an _A suffix
      select(sig, artif_norm_exp = norm_exposure, artif_exp=exposure)
    # the two frames are merged that each signature has a single row
    df_fin <- df %>% filter(!sig=='total',
                            !sample == paste0('Artif', filter_patt)) %>%
      left_join(df_AC_A, by='sig') %>%
      # if artifact exposure is NA -> signature vanished when artifacts were
      # included... a column for that inforamtion is created
      mutate(artif = ifelse(str_detect(sig, 'AC'),
                            ifelse(is.na(artif_norm_exp),
                                   T, F),
                            NA)
      ) %>%
      select(sig, 
             norm=norm_exposure, 
             abs=exposure,
             conf=confidence,
             artif,
             artif_exp,
             artif_norm_exp)
    return(list(df_signatures=df_fin, 
                total_mutations=n_tot))
  }
}
get_tcn <- function(G2_mutated_genes, cnv_file, G2_position, G2_chromosome){
  if(!is.na(G2_mutated_genes)){
    df <- tibble(
      pos=G2_position %>% str_split(', ') %>% unlist(),
      chr=G2_chromosome %>% str_split(', ') %>% unlist()
    ) %>%
      apply(.,1, function(row){
        POS <- as.numeric(row[['pos']])
        CHR<- row[['chr']]
        chr_col <- names(
          cnv_file[which(str_detect(names(cnv_file), 'chromosome'))])
        ret <- cnv_file %>% 
          filter(as.numeric(start)<=POS&
                   as.numeric(end)>=POS&
                   select(.,all_of(chr_col)) == CHR) %>%
          pull(tcnMean)
        if(length(ret)==0){
          return("locus not found in cnv file")
        } else {
          return(ret)
        }
      }) %>% 
      as.numeric() %>%
      return()
  } else {
    return(NA)
  }    
}
get_LOH_info <- function(G2_mutated_genes, cnv_file, G2_position, 
                         G2_chromosome){
  if(!is.na(G2_mutated_genes)){
    patt <- ifelse('CNA.type' %in% names(cnv_file),
                   'CNA.type',
                   ifelse('type' %in% names(cnv_file),
                          'type',
                          ifelse('GNL' %in% names(cnv_file),
                                 'GNL',
                                 'error')))
    df <- tibble(
      pos=G2_position %>% str_split(', ') %>% unlist(),
      chr=G2_chromosome %>% str_split(', ') %>% unlist()
    ) %>% 
      apply(.,1, function(row){
        POS <- as.numeric(row[['pos']])
        CHR<- row[['chr']]
        chr_col <- names(
          cnv_file[which(str_detect(names(cnv_file), 'chromosome'))])
        ret <- cnv_file %>% 
          filter(as.numeric(start)<=POS&
                   as.numeric(end)>=POS&
                   select(.,all_of(chr_col)) == CHR) %>%
          .[[patt]]
        if(length(ret)==0){
          return("locus not found in cnv file")
        } else {
          return(ret %>% str_detect('LOH') %>% ifelse('LOH', 'HZ') %>% 
                   str_replace_na('HZ'))
        }
      }) %>%
      return()
  } else {
    return(NA)
  } 
}
get_match <- function(CHR, POS, temp){
  n <- temp %>% rowwise() %>% filter(chr==CHR&between(POS, start, end)) %>%
    pull(gene) 
  if(isTRUE(nchar(n)!=0)){return(n)}else{return(NA)}
}
extract_all_somatic_variants_of_sample <- function(opt, mut_type){
  if(nrow(opt[[paste0("somatic_", mut_type)]])==0){
    return(NULL)
    } else {
      if(mut_type=="snv"){
      calc_maf <- calc_maf_snv
      chr_col <- names(opt[[paste0("somatic_", mut_type)]]) %>% 
        .[which(str_detect(., 'CHROM|#CHROM|i#CHROM'))]
      info <- 'INFO'
    } else {
      chr_col <- '#CHROM'
      info <- names(opt[[paste0("somatic_", mut_type)]]) %>% 
        .[which(str_detect(., paste0('TUMOR|tumor|metastasis|','sample_', 
                                     opt$pid, '-', '[MT]')))] %>% .[1]
      calc_maf <- calc_maf_indel
    }
    variant_table <- opt[[paste0("somatic_", mut_type)]] %>% 
      dplyr::rename(BASE=all_of(info),
                    CHROM=all_of(chr_col)) %>% 
      rowwise() %>%
      mutate(
        maf=calc_maf(BASE),
        gene=get_match(CHROM, POS, opt$ref_gen_somatic)) %>%
      filter(!is.na(gene))
    if(nrow(variant_table)==0){
      return(NULL)
    } else {
      return(variant_table %>%    select(gene,
                                         class=EXONIC_CLASSIFICATION,
                                         site=ANNOVAR_FUNCTION,
                                         chr=CHROM,
                                         pos=POS,
                                         ref=REF,
                                         alt=ALT,
                                         maf))
    }
  }
}
verify_germline_mutations <- function(class, pos, chr, somatic_snv, 
                                      somatic_indel){
  if(!str_detect(class, 'snv')){
    ff <- somatic_indel %>%
      filter(POS==pos, select(.,1)==chr)
  } else {
    ff <- somatic_snv %>%
      filter(POS==pos, select(.,1)==chr)
  }
  rel <- ifelse(nrow(ff)!=0,
                'mutation found in somatic file -> will be excluded',
                'clear') %>% return()
}
get_HRD_LST <- function(file){
  if(is.na(file)){
    return(NA)
  } else {
    hrd_lst <- read_delim(file, delim ="\t", escape_double = FALSE,
                          col_names = TRUE, trim_ws = TRUE, 
                          col_types = cols(.default = "c"))
    # pre define hrd and lst as NA to avoid problems
    hrd <- 'not found'
    lst <- 'not found'
    tai <- 'not found'
    hrd_patt=c('numberHRDSmoothCentromereReduced',
               'numberHRDSmoothReduced',
               'numberHRDSmooth',
               'numberHRD') %>% 
      sapply(function(x){x %in% names(hrd_lst)}) %>%
      .[.==T] %>% 
      names() %>% 
      first()
    lst_patt=c('LSTCentromereReduced',
               'LSTReduced',
               'LST') %>% 
      sapply(function(x){x %in% names(hrd_lst)}) %>%
      .[.==T] %>% 
      names() %>% 
      first()
    # now we just choose the right columns... if the prio 1 columns does not
    # exist (old data), we just use the second prio column
    hrd <- hrd_lst[[hrd_patt]] %>% as.numeric()
    lst <- hrd_lst[[lst_patt]] %>% as.numeric()
    tai <- hrd_lst[["TAI"]] %>% as.numeric()
    if(length(tai)==0){
      tai <- "not_found"
    }
    # we paste them and return them
    return(c(hrd, lst, tai, hrd_patt, lst_patt) %>% 
             set_names(nm=c('HRD', 'LST', "TAI", "colname_HRD", "colname_LST")))
  }
}
get_homdels <- function(cnv_file, gene_annotation){
  patt <- ifelse('CNA.type' %in% names(cnv_file),
                 'CNA.type',
                 ifelse('type' %in% names(cnv_file),
                        'type',
                        ifelse('GNL' %in% names(cnv_file),
                               'GNL',
                               'error')))
  homdels_f <- cnv_file %>% 
    filter(select(.,all_of(patt))=='HomoDel',
           !str_detect(as.character(tcnMeanRaw), '0\\.000'))
  if(nrow(homdels_f)==0){return(NA)}else{
    homdels <- homdels_f %>% 
      apply(.,1,function(cnv){
        gene_annotation %>%
          filter(
            chr==cnv[[which(str_detect(names(cnv), 'chromosome'))]],
            !as.numeric(start)>as.numeric(cnv[['end']]),
            !as.numeric(end)<as.numeric(cnv[['start']]),
          ) %>%
          .$gene %>% 
          return()
      }) %>% 
      compact() %>% 
      unlist()
    if(length(homdels)!=0){
      return(homdels)
    } else {
      return(NA)
    }
  }
}
get_sex <- function(pid, file, seqm){
  if(!is.na(file)){
    if(seqm=='WES'){
      str_split(file, 'segment') %>% map_chr(.,1) %>%
        list.files(pattern = 'Identified') %>% 
        str_match('Male|Female') %>% 
        tolower() %>% as.character() %>%
        return()
    } else {
      str_split(file, paste0(pid, '_comb_pro')) %>% map_chr(.,1) %>%
        list.files(pattern = 'sex', full.names = T) %>% 
        read_lines() %>% as.character() %>%
        return()
    }
  } else {
    return(NA)
  }
}
get_methylation <- function(opt){
  meth_patt <- pull(opt$ref_gen_somatic, gene) %>% 
    lapply(function(x){return(paste0(x,";|", x, "$"))}) %>% 
    unlist() %>% paste(collapse = "|")
  meth_match <- pull(opt$ref_gen_somatic, gene) %>% paste(collapse = "|")
  
  opt$methylation_file %>% left_join(opt$promoter, 
                                     by=c("sample_id"="probeID")) %>%
    filter(!is.na(gene)) %>%
    filter(str_detect(gene, meth_patt)) %>%
    mutate(gene=str_match(gene, meth_match)) %>%
    select(2, gene) %>% 
    set_names(nm=c("beta", "gene")) %>%
    group_by(gene) %>%
    summarize(beta_mean=mean(beta)) %>%
    select(gene, beta=beta_mean) %>%
    return()
}
extract_all_homdels_of_sample <- function(cnv_file, gene_annotation){
  patt <- ifelse('CNA.type' %in% names(cnv_file),
                 'CNA.type',
                 ifelse('type' %in% names(cnv_file),
                        'type',
                        ifelse('GNL' %in% names(cnv_file),
                               'GNL',
                               'error')))
  homdels_f <- cnv_file %>% 
    filter(select(.,all_of(patt))=='HomoDel',
           !str_detect(as.character(tcnMeanRaw), '0\\.000'))
  if(nrow(homdels_f)==0){return(NULL)}else{
    homdels <- homdels_f %>% 
      apply(.,1,function(cnv){
        gene_annotation %>%
          filter(
            chr==cnv[[which(str_detect(names(cnv), 'chromosome'))]],
            !as.numeric(start)>as.numeric(cnv[['end']]),
            !as.numeric(end)<as.numeric(cnv[['start']]),
          ) %>%
          .$gene %>% 
          return()
      }) %>% 
      compact() %>% 
      unlist()
    if(length(homdels)!=0){
      return(tibble(gene=homdels,
                    class='homdel',
                    site='',
                    chr = '',
                    pos = '',
                    ref= '',
                    alt= '',
                    maf='',
                    tcn = 0,
                    zygost='HOMDEL',
                    aff_cp=0,
                    wt_cp=0,
                    pre_info='homdel -> all copies affected',)
      )
    } else {
      return(NULL)
    }
  }
}
prepare_somatic_variant_table <- function(opt){
  df_som_pre <- rbind(extract_all_somatic_variants_of_sample(opt,
                                        'snv'),
                      extract_all_somatic_variants_of_sample(opt,
                                        'indel'))
  if(!is.null(df_som_pre)){
    df_som_pre %>% 
      rowwise() %>%
      mutate(
        tcn=get_tcn(gene, opt$somatic_cnv, pos, chr),
        zygost=get_LOH_info(gene, opt$somatic_cnv, pos, chr),
        aff_cp = aff_som_copies(maf, tcn, opt$purity, chr, opt$sex),
        wt_cp = tcn-aff_cp,
        pre_info = ifelse(isTRUE(str_detect(zygost,'LOH')),
                     ifelse(wt_cp<=0.5,
                       ifelse(wt_cp<=0.4,
                         'somatic-variant -> LOH detected -> left wt-copies' %+%
                           ' <= 0.5 -> all copies affected',
                         'somatic-variant -> LOH detected -> ' %+% 
                           red('left wt-copies <= 0.5') %+% 
                           ' -> all copies affected'),
                       ifelse(wt_cp>0.6,
                             'somatic-variant -> LOH detected -> ' %+%
                               'left wt-copies > 0.5 -> not all affected',
                             'somatic-variant -> LOH detected -> ' %+% 
                               red('left wt-copies > 0.5') %+% 
                               ' -> not all affected')),
                     ifelse(chr=="X"&opt$sex=="male",
                       ifelse(wt_cp<=0.5,
                         "somatic-variant -> chrX & male -> " %+%
                           "left wt-copies <= 0.5 -> all copies affected",
                         'somatic-variant -> chrX & male -> left wt-copies ' %+%
                           '> 0.5 -> not all affected'),
                       'somatic-variant -> no LOH detected -> not all affected')
        )
      ) %>%
      return()
  } else {
    return(NULL)
  }
}
#####################
## primary methods ##
#####################
load_input_data <- function(opt){
  if(!is.null(opt$rna_bam)){
    if(opt$rna_bam=="NA"){
      opt$rna_bam <- NULL
    }
  }
  opt$gencode <- read_tsv(opt$ref_gen,
    col_types = cols(.default = "c"))
  opt$ref_gen_germline_topart_relevant <- read_tsv(opt$germline_genes,
    col_types = cols(.default = "c")) %>%
    pull(gene)
  opt$ref_gen_somatic_topart_relevant <- read_tsv(opt$somatic_genes,
    col_types = cols(.default = "c")) %>%
    pull(gene)
  if(isTRUE(!is.null(opt$all_genes)&opt$all_genes==T)){
    opt$ref_gen_somatic <- opt$gencode
    opt$ref_gen_germline <- opt$gencode
    opt$promoter <- opt$gencode %>% extract_promoters_from_gencode()
  } else {
    opt$ref_gen_somatic <- opt$gencode %>% 
      filter(gene %in% opt$ref_gen_somatic_topart_relevant)
    opt$ref_gen_germline <- opt$gencode %>% 
      filter(gene %in% opt$ref_gen_germline_topart_relevant)
    opt$promoter <- opt$ref_gen_somatic %>% extract_promoters_from_gencode()
  }
  opt$g1_genes <- pull(opt$ref_gen_somatic, gene) %>% 
    c("XXXXXX", ., "XXXXXX") %>% paste(collapse = ", ")
  opt$g2_genes <- pull(opt$ref_gen_germline, gene) %>% 
    c("XXXXXX", ., "XXXXXX") %>% paste(collapse = ", ")
  opt$cohort_analysis <- ifelse(!is.null(opt$cohort_analysis),
                                as.logical(opt$cohort_analysis),
                                FALSE)
  opt$pid <- ifelse(!is.null(opt$pid),
                    opt$pid,
                    opt$bam_file %>% 
                      str_split('view-by-pid/') %>% 
                      unlist() %>% 
                      last() %>% 
                      str_split('/') %>% 
                      unlist() %>% 
                      first())
  opt$seq_method <- ifelse(str_detect(opt$bam_file, 'whole_genome_sequencing'),
                           'WGS',
                           'WES')
  opt$sex <- get_sex(opt$pid, opt$cnv_file, opt$seq_method)
  opt$id <- paste(opt$pid, opt$sample, opt$seq_method, sep='_')
  invisible(capture.output(
    opt$germline_variants <- read.vcfR(opt$charger_file) 
  ))
  opt$somatic_cnv <- read_tsv(opt$cnv_file, col_types = cols(.default = "c"))
  opt$somatic_snv <- read_tsv(opt$somatic_snv, col_types = cols(.default = "c"))
  opt$somatic_indel <- read_tsv(opt$somatic_indel, 
                                col_types = cols(.default = "c")) 
  return(opt)
}
estimate_signature_score <- function(opt){
  signatures = get_mutational_signatures(opt$yapsa_file, opt$seq_method)
  opt$total_mutations <- signatures[['total_mutations']]
  opt$df_signatures <- signatures[['df_signatures']]
  AC3_data <- signatures[['df_signatures']] %>% 
    filter(sig=='AC3')
  if(as.numeric(signatures[['total_mutations']]) >=25){
    if(nrow(AC3_data)!=0){
      if(AC3_data$conf=='high'){
        opt$signature_score <- 2
        opt$signature_scoring <- paste(
          'signature analysis based on', 
          signatures[['total_mutations']],
          'snvs -> AC3 detected -> zero not in confidence interval -> 2 pts')
        opt$signature_add <- ifelse(
          isTRUE(AC3_data$artif),
          red('AC3 disappears when artifact signatures are included!'),
          '')
      } else {
        opt$signature_score <- 1
        opt$signature_scoring <- paste(
          'signature analysis based on', 
          signatures[['total_mutations']], 
          'snvs -> AC3 detected -> zero in confidence interval -> 1 pts')
        opt$signature_add <- ifelse(
          isTRUE(AC3_data$artif),
          'AC3 disappears when artifact signatures are included!',
          '')
      }
      AC3_out <- AC3_data %>% 
        as.character() %>% 
        set_names(nm=paste0('AC3_', names(AC3_data))) %>% .[-1]
    } else {
      AC3_out <- NULL
      opt$signature_score <- 0
      opt$signature_scoring <- paste('signature analysis based on', 
                                     signatures[['total_mutations']], 
                                     'snvs -> no AC3 detected -> 0 pts')
      opt$signature_add <- ''
    }
  } else {
    opt$signature_score <- 0
    opt$signature_scoring <- paste('signature analysis based on', 
                                   signatures[['total_mutations']], 
                                   'snvs -> 0 pts')
    if(nrow(AC3_data)==0){
      AC3_out <- NULL
      opt$signature_add <- ''
    } else {
      opt$signature_add <- paste0('but ', 
                                  AC3_data$conf, 
                                  '-confidence AC3 detection')
      AC3_out <- AC3_data %>% 
        as.character() %>% 
        set_names(nm=paste0('AC3_', names(AC3_data))) %>% .[-1]
    }
  }
  return(c(opt, AC3_out))
}
estimate_rearrangement_score <- function(opt){
  HRD_LST_raw = get_HRD_LST(opt$hrd_file)
  opt$rearrangement_score <- case_when(
    sum(as.numeric(HRD_LST_raw[1:2])) > 20 ~ 2,
    sum(as.numeric(HRD_LST_raw[1:2])) > 10 ~ 1,
    TRUE ~ 0
  )
  opt$rearrangement_scoring <- paste(
    'HRD:', HRD_LST_raw[1], 
    'LST:', HRD_LST_raw[2], 
    ifelse(sum(as.numeric(HRD_LST_raw[1:2]))%in% c(9,10,11,12,19,20,21,22),
           red(paste('-> sum:', sum(as.numeric(HRD_LST_raw[1:2])))),
           paste('-> sum:', sum(as.numeric(HRD_LST_raw[1:2])))), 
    '->', opt$rearrangement_score, 'pts')
  pp <- opt$cnv_file %>% 
    str_match('comb_pro_extra(\\d\\.*\\d*)_(\\d\\.*\\d*).txt$') %>% 
    .[2:3] %>% as.numeric() %>% set_names(c('ploidy', 'purity'))
  return(c(opt, HRD_LST_raw, pp))
}
estimate_germline_score <- function(opt){
  samples_used_by_charger <- colnames(opt$germline_variants@gt)[2:3] %>% 
    set_names(nm=c('charger_control_sample', 'charger_tumor_sample'))
  df_germ_raw <- opt$germline_variants@gt %>%
    as_tibble() %>%
    filter(CharGer_Classification %in% c('Pathogenic','Likely Pathogenic'))
  if(nrow(df_germ_raw)==0){
    opt$germline_score <- 0
  } else {
    df_germ_pre <- df_germ_raw %>%
    rowwise() %>%
      mutate(
        pos=ifelse(str_detect(HGVSg, 'del'),
                   str_match(HGVSg, 'g\\.\\d*') %>% str_split('\\.') %>% 
                     map_chr(.,2) %>% as.numeric() %>% -1 %>% as.character(),
                   str_match(HGVSg, 'g\\.\\d*') %>% str_split('\\.') %>% 
                     map_chr(.,2) %>% as.character()),
        chr=HGVSg %>% str_match('^\\d*') %>% unlist() %>% as.character(),
        gene=get_match(chr, pos, opt$ref_gen_germline)) %>%
      filter(!is.na(gene)) 
    if(nrow(df_germ_pre)!=0){ 
      opt$df_germ_unfiltered <- df_germ_pre %>% 
        mutate(
          ref = case_when(
            str_detect(HGVSg,'ins') ~ 'R',
            str_detect(HGVSg,'del') ~ str_match(HGVSg, 'del[:upper:]*') %>% 
              str_replace('del', '') %>% paste0('R',.),
            TRUE ~ str_match(HGVSg, '[:upper:]>[:upper:]') %>% 
              str_split('>') %>% 
              unlist() %>% first(),
          ),
          alt = case_when(
            str_detect(as.character(HGVSg),'ins') ~ str_match(
              HGVSg, 
              'ins[:upper:]*') %>% 
              str_replace('ins', '')%>% paste0('R',.),
            str_detect(as.character(HGVSg),'del') ~ 'R',
            TRUE ~ str_match(HGVSg, '[:upper:]>[:upper:]') %>% 
              str_split('>') %>% unlist() %>% last()
          ),
          class=case_when(
            str_detect(HGVSg,'ins') ~ 'germline_insertion',
            str_detect(HGVSg,'del') ~ 'germline_deletion',
            TRUE ~ 'germline_snv'
          ),
          maf=as.numeric(Tumor_VAF),
          rel=verify_germline_mutations(class, pos, chr, opt$somatic_snv, 
                                        opt$somatic_indel),
          tcn=get_tcn(HUGO_Symbol, opt$somatic_cnv, pos, chr),
          zygost=get_LOH_info(HUGO_Symbol, opt$somatic_cnv, pos, chr),
          aff_cp = aff_germ_copies(maf, tcn, opt$purity, chr, opt$sex, 1),
          wt_cp = tcn-aff_cp,
          pre_info = ifelse(
            isTRUE(zygost=='LOH'),
            ifelse(maf>=0.5,
              ifelse(maf>=0.6,
                'germline-variant -> LOH detected -> ' %+% 
                  'AF >= 0.5 -> all copies affected',
                'germline-variant -> LOH detected -> ' %+% red('AF >= 0.5') %+% 
                  ' -> all copies affected'),
              ifelse(maf<=0.4,
                     'germline-variant -> LOH detected ->' %+% 
                       ' AF < 0.5 -> variant lost in tumor',
                     'germline-variant -> LOH detected -> ' %+% red('AF < 0.5') 
                     %+% ' -> variant lost in tumor')),
            ifelse(maf>=0.8,
                   'germline-variant -> no LOH detected -> ' %+% 
                     red('!! but AF >= 0.8 ??') %+% ' -> not all affected',
                   'germline-variant -> no LOH detected -> not all affected')),
          site = 'germ_NA',
      ) %>%
        select(gene,
               class,
               site,
               chr,
               pos,
               ref,
               alt,
               maf,
               tcn,
               zygost,
               aff_cp,
               wt_cp,
               ACMG_class=CharGer_Classification,
               rel,
               pre_info
        ) 
      opt$G2_raw_output <- apply(opt$df_germ_unfiltered %>% 
                                   select(G2_raw_mutated_genes=gene, 
                                          G2_ACMG_class=ACMG_class, 
                                          G2_reliability=rel), 
                                 2, 
                                 function(x)paste(x,collapse = ', '))
      opt$df_germ <- opt$df_germ_unfiltered  %>% filter(rel=='clear')%>% 
        select(-rel, -ACMG_class)
      opt$df_germ_topart_relevant <- opt$df_germ %>% 
        filter(gene %in% opt$ref_gen_germline_topart_relevant)
      if(nrow(opt$df_germ)==0){
      } else {
        if(nrow(opt$df_germ_topart_relevant)==0){
          opt$germline_score <- 0
          opt$df_germ_topart_relevant <- NULL
        } else {
          opt$germline_score <- 1
          opt$G2_final_output <- apply(
            opt$df_germ_unfiltered %>% 
              filter(gene %in% opt$ref_gen_germline_topart_relevant) %>% 
              select(G2_final_mutated_genes=gene, 
                     G2_final_ACMG_class=ACMG_class), 
            2, 
            function(x)paste(x,collapse = ', '))
        }
      }  
    } else {
      opt$germline_score <- 0
    }
  }  
  return(opt)
}
estimate_somatic_score_for_each_gene <- function(GENE, df_all_mutations, opt){
  #print(GENE)
  pre_df_gene <- df_all_mutations %>% filter(gene==GENE) %>%
    rownames_to_column('mut_id') %>%  mutate(mut_id=paste0('m', mut_id))
  df_gene <- pre_df_gene %>% 
    # all germline variants which are lost in the tumor are excluded
    filter(!str_detect(pre_info, 'variant lost in tumor'),
           !str_detect(zygost, "locus not found"))
  if(nrow(df_gene)==0){
    return(list(pre_df_gene, NULL))
  } else {
    log_obj <- list(gene=GENE, n_mut=nrow(df_gene))
    all_comb <- NULL
    # at first check if any of the mutations already affects all alleles -> 2pts            
    if(str_detect(paste(df_gene$pre_info, collapse = ' '), 
                  'all copies affected')){
      concern_info <- df_gene %>% 
        filter(str_detect(pre_info, "all copies affected")) %>% .[1,]
      class <- get_classification(concern_info)
      origin <- ifelse(str_detect(concern_info$class, "germ"),
                       "germ-",
                       "som-")
      log_obj[['score']] <- 2
      log_obj[['info']] <- 'all-aff in pre-evaluation -> 2 pts' %>% 
        paste(green('clear'),.,sep=': ')                  
      log_obj[['tcn']] <- concern_info$tcn
      log_obj[['aff_cp']] <- concern_info$aff_cp %>% as.numeric() 
      log_obj[['class']] <- paste0(origin, class)
    } else if(nrow(df_gene)==1){
      # now check if that mutation is germline...                 
      if(str_detect(df_gene$class, 'germline')){
        log_obj[['score']] <- 0
        log_obj[['info']] <- paste(
          green('clear'),
          'germline mutation has no influence to somatic score -> 0 pts',
          sep=': ')
        log_obj[['tcn']] <- df_gene$tcn
        log_obj[['aff_cp']] <- df_gene$aff_cp %>% as.numeric()
        log_obj[['class']] <- paste0("germ-", get_classification(df_gene))
      } else {
        log_obj[['score']] <- 1
        log_obj[['info']] <- 'affects not all copies -> 1 pts' %>% 
          paste(green('clear'),.,sep=': ')
        log_obj[['tcn']] <- df_gene$tcn
        log_obj[['aff_cp']] <- df_gene$aff_cp %>% as.numeric()
        log_obj[['class']] <- paste0("som-", get_classification(df_gene))
      }
    } else {  
      all_comb <- phase(opt,
                        df_gene)
      if(nrow(all_comb)==1){
        log_obj[['score']] <- all_comb$score
        log_obj[['info']] <- all_comb$info
        log_obj[['tcn']] <- all_comb$tcn
        log_obj[['aff_cp']] <- all_comb$tcn-all_comb$wt_cp %>% as.numeric()
        log_obj[['class']] <- paste("comb:", all_comb$comb)
      } else if(max(all_comb$score)==2){
        most_important_comb <- all_comb %>% 
          filter(wt_cp==min(all_comb$wt_cp)) %>%
          filter(rownames(.)==1)
        log_obj[['score']] <- 2
        log_obj[['info']] <- paste(
          green('clear'),
          'at least one of the combinations of muts affects all copies',
          sep=': ') 
        log_obj[['tcn']] <- most_important_comb$tcn
        log_obj[['aff_cp']] <- 
          most_important_comb$tcn-most_important_comb$wt_cp %>% 
          as.numeric()
        log_obj[['class']] <- paste("comb:", most_important_comb$comb)
      } else if(str_count(paste(all_comb$status, collapse=' '),
                          'at diff')==nrow(all_comb)&
                nrow(all_comb)>=mean(as.numeric(df_gene$tcn))){
        most_important_comb <- all_comb %>% 
          filter(wt_cp==min(all_comb$wt_cp)) %>%
          filter(rownames(.)==1)
        log_obj[['score']] <- 1
        log_obj[['info']] <- paste(
          'unclear: 2 pts are very likely because all',nrow(all_comb),
          'mutations are at different copies and tcn=', 
          mean(as.numeric(df_gene$tcn)),
          'but as said in the meeting we are not sure if',
          'tumor is subclonal and therefore only 1 pts')
        log_obj[['tcn']] <- most_important_comb$tcn
        log_obj[['aff_cp']] <- 
          most_important_comb$tcn-most_important_comb$wt_cp %>% 
          as.numeric()
        log_obj[['class']] <- paste("comb:", most_important_comb$comb)
      } else if(str_detect(paste(all_comb$status, collapse=' '),
                           'unclear: muts at different copies and tcn=2')){
        most_important_comb <- all_comb %>% 
          filter(str_detect(status, 
                            'unclear: muts at different copies and tcn=2')) %>%
          filter(wt_cp==min(all_comb$wt_cp)) %>%
          filter(rownames(.)==1)
        log_obj[['score']] <- 1
        log_obj[['info']] <- paste(
          "unclear: 2 pts are very likely because two mutations were detected",
          "at different reads and the tcn is 2 but as said in the meeting", 
          "we are not sure if tumor is subclonal and therefore only 1 pts")
        log_obj[['tcn']] <- most_important_comb$tcn
        log_obj[['aff_cp']] <- 
          most_important_comb$tcn-most_important_comb$wt_cp %>% 
          as.numeric()
        log_obj[['class']] <- paste("comb:", most_important_comb$comb)
      } else if(nrow(all_comb)>3){
        most_important_comb <- all_comb %>% 
          filter(wt_cp==min(all_comb$wt_cp)) %>%
          filter(rownames(.)==1)
        log_obj[['score']] <- 1
        log_obj[['info']] <- 'unclear: very high number of mutations... ' %+% 
          'probably 1 pts but please check the results above'
        log_obj[['tcn']] <- most_important_comb$tcn
        log_obj[['aff_cp']] <- 
          most_important_comb$tcn-most_important_comb$wt_cp %>% 
          as.numeric()
        log_obj[['class']] <- paste("comb:", most_important_comb$comb)
      } else {
        most_important_comb <- all_comb %>% 
          filter(wt_cp==min(all_comb$wt_cp)) %>%
          filter(rownames(.)==1)
        log_obj[['score']] <- 1
        log_obj[['info']] <- 'wt-copies left' %>% 
          paste(green('clear'),.,sep=': ')
        log_obj[['tcn']] <- most_important_comb$tcn
        log_obj[['aff_cp']] <- 
          most_important_comb$tcn-most_important_comb$wt_cp %>% 
          as.numeric()
        log_obj[['class']] <- paste("comb:", most_important_comb$comb)            
      }
      ## add methylation info
    } 
    df_reduced <- log_obj %>% as.data.frame() 
    return(list(pre_df_gene, list(df_reduced, all_comb)))
  }
}
estimate_somatic_score <- function(opt){
  opt$df_homdels <- 
    extract_all_homdels_of_sample(
      opt$somatic_cnv, 
      opt$ref_gen_somatic)
  opt$df_som <- 
    prepare_somatic_variant_table(opt)
  if(!is.null(opt$df_germ)|!is.null(opt$df_som)|!is.null(opt$df_homdels)){
    df_all_mutations <- list(opt$df_germ, opt$df_som, opt$df_homdels) %>% 
      compact() %>%
      Reduce(function(x,y)rbind(x,y),.) %>% 
      as_tibble() %>% mutate(purity=opt$purity)
    result_list <- 
      sapply(
        unique(df_all_mutations$gene), 
        estimate_somatic_score_for_each_gene, 
        df_all_mutations, 
        opt, 
        simplify = F)
    opt$pre_scoring <- lapply(result_list, nth, n=1) %>% 
      bind_rows()
    opt$pre_scoring_topart_relevant <- opt$pre_scoring %>%
      filter(gene %in% opt$ref_gen_somatic_topart_relevant)
    final_result_list <- lapply(result_list, nth, n=2) %>% compact()
    if(length(final_result_list)!=0){
      topart_relevant_list_positions <- 
        names(final_result_list)[which(names(final_result_list) %in% 
                                         opt$ref_gen_somatic_topart_relevant)]
      final_result_list_topart_relevant <- 
        final_result_list[topart_relevant_list_positions]
      opt$final_scoring <- lapply(final_result_list, first) %>% 
        compact() %>% bind_rows() %>% 
        t() %>% as.data.frame() %>% rownames_to_column('names') %>% 
        apply(.,1,correct_space_in_printed_tables) %>% as_tibble() %>% 
        set_names(nm=as.character(.[1,])) %>% .[c(2:nrow(.)),]
      opt$G1_final_output <- apply(lapply(final_result_list, first) %>% 
                                     bind_rows() %>%
                                     select(G1_final_mutated_genes=gene,
                                            G1_final_score=score,
                                            G1_final_info=info),
                                   2, 
                                   function(x)paste(x,collapse = ', '))
      opt$final_phasing_info <-  lapply(final_result_list, last) %>% 
        compact() %>% 
        Map(cbind,. , name=names(.)) %>% 
        bind_rows()
      if(length(final_result_list_topart_relevant)!=0){
        opt$phasing_info_topart_relevant <- 
          lapply(final_result_list_topart_relevant, last)
        opt$final_scoring_topart_relevant <- 
          lapply(final_result_list_topart_relevant, first) %>% 
          compact() %>% 
          bind_rows() %>% 
          t() %>% 
          as.data.frame() %>% 
          rownames_to_column('names') %>% 
          apply(.,1,correct_space_in_printed_tables) %>% 
          as_tibble() %>% 
          set_names(nm=as.character(.[1,])) %>% .[c(2:nrow(.)),]
        opt$G1_final_output_topart_relevant <- apply(
          lapply(final_result_list_topart_relevant, first) %>% 
            bind_rows() %>%
            select(G1_final_mutated_genes=gene,
                   G1_final_score=score,
                   G1_final_info=info),
          2, 
          function(x)paste(x,collapse = ', '))
        opt$final_phasing_info_topart_relevant <- 
          lapply(final_result_list_topart_relevant, last) %>% 
          compact() %>% 
          Map(cbind,. , name=names(.)) %>% 
          bind_rows()
        opt$somatic_score <- 
          max(as.numeric(opt$final_scoring_topart_relevant[['score']]))
      } else {
        opt$somatic_score <- 0
        opt$final_scoring_topart_relevant <- NULL
        opt$G1_final_output_topart_relevant <- NULL
        opt$final_phasing_info_topart_relevant <- NULL
      }
    } else {
      opt$somatic_score <- 0
    }  
  } else {
    opt$somatic_score <- 0
  }
  return(opt)
}  
estimate_topart_score <- function(opt){
  opt$TOP_ART_score <- sum(opt$signature_score,
                           opt$rearrangement_score,
                           opt$somatic_score,
                           opt$germline_score)
  opt$biomarker <- ifelse(opt$TOP_ART_score>2,
                          T,
                          F)
  return(opt)
}
analyse_methylation <- function(opt){
  # if a methylation file is provided ...
  if(!is.na(opt$meth_file)&opt$meth_file!="NA"){
    # ... load methylation file
    opt$methylation_file <- read_feather(opt$meth_file)
    # then extract required information
    opt$methylation_table <- get_methylation(opt)
    opt$methylation_table_topart_relevant <- opt$methylation_table %>%
      filter(gene %in% opt$ref_gen_somatic_topart_relevant)
  } else {
    # ... otherwise set methylation result to NULL
    opt$methylation_table <- NULL
  }
  return(opt)
}
analyse_gene_expression <- function(opt){
  if(!is.null(opt$rna_bam)){
    feature_patt <- pull(opt$ref_gen_somatic, gene)
    opt$feature_counts <- opt$rna_bam %>% 
      str_match("(.+/).*$") %>% 
      unlist() %>% 
      last() %>%
      paste0(., "featureCounts") %>% 
      list.files(., pattern="fpkm_tpm.featureCounts.tsv$", full.names = T) %>%
      read_tsv(col_types = cols(.default = "c")) %>%
      filter(name %in% feature_patt) %>%
      select(gene=name, TPM)
    opt$feature_counts_topart_relevant <- opt$feature_counts %>%
      filter(gene %in% opt$ref_gen_somatic_topart_relevant)
  } else {
    opt$feature_counts <- NULL
  }
  return(opt)
}
export_data <- function(opt){
  if(opt$cohort_analysis==T){
    if(isTRUE(opt$all_genes==T)){
      outdir_all_genes <- file.path(opt$main_output_path, "hg19_genes")
      dir.create(outdir_all_genes)
      store_table(opt, "final_phasing_info", outdir_all_genes, 
                  'phasing_info.tsv')
      store_table(opt, "final_scoring", outdir_all_genes, 
                  'evaluation_per_gene.tsv')
      store_table(opt, "pre_scoring", outdir_all_genes, 
                  'evaluation_per_variant.tsv')
      store_table(opt, "methylation_table", outdir_all_genes, 
                  'methylation_per_gene.tsv')
      store_table(opt, "feature_counts", outdir_all_genes, 
                  'expression_per_gene.tsv')
    }
    outdir_topart_genes <- file.path(opt$main_output_path, "TOP-ART_genes")
    dir.create(outdir_topart_genes)
    store_table(opt, "final_phasing_info_topart_relevant", outdir_topart_genes, 
                'phasing_info.tsv')
    store_table(opt, "final_scoring_topart_relevant", outdir_topart_genes, 
                'evaluation_per_gene.tsv')
    store_table(opt, "pre_scoring_topart_relevant", outdir_topart_genes, 
                'evaluation_per_variant.tsv')
    store_table(opt, "methylation_table_topart_relevant", outdir_topart_genes, 
                'methylation_per_gene.tsv')
    store_table(opt, "feature_counts_topart_relevant", outdir_topart_genes, 
                'expression_per_gene.tsv')
    store_table(opt, "df_signatures", opt$main_output_path, 
                'mutational_signatures.tsv')
    final_output <- opt %>% lapply(function(obj){
      if(!class(obj) %in% c("character", "numeric", 
                            "logical", "integer")|length(obj)>1){
        return(NULL)
      } else {
        return(obj)
      }
    }) %>% 
      compact() %>% 
      unlist() %>% 
      as.data.frame(nm='val') %>% 
      rownames_to_column('var') %>%
      mutate(ID=opt$id)
    write_tsv(final_output, 
              file=file.path(opt$main_output_path, 'full_data.tsv'))
    write_tsv(tibble(paste("all steps succesfully finished at script version:", 
                           opt$version)), 
              file=file.path(opt$main_output_path, paste0('success_', opt$id, 
                                                          "_V",opt$version, 
                                                          ".txt")),
              col_names=F)
  }
  return(opt) 
}
output_message <- function(opt){
  opt %>%
    prepare_germline_output_message() %>%
    prepare_somatic_output_message() %>%
    print_full_output_message()
}
###################
## main function ##
###################
main <- function(input){
  opt <- input                     %>% 
    as.list()                      %>%
    load_input_data()              %>%
    analyse_methylation()          %>%
    analyse_gene_expression()      %>%
    estimate_signature_score()     %>%
    estimate_rearrangement_score() %>%
    estimate_germline_score()      %>%
    estimate_somatic_score()       %>%
    estimate_topart_score()        %>%
    export_data()                  %>%
    output_message()
}