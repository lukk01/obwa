/*
  This work is published under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007.
  
  Copyright (c) 2013 by Margus Lukk <margus.lukk@cruk.cam.ac.uk>
  
  This work uses code from bwa-0.6.1 published under the same licence by Heng Li <lh3lh3@gmail.com>
*/
#include "bwtaln.h" // needed for gap_opt_t
#include "bntseq.h" // needed for bntseq_t

#ifndef OBWAPE_H
#define OBWAPE_H

#ifdef __cplusplus
extern "C" {
#endif

  void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2);
  void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);
  void obwape_bwa_aln_core_inner(bwt_t *bwt, bwa_seq_t *seqs, int n_seqs, const gap_opt_t *opt);
  bwa_seq_t *obwape_copy_seqs(bwa_seq_t *seqsFrom,int n_seqs,gap_opt_t *opt);
  int parse_aln(gap_opt_t *opt, char *args);
  bntseq_t *bwa_open_nt(const char *prefix);
  int obwape(int argc, char *argv[]);
  
#ifdef __cplusplus
}
#endif

#endif
