/*
  This work is published under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007.

  Copyright (c) 2013 by Margus Lukk <margus.lukk@cruk.cam.ac.uk>
  
  This work uses code from bwa-0.6.1 published under the same licence by Heng Li <lh3lh3@gmail.com>
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "utils.h"
#include "obwape.h"
#include "bwase.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.6.1-r104"
#endif
#ifndef BWASE_VERSION
#define BWASE_VERSION "0.11"
#endif

void bwa_print_sam_SQ(const bntseq_t *bns);
void bwa_print_sam_PG();

void obwase_bwa_cal_pac_pos(const bntseq_t *bns, bwt_t *bwt, int n_seqs, bwa_seq_t *seqs, int max_mm, float fnr)
{
  // adapted from bwa_cal_pac_pos from bwase.c
  // Modifications:
  // 1. Used to take const char *prefix for loading indexes, now takes pre-loaded index "bwt_t *bwt" instead.
  int i, j, strand, n_multi;
  /*
    char str[1024];
    bwt_t *bwt;
    // load forward SA
    strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
    strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
  */
  for (i = 0; i != n_seqs; ++i) {
    bwa_seq_t *p = &seqs[i];
    bwa_cal_pac_pos_core(bns, bwt, p, max_mm, fnr);
    for (j = n_multi = 0; j < p->n_multi; ++j) {
      bwt_multi1_t *q = p->multi + j;
      q->pos = bwa_sa2pos(bns, bwt, q->pos, p->len, &strand);
      q->strand = strand;
      if (q->pos != p->pos)
        p->multi[n_multi++] = *q;
    }
    p->n_multi = n_multi;
  }
  // bwt_destroy(bwt);
}

int obwase_core(char *prefix, char *fn_fa, gap_opt_t *opt, int n_occ, int skipheader) {
  // The function has been adapted/merged from bwa_sai2sam_se_core in bwase.c and bwa_aln_core in bwtaln.c

  // options from bwa_sai2sam_se_core
  extern bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);
  int i, n_seqs, tot_seqs = 0;
  // int m_aln;
  bwt_aln1_t *aln = 0;
  bwa_seq_t *seqs, *seqsC;
  bwa_seqio_t *ks;
  clock_t t;
  bntseq_t *bns, *ntbns = 0;
  // FILE *fp_sa;
  // gap_opt_t opt;
  // additional options from bwa_aln_core
  bwt_t *bwt;

  // initialization                                                                                                                                             
  bwase_initialize();
  bns = bns_restore(prefix);
  srand48(bns->seed);
  // fp_sa = xopen(fn_sa, "r");

  // load index as in bwa_aln_core
  { // load BWT
    char *str = (char*)calloc(strlen(prefix) + 10, 1);
    strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
    // following line is from bwa_cal_pac_pos (commented out in obwase_bwa_cal_pac_pos)
    strcpy(str, prefix); strcat(str, ".sa");   bwt_restore_sa(str,bwt);
    free(str);
  }

  // m_aln = 0;
  // fread(&opt, sizeof(gap_opt_t), 1, fp_sa);
  if (!(opt->mode & BWA_MODE_COMPREAD)) // in color space; initialize ntpac
    ntbns = bwa_open_nt(prefix);
  if(!skipheader) {
    bwa_print_sam_SQ(bns);
    bwa_print_sam_PG();
  }

  // set ks
  ks = bwa_open_reads(opt->mode, fn_fa);
  // core loop
  while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual)) != 0) {
    tot_seqs += n_seqs;
    t = clock();

    /**  ************************ **/
    // duplicate reads
    // take a copy of seqs
    seqsC = obwape_copy_seqs(seqs,n_seqs,opt);
    
    // do alignment
    fprintf(stderr, "[bwa_aln_core] calculate SA coordinate ... ");
    obwape_bwa_aln_core_inner(bwt, seqsC, n_seqs,opt);
    fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

    /*
    // read alignment
    for (i = 0; i < n_seqs; ++i) {
      bwa_seq_t *p = seqs + i;
      int n_aln;
      fread(&n_aln, 4, 1, fp_sa);
      if (n_aln > m_aln) {
	m_aln = n_aln;
	aln = (bwt_aln1_t*)realloc(aln, sizeof(bwt_aln1_t) * m_aln);
      }
      fread(aln, sizeof(bwt_aln1_t), n_aln, fp_sa);
      bwa_aln2seq_core(n_aln, aln, p, 1, n_occ);
    }
    */
    for (i = 0; i < n_seqs; ++i) {
      bwa_seq_t *p = seqs + i;
      bwa_seq_t *px = seqsC + i;
      int n_aln = px->n_aln;
      bwa_aln2seq_core(n_aln, px->aln, p, 1, n_occ);
    }

    fprintf(stderr, "[bwa_aln_core] convert to sequence coordinate... ");
    // bwa_cal_pac_pos(bns, prefix, n_seqs, seqs, opt.max_diff, opt.fnr); // forward bwt will be destroyed here
    obwase_bwa_cal_pac_pos(bns, bwt, n_seqs, seqs, opt->max_diff, opt->fnr); // forward bwt will be destroyed here                                              
    fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

    fprintf(stderr, "[bwa_aln_core] refine gapped alignments... ");
    bwa_refine_gapped(bns, n_seqs, seqs, 0, ntbns);
    fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

    fprintf(stderr, "[bwa_aln_core] print alignments... ");
    for (i = 0; i < n_seqs; ++i)
      bwa_print_sam1(bns, seqs + i, 0, opt->mode, opt->max_top2);
    fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

    bwa_free_read_seq(n_seqs, seqs);
    fprintf(stderr, "[bwa_aln_core] %d sequences have been processed.\n", tot_seqs);
  }

  // destroy
  bwa_seq_close(ks);
  if (ntbns) bns_destroy(ntbns);
  bns_destroy(bns);
  // fclose(fp_sa);
  free(aln);

  // additionally, free variables incorporated from bwa_aln_core
  bwt_destroy(bwt);

  return 0;
}

int parse_samse(gap_opt_t *opt, char *bwa_rg_line, char *bwa_rg_id, char *args, int *n_occ) {

  extern int bwa_set_rg(const char *s);
  // tokenise arguments                                                                                                                                               
  enum { kMaxArgs = 64 };
  int argc = 0;
  char *argv[kMaxArgs];
  char *p2 = strtok(args, " ");
  while (p2 && argc < kMaxArgs) {
    argv[argc++] = p2;
    p2 = strtok(0, " ");
  }
  
  int c;
  *n_occ = 3;

  while ((c = getopt(argc, argv, "hn:f:r:")) >= 0) {
    switch (c) {
    case 'h': break;
    case 'r':
      if (bwa_set_rg(optarg) < 0) {
	fprintf(stderr, "[%s] malformated @RG line\n", __func__);
	return 1;
      }
      break;
    case 'n': *n_occ = atoi(optarg); break;
    case 'f': xreopen(optarg, "w", stdout); break;
    default: return 1;
    }
  }

  return 0;
}

int print_bwase_usage(gap_opt_t *opt) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Program:  obwase combines \'aln\' and \'samse\' from bwa (version %s) into one command.\n",PACKAGE_VERSION);
  fprintf(stderr, "          If executed with <prefix> only, returns the sam header.\n");
  fprintf(stderr, "Version:  %s\n",BWASE_VERSION);
  fprintf(stderr, "Contact:  Margus Lukk <margus.lukk@cruk.cam.ac.uk>\n\n");
  fprintf(stderr, "Usage:    obwase [options] <prefix> [<in1.fq>]\n");
  fprintf(stderr, "Options:  aln      takes bwa aln optsions\n");
  fprintf(stderr, "          bwase    bwa bwase optsions\n");
  fprintf(stderr, "          -s       skip sam header in output sam file\n\n");
  fprintf(stderr, "Examples: obwase hg19.fa in1.fq\n");
  fprintf(stderr, "          obwase aln \"-e 2 -t 5\" bwase \"-n 30\" hg19.fa in1.fq\n\n");
  fprintf(stderr, "Options for \'aln\' (bwa version %s):\n\n",PACKAGE_VERSION);
  fprintf(stderr, "         -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n",
          BWA_AVG_ERR, opt->fnr);
  fprintf(stderr, "         -o INT    maximum number or fraction of gap opens [%d]\n", opt->max_gapo);
  fprintf(stderr, "         -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
  fprintf(stderr, "         -i INT    do not put an indel within INT bp towards the ends [%d]\n", opt->indel_end_skip);
  fprintf(stderr, "         -d INT    maximum occurrences for extending a long deletion [%d]\n", opt->max_del_occ);
  fprintf(stderr, "         -l INT    seed length [%d]\n", opt->seed_len);
  fprintf(stderr, "         -k INT    maximum differences in the seed [%d]\n", opt->max_seed_diff);
  fprintf(stderr, "         -m INT    maximum entries in the queue [%d]\n", opt->max_entries);
  fprintf(stderr, "         -t INT    number of threads [%d]\n", opt->n_threads);
  fprintf(stderr, "         -M INT    mismatch penalty [%d]\n", opt->s_mm);
  fprintf(stderr, "         -O INT    gap open penalty [%d]\n", opt->s_gapo);
  fprintf(stderr, "         -E INT    gap extension penalty [%d]\n", opt->s_gape);
  fprintf(stderr, "         -R INT    stop searching when there are >INT equally best hits [%d]\n", opt->max_top2);
  fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
  // fprintf(stderr, "         -f FILE   file to write output to instead of stdout\n"); // deprecated as sai files do not get written
  fprintf(stderr, "         -B INT    length of barcode\n");
  //              fprintf(stderr, "         -c        input sequences are in the color space\n");
  fprintf(stderr, "         -L        log-scaled gap penalty for long deletions\n");
  fprintf(stderr, "         -N        non-iterative mode: search for all n-difference hits (slooow)\n");
  fprintf(stderr, "         -I        the input is in the Illumina 1.3+ FASTQ-like format\n");
  fprintf(stderr, "         -b        the input read file is in the BAM format\n");
  fprintf(stderr, "         -0        use single-end reads only (effective with -b)\n");
  fprintf(stderr, "         -1        use the 1st read in a pair (effective with -b)\n");
  fprintf(stderr, "         -2        use the 2nd read in a pair (effective with -b)\n");
  fprintf(stderr, "         -Y        filter Casava-filtered sequences\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options for \'samse\' (bwa version %s):\n\n",PACKAGE_VERSION);
  fprintf(stderr, "         -n max_occ\n");
  fprintf(stderr, "         -f out.sam\n");
  fprintf(stderr, "         -r RG_line\n");
  fprintf(stderr, "\n");
  return 1;
}

int print_sam_header(char *prefix) {
  bntseq_t *bns = bns_restore(prefix);

  bwa_print_sam_SQ(bns);
  bwa_print_sam_PG();

  bns_destroy(bns);

  return 0;
}

int obwase(int argc, char *argv[]) {
  int i, j, skipheader = 0;
  // vars for aln

  gap_opt_t *opt;

  // vars for samse
  int n_occ = 3;
  
  extern char *bwa_rg_line, *bwa_rg_id;
  
  // install opt
  opt = gap_init_opt();
  // install opt and popt
  if (opt->fnr > 0.0) {
    int i, k;
    for (i = 17, k = 0; i <= 250; ++i) {
      int l = bwa_cal_maxdiff(i, BWA_AVG_ERR, opt->fnr);
      if (l != k) { 
	// fprintf(stderr, "[bwa_aln] %dbp reads: max_diff = %d\n", i, l); 
      }
      k = l;
    }
  }

  // parse options

  int c;

  while ((c = getopt(argc, argv, "s")) >= 0) {
    switch (c) {
    case 's':
      skipheader = 1;
      break;
    default: return 1;
    }
  }

  j=optind;
  for(i = optind; i < argc; i++ ) {
    if(strcmp(argv[i],"aln") == 0) {
      if( (i+1) < argc) {
        j=i+2;
        i++;
        parse_aln(opt,argv[i]);
      }
      else print_bwase_usage(opt);
      continue;
    }
    if(strcmp(argv[i],"samse") == 0) {
      if( (i+1) < argc) {
        j=i+2;
        i++;
        parse_samse(opt,bwa_rg_line,bwa_rg_id,argv[i],&n_occ);
      }
      else print_bwase_usage(opt);
      continue;
    }
  }  

  if( (j+1 > argc) || (j+2 < argc) ) {
    print_bwase_usage(opt);
  }

  if(j+1 == argc) {
    if(!skipheader) {
      print_sam_header(argv[j]);
    }
    else {
      print_bwase_usage(opt);
    }
  }
  
  if(j+2 == argc) {
    // aln + sampe
    obwase_core(argv[j],argv[j+1],opt,n_occ,skipheader);
  }
  
  // free vars for sampe

  free(opt);
  free(bwa_rg_line); free(bwa_rg_id);

  return 0;
}
