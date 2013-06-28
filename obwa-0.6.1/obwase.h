/*
  This work is published under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007.
  
  Copyright (c) 2013 by Margus Lukk <margus.lukk@cruk.cam.ac.uk>
  
  This work uses code from bwa-0.6.1 published under the same licence by Heng Li <lh3lh3@gmail.com>
*/
#include "bwtaln.h" // needed for gap_opt_t

#ifndef OBWASE_H
#define OBWASE_H

#ifdef __cplusplus
extern "C" {
#endif

  int obwase(int argc, char *argv[]);
  
#ifdef __cplusplus
}
#endif

#endif
