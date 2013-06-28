/*
  This work is published under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007.

  Copyright (c) 2013 by Margus Lukk <margus.lukk@cruk.cam.ac.uk>

  This work uses code from bwa-0.6.1 published under the same licence by Heng Li <lh3lh3@gmail.com>
*/
#include <stdio.h>
#include <string.h>
#include "main.h"
#include "utils.h"
#include "obwase.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.6.1-r104"
#endif

void bwa_print_sam_PG()
{
	printf("@PG\tID:bwa\tPN:bwa\tVN:%s\n", PACKAGE_VERSION);
}

int main(int argc, char *argv[]) {
  return obwase(argc,argv);
}
