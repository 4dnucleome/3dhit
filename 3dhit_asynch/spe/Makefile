include ../make.debug

PROGRAMS_spu = worker_spe
SPU_COMPILER = gcc
LIBRARY_embed = worker_spe.a
CPPFLAGS = $(DEBFLAGS)
CFLAGS = -ffast-math
CC_OPT_LEVEL = -O3 -ftree-vectorize
IMPORTS = -lm

 #-fdump-tree-vect-details

include $(CELL_TOP)/buildutils/make.footer
