include make.debug

DIRS = spe

PROGRAM_ppu = 3dhit
PPU_COMPILER = gcc
CPPFLAGS = $(DEBFLAGS)
CFLAGS = -ffast-math
CC_OPT_LEVEL = -O3
VPATH = modules:.
IMPORTS = spe/worker_spe.a -lm -lmisc -lspe2 -lpthread 

include $(CELL_TOP)/buildutils/make.footer
