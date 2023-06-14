CLIQUER_DIR = ~/libraries/cliquer

# Cliquer source files
CLSRC := $(CLIQUER_DIR)/cliquer.c $(CLIQUER_DIR)/graph.c $(CLIQUER_DIR)/reorder.c

# Cliquer objects
CLOBS := $(CLIQUER_DIR)/cliquer.o $(CLIQUER_DIR)/graph.o $(CLIQUER_DIR)/reorder.o

CURDIR_SRC := $(CURDIR)/haspac.c $(CURDIR)/hsp_utils.c


CC := gcc
CFLAGS := -Wall -I $(CLIQUER_DIR)

haspac: $(CLOBS)
	$(CC) $(CFLAGS) $(CURDIR_SRC) $(CLSRC) -o build/$@

test : $(CLOBS)
	$(CC) $(CFLAGS) $(CURDIR)/test.c $(CURDIR)/hsp_utils.c $(CLSRC) -o build/$@
