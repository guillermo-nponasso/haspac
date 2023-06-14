CLIQUER_DIR = ~/libraries/cliquer

# Cliquer source files
CLSRC := $(CLIQUER_DIR)/cliquer.c $(CLIQUER_DIR)/graph.c $(CLIQUER_DIR)/reorder.c

# Cliquer objects
CLOBS := $(CLIQUER_DIR)/cliquer.o $(CLIQUER_DIR)/graph.o $(CLIQUER_DIR)/reorder.o


CC := gcc
CFLAGS := -Wall -I $(CLIQUER_DIR)

haspac: $(CLOBS)
	$(CC) $(CFLAGS) $(CURDIR)/haspac.c $(CLSRC) -o build/$@
