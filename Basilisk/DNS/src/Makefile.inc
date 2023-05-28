CC := $(if $(filter $(COMPILEMPI),1),mpicc,gcc)

.PHONY: clean all cleansim

DEPS := $(shell find $(BASILISK_USER_SOURCE) -maxdepth 1 -name "*.h" -printf "%f\n")

# Uncomment to not remove the header files after compilation:
#.SECONDARY: $(DEPS)

all: $(SOURCE)

$(SOURCE): _$(SOURCE).c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

_%.c: %.c $(DEPS)
	qcc -source -grid=$(GRID) -D_MPI=$(COMPILEMPI) $<

%.h: $(BASILISK_USER_SOURCE)/%.h
	ln -s $< $@

clean: cleansim
	rm -f $(SOURCE) _$(SOURCE).c

cleansim:
	rm -fr [0-9][0-9][0-9][0-9][0-9][0-9] log* restart out-* *.pvd .qcc* *.h $(BASILISK_CLEAN_APP)
