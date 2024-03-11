CFLAGS += -O2 -I../../src/ -I$(BASILISK)/../thirdParty/cJSON-1.7.15
LDFLAGS += -lm -lcjson -L$(BASILISK)/../thirdParty/cJSON-1.7.15/build -Wl,-rpath $(BASILISK)/../thirdParty/cJSON-1.7.15/build
SOURCE := bin2ascii.c
ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
INST = $(shell sudo ln -sf $(ROOT_DIR)/bin/* /usr/local/bin/)
UNINST = $(shell ls $(ROOT_DIR)/bin/* | xargs -n 1 find -L /usr/local/bin -samefile | xargs sudo rm -f )

all:
	make wclean
	make bin2ascii ascii2bin
	make wclean
	rm -rf *.dSYM

bin2ascii: bin2ascii.c
	$(QCC) $(CFLAGS) -D_MPI=0 -o $@ $^ $(LDFLAGS)

ascii2bin: ascii2bin.c
	$(QCC) $(CFLAGS) -D_MPI=0 -o $@ $^ $(LDFLAGS)

winstall:
	echo $(ROOT_DIR)
	echo $(INST)

wuninstall:
	echo $(ROOT_DIR)
	echo $(UNINST)

wclean:
	rm -f *.exec _*.c rm -rf .qcc*
	make clean
	
include $(BASILISK)/Makefile.defs
