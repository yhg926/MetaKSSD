CC=gcc
CFLAGS= -std=gnu11 -Wall -Wno-format-overflow -Wno-unused-result -O3 -ggdb -lz -fopenmp
SOURCE="./src"
BIN="./bin"
PRONAME="metakssd"

all:
	$(CC) $(CFLAGS)  $(SOURCE)/*.c -o $(BIN)/$(PRONAME) -lm
alert:
	$(CC) $(CFLAGS) -DCOMPONENT_SZ=8 $(SOURCE)/*.c -o $(BIN)/$(PRONAME)_CSZ8 -lm
strange:
	$(CC) $(CFLAGS) -DCTX_SPC_USE_L=10 $(SOURCE)/*.c -o $(BIN)/$(PRONAME)_strange -lm
16S:
	$(CC) $(CFLAGS) -DMIN_KM_S=1 $(SOURCE)/*.c -o $(BIN)/$(PRONAME)_16S -lm

