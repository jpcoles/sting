CC=gcc
CFLAGS=-Wall -O2 -g
LDFLAGS=

all: sting

sting: sting.o ic.o fg.o log.o
	$(CC) $(LDFLAGS) $^ -o $@

sting.o: sting.c sting.h
	$(CC) $(CFLAGS) -c $< -o $@

ic.o: ic.c ic.h
	$(CC) $(CFLAGS) -c $< -o $@

fg.o: fg.c fg.h
	$(CC) $(CFLAGS) -c $< -o $@

log.o: log.c log.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f *.o sting
