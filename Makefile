TARGET = loadbal

CC = mpicc
CFLAGS = -Wall -O2
CFLAGS += -I ./include
LDFLAGS = -lm

SRCS = $(wildcard *.c)
#INCLUDE = $(wildcard include/*.h)

all : $(SRCS) .depend
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

.depend :
	$(CC) $(CFLAGS) -M $(SRCS) > .depend

clean :
	rm -f *.o $(TARGET)
	rm -f .depend
	rm -f .*.swp
	rm -f include/.*.swp

countlines:
	@find . -name '*.c' -o -name '*.h' | xargs wc -l |\
		gawk -e 'END {print $$1}'
