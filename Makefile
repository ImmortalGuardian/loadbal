TARGET = loadbal

CC = mpicc
CFLAGS = -Wall -O2
CFLAGS += -I ./include

SRCS = $(wildcard *.c)
#INCLUDE = $(wildcard include/*.h)

all : $(SRCS) .depend
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET)

.depend :
	$(CC) $(CFLAGS) -M $(SRCS) > .depend

clean :
	rm -f *.o $(TARGET)
	rm -f .depend
	rm -f .*.swp
	rm -f include/.*.swp
