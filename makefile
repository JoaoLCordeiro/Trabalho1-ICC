CFLAGS  = -Wall -g 
LDFLAGS = -lmatheval -lm
OBJS = main.o utils.o libGeral.o libNP.o libNM.o libNI.o libDefine.h

all: newtonPC

newtonPC: $(OBJS)
	gcc $(CFLAGS) $(OBJS) $(LDFLAGS)

main.o: main.c
	gcc $(CFLAGS) -c main.c

utils.o: utils.c
	gcc $(CFLAGS) -c utils.c

libGeral.o: libGeral.c
	gcc $(CFLAGS) -c libGeral.c

libNP.o: libNP.c
	gcc $(CFLAGS) -c libNP.c

libNM.o: libNM.c
	gcc $(CFLAGS) -c libNM.c

libNI.o: libNI.c
	gcc $(CFLAGS) -c libNI.c

clean:
	-rm -f *~ *.o