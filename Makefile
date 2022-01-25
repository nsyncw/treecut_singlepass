CCOMP = gcc
#CFLAGS = -O4 -DNDEBUG -DEXCESS_TYPE_LONG -DPRINT_STAT -DCHECK_SOLUTION -Wall -lm
CFLAGS = -g -DPRINT_FLOW -DEXCESS_TYPE_LONG -DPRINT_STAT -DCHECK_SOLUTION -Wall -lm

all: hi_treem 
hi_treem: main.c hi_treem.c parser_treem.c timer.c
	$(CCOMP) $(CFLAGS) -o main_treem main.c libm.so 
clean: 
	rm -f main_treem *~
