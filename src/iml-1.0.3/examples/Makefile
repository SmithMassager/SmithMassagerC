.PHONY: all clean

TGTS := #empty
TGTS += exam-certsolve
TGTS += exam-nonsingsolv
TGTS += exam-nullspace


all: $(TGTS)

%: %.c
	$(CC) $^ -L../lib -liml -lcblas -lopenblas -lgmp -o $@

clean:
	rm $(TGTS)
