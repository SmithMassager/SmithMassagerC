.PHONY: all clean

TGTS := #empty
TGTS += test-largeentry
TGTS += test-smallentry


all: $(TGTS)

%: %.c
	$(CC) $^ -L../lib -liml -lcblas -lopenblas -lgmp -o $@

clean:
	rm $(TGTS)
