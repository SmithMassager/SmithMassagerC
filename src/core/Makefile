
.PHONY: core

#OBJDIR  := objs/
OBJS    := #empty
OBJS    += $(OBJDIR)arith_utils.o
OBJS    += $(OBJDIR)openblas_residue.o
OBJS    += $(OBJDIR)basis.o
OBJS    += $(OBJDIR)dbg_print.o
OBJS    += $(OBJDIR)mpz_matrix.o
OBJS    += $(OBJDIR)reconstruct.o
OBJS    += $(OBJDIR)rns_matrix.o
OBJS    += $(OBJDIR)thresholds.o
OBJS    += $(OBJDIR)timer.o
OBJS    += $(OBJDIR)pl_matrix.o
OBJS    += $(OBJDIR)basic.o
OBJS    += $(OBJDIR)projection.o
OBJS    += $(OBJDIR)fastIntCert.o
OBJS    += $(OBJDIR)specialIntCert.o
OBJS    += $(OBJDIR)fmpz_conv.o
OBJS    += $(OBJDIR)fmpz_mat_helper.o
OBJS    += $(OBJDIR)indexMassager.o
OBJS    += $(OBJDIR)largestInvariantFactor.o
OBJS    += $(OBJDIR)smithMassager.o
OBJS    += $(OBJDIR)genHerm.o
OBJS    += $(OBJDIR)main.o

core: $(OBJS)

$(OBJDIR)%.o : %.c
	$(CC) $< $(CFLAGS) -c -o $@

clean:
	rm $(OBJS)
