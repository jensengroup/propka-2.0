FC = gfortran

FCFLAGS = -O2

propka2.0: propka2.0_2008-11-12.o
	${FC} ${FCFLAGS} $^ -o $@

propka2.0_2008-11-12.o: propka2.0_2008-11-12.f
	${FC} ${FCFLAGS} -c $^

clean:
	rm -f *.o
