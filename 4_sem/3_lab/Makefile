mmake:
	gcc -Wall -fopenmp -o m1 midpoint.c -lm
mmkmake:
	gcc -Wall -fopenmp -o m2 mmk.c -lm

rmid: mmake
	./m1
	make clean
rmmk: mmkmake
	./m2
	make clean
graf:
	gnuplot graf0.gnu
	gnuplot graf1.gnu
clean:
	rm m1 m2