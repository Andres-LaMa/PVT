mmake:
	gcc -Wall -fopenmp -o m1 midpoint.c -lm
mmkmake:
	gcc -Wall -fopenmp -o m2 mmk.c -lm

rmid: mmake
	./m1
rmmk: mmkmake
	./m2
graf:
	gnuplot graf.gnu
clean:
	rm m1 m2