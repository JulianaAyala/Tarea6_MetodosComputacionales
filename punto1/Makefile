all: plot_solutions.py
		python plot_solutions.py

datosxy: volterra_lotka.x
		./volterra_lotka.x

volterra_lotka: volterra_lotka.c
		cc volterra_lotka.c -o volterra_lotka.x -lm

clean:
	rm *.x *.dat *.png