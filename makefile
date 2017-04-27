zedov.pdf : rad_dat.txt dens_dat.txt plot.py
	python plot.py

rad_dat.txt dens_dat.txt : a.out
	./a.out

a.out : Zedov.c
	gcc -lm Zedov.c

clean:
	rm -rf a.out rad_dat.txt dens_dat.txt Densidad.pdf
