F77 = gfortran
aggdriv.exe: aggdriv.o section.o secstrt.h
	$(F77) aggdriv.o section.o -o aggdriv.exe
aggdriv.o: aggdriv.f secstrt.h
	$(F77) aggdriv.f -c -o aggdriv.o
section.o: section.f secstrt.h
	$(F77) section.f -c -o section.o
clean:
	rm -r *.o *.exe *.out
run: output.out
output.out: aggdriv.exe
	./aggdriv.exe>output.out
