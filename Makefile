OPTF90 = -g 
OPTLD  = -g

F90    = gfortran

STOR   = csr_stor_class#ell_stor_class

all:    main.o contxt_class.o $(STOR).o \
		blas_class.o domain_class.o matrix_class.o
	$(F90) $(OPTLD) -o ns_fd.out main.o \
		contxt_class.o domain_class.o \
		$(STOR).o \
		blas_class.o matrix_class.o function_class.o

clean:
	rm -f *.o *.mod ns_fd.out

main.o: main.f90 contxt_class.mod blas_class.mod $(STOR).mod \
		domain_class.mod matrix_class.mod 
	$(F90) $(OPTF90) -c main.f90

contxt_class.o contxt_class.mod: contxt_class.f90
	$(F90) $(OPTF90) -c contxt_class.f90

funtion_class.o function_class.mod: function_class.f90
	$(F90) $(OPTF90) -c function_class.f90

ell_stor_class.o ell_stor_class.mod: ell_stor_class.f90 contxt_class.mod
	$(F90) $(OPTF90) -c ell_stor_class.f90

csr_stor_class.o csr_stor_class.mod: csr_stor_class.f90 contxt_class.mod
	$(F90) $(OPTF90) -c csr_stor_class.f90

#my_stor_class.o my_stor_class.mod: my_stor_class.f90 contxt_class.mod
#	$(F90) $(OPTF90) -c my_stor_class.f90

blas_class.o blas_class.mod: blas_class.f90 contxt_class.mod $(STOR).mod
	$(F90) $(OPTF90) -c blas_class.f90

matrix_class.o matrix_class.mod: matrix_class.f90 contxt_class.mod \
		$(STOR).mod domain_class.mod function_class.mod
	$(F90) $(OPTF90) -c matrix_class.f90

domain_class.o domain_class.mod: domain_class.f90 contxt_class.mod 
	$(F90) $(OPTF90) -c domain_class.f90

