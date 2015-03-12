default:
	nvcc -gencode arch=compute_35,code=sm_35 -ccbin gcc -c cuda_api.cu
	ftn -eZ -hsystem_alloc cuda_api.o -c gpu_interfaces.F90 
	ftn -eZ -hsystem_alloc -c RIMP2.GPU.F90
#	ftn -eZ -hsystem_alloc -c RIMP2.F90
#	ftn -eZ -hsystem_alloc -c generate.F90
	ftn -eZ -hsystem_alloc -o RIMP2.GPU.async.x cuda_api.o gpu_interfaces.o RIMP2.GPU.o
#	ftn -eZ -hsystem_alloc -o RIMP2.x RIMP2.o
#	ftn -eZ -hsystem_alloc -o generate.x generate.o
clean:
	rm *.o
	rm *.mod
	rm *.i
#	rm RIMP2.x
	rm RIMP2.GPU.x
#	rm generate.x
