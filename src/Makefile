#
#=======================================#
#   Makefile to build  m c c i   2 . 0  #
#=======================================#
#
       DEST   = ./
#
#------PGF90 on walton.ichec.ie------------------
# 
#      FC     = mpif90
#      FFLAGS = -tp=k8-64 -i4 -r8 -fpic -fastsse -O4
#      LFLAGS = -tp=k8-64 -fpic -fastsse -O4 
#      LFLAGS = -tp=k8-64 -fpic -L/usr/pgi/linux86-64/6.0/libso -fastsse -O4 
#      FFLAGS = -tp=k8-64 -i4 -r8 -fpic -mcmodel=medium -mp -fastsse -O4 -Minline=levels:8 -I/opt/packages/gcc-compat/acml/include 
#      LFLAGS = -tp=k8-64 -mp -fastsse -O4 -Minline=levels:8 -L/opt/packages/pgi-compat/acml/lib -fpic -mcmodel=medium -lblacs -lacml -lscalapack
# 
#      LIBS   = ./libacml.a 
# 
#------F90 on sfi500.tyndall.ie----------------
#   
#        FC     = /opt/mpich2/gnu/bin/mpif90
#        FFLAGS =  -i4 -r8 -O4 
#        LFLAGS = -llapack -lblas
#        LIBS   =
#
#------F90 on rocks-hpc-------------------------

         FC     =  mpif90
         FFLAGS = -i4 -r8 -O4 
         LFLAGS = -llapack -lblas
         LIBS   =  

  
#------PGF90 on sfi400.tyndall.ie----------------
#  
#       FC     = /opt/mpich-1.2.7/test_pgi6_64/bin/mpif90
#       FFLAGS = -tp=k8-64 -i4 -r8 -fastsse -O4 -fpic 
#       LFLAGS = -tp=k8-64 -lacml -fastsse -O4 -fpic -L/usr/pgi/linux86-64/6.0/libso 
#       FFLAGS = -tp=k8-64 -i4 -r8 -fastsse -O4 -mcmodel=medium
#       LFLAGS = -tp=k8-64 -lacml -fastsse -O4 -mcmodel=medium
#       LIBS   =
# 
#------PGF90 on sfi300.tyndall.ie----------------
# 
#      FC     = mpif90
#      FFLAGS = -i4 -r8 -fastsse -O4 -fpic 
#      LFLAGS = -llapack -lblas -fastsse -O4 -fpic -L/usr/local/pgi/linux86/5.2/lib
#      LIBS   =
# 
#------PATHF90 on sfi400.tyndall.ie--------------
#
#        FC     = mpif90
#        FFLAGS = -m64 -shlib -i4 -r8 -O3 -OPT:Ofast
#        LFLAGS = -m64 -shlib -L/opt/acml/pathscale64/lib -lacml
#        LIBS   =
#
#------F95 on UNIX (Sun)--------------------------
#
#      FC     = f90
#      FFLAGS = -fast -xarch=v9
#      LFLAGS =
#      LAPACK = lapack_sun_f90.a
#      LIBS   = blas_sun_f90.a
#
#-------------------------------------------------

OBJS =  mcci.o get_int.o branch.o singles.o doubles.o reorder.o swap.o\
        h_s_sparse.o random.o chk_list.o sym_init.o prune.o init.o\
        gen_seed.o energy.o sch0.o sch1.o sch1_0.o sch1_1.o sch1_2.o\
        sch2.o sch2_0.o sch2_1.o sch2_2.o sch2_3.o sch2_4.o sch2_5.o\
        get_case1.o get_case2.o restart.o muHnu.o mxv_sparse.o \
        wxmxv_sparse.o h_move.o s_move.o timer.o ck.o dump.o   \
        ldump.o genealogy.o davidson.o h_s_reduced.o init_bk.o exc.o \
        get_int_TM.o get_int_bTM.o precision.o read_mcci_in.o commonarrays.o \
	dyn_par.o read_params.o

mcci:   $(OBJS)
	$(FC) $(LFLAGS) -o $(DEST)$@ $(OBJS) $(LAPACK) $(LIBS) 

clean:
	rm -f *.o *.mod mcci
 
mcci.o:			dyn_par.o commonarrays.o read_params.o
get_int.o:              dyn_par.o
h_s_sparse.o:           dyn_par.o
h_s_reduced.o:          
muHnu.o:		dyn_par.o
reorder.o:		dyn_par.o
swap.o:   		
mxv_sparse.o:		
wxmxv_sparse.o:		
h_move.o:		
s_move.o:		
branch.o:		dyn_par.o
genealogy.o:		dyn_par.o
prune.o: 		dyn_par.o
singles.o:      	dyn_par.o
doubles.o:	        dyn_par.o
chk_list.o:		
sym_init.o:		
cmminit.o:		
energy.o:  		
sch0.o:			
sch1.o:			
sch1_0.o:		
sch1_1.o:		
sch1_2.o:		
sch2.o:	   		
sch2_0.o:		
sch2_1.o:		
sch2_2.o:		
sch2_3.o:		
sch2_4.o:		
sch2_5.o:		
get_case1.o:		
get_case2.o:		
restart.o:		
dump.o:  		
ldump.o:  		
davidson.o:		dyn_par.o
h_s_reduced.o:		
init_bk.o:		
get_int_TM.o:		dyn_par.o
get_int_bTM.o:		dyn_par.o
precision.o:            
exc.o:                  dyn_par.o
read_mcci_in.o:		precision.o dyn_par.o
commonarrays.o:         
dyn_par.o:
read_params.o:          precision.o dyn_par.o commonarrays.o

.SUFFIXES:
.SUFFIXES: .o .f .f90

.f.o:
	$(FC) $(FFLAGS) -c $*.f

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90



