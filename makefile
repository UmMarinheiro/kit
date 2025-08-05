#detecta se o sistema é de 32 ou 64 bits
BITS_OPTION = -m64

#### define o compilador
CPPC = g++
#############################

#### opcoes de compilacao e includes
CCOPT = $(BITS_OPTION) -O3 -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++0x
CONCERTINCDIR = $(CONCERTDIR)/include
CCFLAGS = $(CCOPT)
#############################

#### flags do linker
#############################

#### diretorios com os source files e com os objs files
SRCDIR = src
OBJDIR = obj
TESTDIR = tests
#############################

#### lista de todos os srcs e todos os objs
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))
#############################

#### regra principal, gera o executavel
tsp: $(OBJS) 
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(BITS_OPTION) $(OBJS) -o $@ $(CCLNFLAGS)
############################

#inclui os arquivos de dependencias
-include $(OBJS:.o=.d)

#regra para cada arquivo objeto: compila e gera o arquivo de dependencias do arquivo objeto
#cada arquivo objeto depende do .c e dos headers (informacao dos header esta no arquivo de dependencias gerado pelo compiler)
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CPPC) $(CCFLAGS) -c $< -o $@
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CPPC) -std=c++0x  -MM $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp #proximas tres linhas colocam o diretorio no arquivo de dependencias (g++ nao coloca, surprisingly!)
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

#delete objetos e arquivos de dependencia
clean:
	@echo "\033[31mcleaning obj directory \033[0m"
	@rm tsp -f $(OBJDIR)/*.o $(OBJDIR)/*.d $(TESTDIR)/*.txt

runall: tsp
	./tsp instances/dantzig42.tsp > ${TESTDIR}/dantzig42.txt;
	./tsp instances/swiss42.tsp > ${TESTDIR}/swiss42.txt;
	./tsp instances/att48.tsp > ${TESTDIR}/att48.txt;
	./tsp instances/gr48.tsp > ${TESTDIR}/gr48.txt;
	./tsp instances/hk48.tsp > ${TESTDIR}/hk48.txt;
	./tsp instances/eil51.tsp > ${TESTDIR}/eil51.txt;
	./tsp instances/berlin52.tsp > ${TESTDIR}/berlin52.txt;
	./tsp instances/brazil58.tsp > ${TESTDIR}/brazil58.txt;
	./tsp instances/st70.tsp > ${TESTDIR}/st70.txt;
	./tsp instances/eil76.tsp > ${TESTDIR}/eil76.txt;
	./tsp instances/pr76.tsp > ${TESTDIR}/pr76.txt;
	./tsp instances/gr96.tsp > ${TESTDIR}/gr96.txt;
	./tsp instances/rat99.tsp > ${TESTDIR}/rat99.txt;
	./tsp instances/kroA100.tsp > ${TESTDIR}/kroA100.txt;
	./tsp instances/kroB100.tsp > ${TESTDIR}/kroB100.txt;
	./tsp instances/kroC100.tsp > ${TESTDIR}/kroC100.txt;
	./tsp instances/kroD100.tsp > ${TESTDIR}/kroD100.txt;
	./tsp instances/kroE100.tsp > ${TESTDIR}/kroE100.txt;
	./tsp instances/rd100.tsp > ${TESTDIR}/rd100.txt;
	./tsp instances/eil101.tsp > ${TESTDIR}/eil101.txt;
	./tsp instances/lin105.tsp > ${TESTDIR}/lin105.txt;
	./tsp instances/pr107.tsp > ${TESTDIR}/pr107.txt;

rebuild: clean tsp
