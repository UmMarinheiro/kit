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

#### quantas vezes cada teste vai
FREQUENCY = 10

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
	./tsp instances/a280.tsp > $(TESTDIR)/a280.txt;  	
	./tsp instances/kroA150.tsp > $(TESTDIR)/kroA150.txt;  
	./tsp instances/kroA200.tsp > $(TESTDIR)/kroA200.txt;  
	./tsp instances/att48.tsp > $(TESTDIR)/att48.txt;  
	./tsp instances/kroB100.tsp > $(TESTDIR)/kroB100.txt;  
	./tsp instances/kroB150.tsp > $(TESTDIR)/kroB150.txt; 
	./tsp instances/bayg29.tsp > $(TESTDIR)/bayg29.txt; 
	./tsp instances/kroB200.tsp > $(TESTDIR)/kroB200.txt; 
	./tsp instances/bays29.tsp > $(TESTDIR)/bays29.txt; 
	./tsp instances/kroC100.tsp > $(TESTDIR)/kroC100.txt; 
	./tsp instances/berlin52.tsp > $(TESTDIR)/berlin52.txt; 
	./tsp instances/kroD100.tsp > $(TESTDIR)/kroD100.txt; 
	./tsp instances/bier127.tsp > $(TESTDIR)/bier127.txt; 
	./tsp instances/kroE100.tsp > $(TESTDIR)/kroE100.txt; 
	./tsp instances/brazil58.tsp > $(TESTDIR)/brazil58.txt; 
	./tsp instances/lin105.tsp > $(TESTDIR)/lin105.txt; 
	./tsp instances/brg180.tsp > $(TESTDIR)/brg180.txt; 
	./tsp instances/burma14.tsp > $(TESTDIR)/burma14.txt; 
	./tsp instances/ch130.tsp > $(TESTDIR)/ch130.txt; 
	./tsp instances/pcb442.tsp > $(TESTDIR)/pcb442.txt; 
	./tsp instances/ch150.tsp > $(TESTDIR)/ch150.txt; 
	./tsp instances/pr107.tsp > $(TESTDIR)/pr107.txt; 
	./tsp instances/d198.tsp > $(TESTDIR)/d198.txt; 
	./tsp instances/pr124.tsp > $(TESTDIR)/pr124.txt; 
	./tsp instances/pr136.tsp > $(TESTDIR)/pr136.txt; 
	./tsp instances/dantzig42.tsp > $(TESTDIR)/dantzig42.txt; 
	./tsp instances/pr144.tsp > $(TESTDIR)/pr144.txt; 
	./tsp instances/eil101.tsp > $(TESTDIR)/eil101.txt; 
	./tsp instances/pr152.tsp > $(TESTDIR)/pr152.txt; 
	./tsp instances/eil51.tsp > $(TESTDIR)/eil51.txt; 
	./tsp instances/pr226.tsp > $(TESTDIR)/pr226.txt; 
	./tsp instances/eil76.tsp > $(TESTDIR)/eil76.txt; 
	./tsp instances/pr264.tsp > $(TESTDIR)/pr264.txt; 
	./tsp instances/pr299.tsp > $(TESTDIR)/pr299.txt; 
	./tsp instances/fri26.tsp > $(TESTDIR)/fri26.txt; 
	./tsp instances/pr76.tsp > $(TESTDIR)/pr76.txt; 
	./tsp instances/gil262.tsp > $(TESTDIR)/gil262.txt; 
	./tsp instances/rat195.tsp > $(TESTDIR)/rat195.txt; 
	./tsp instances/gr120.tsp > $(TESTDIR)/gr120.txt; 
	./tsp instances/rat99.tsp > $(TESTDIR)/rat99.txt; 
	./tsp instances/gr137.tsp > $(TESTDIR)/gr137.txt; 
	./tsp instances/rd100.tsp > $(TESTDIR)/rd100.txt; 
	./tsp instances/gr17.tsp > $(TESTDIR)/gr17.txt; 
	./tsp instances/gr202.tsp > $(TESTDIR)/gr202.txt; 
	./tsp instances/si175.tsp > $(TESTDIR)/si175.txt; 
	./tsp instances/gr21.tsp > $(TESTDIR)/gr21.txt; 
	./tsp instances/gr229.tsp > $(TESTDIR)/gr229.txt; 
	./tsp instances/st70.tsp > $(TESTDIR)/st70.txt; 
	./tsp instances/gr24.tsp > $(TESTDIR)/gr24.txt; 
	./tsp instances/swiss42.tsp > $(TESTDIR)/swiss42.txt; 
	./tsp instances/ts225.tsp > $(TESTDIR)/ts225.txt; 
	./tsp instances/gr48.tsp > $(TESTDIR)/gr48.txt; 
	./tsp instances/tsp225.tsp > $(TESTDIR)/tsp225.txt; 
	./tsp instances/gr96.tsp > $(TESTDIR)/gr96.txt; 
	./tsp instances/u159.tsp > $(TESTDIR)/u159.txt; 
	./tsp instances/hk48.tsp > $(TESTDIR)/hk48.txt; 
	./tsp instances/ulysses16.tsp > $(TESTDIR)/ulysses16.txt; 
	./tsp instances/kroA100.tsp > $(TESTDIR)/kroA100.txt; 
	./tsp instances/ulysses22.tsp > $(TESTDIR)/ulysses22.txt; 


rebuild: clean tsp
