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
	@rm tsp -f $(OBJDIR)/*.o $(OBJDIR)/*.d

runall:
	./tsp instances/a280.tsp instances/kroA150.tsp instances/ali535.tsp instances/kroA200.tsp instances/att48.tsp instances/kroB100.tsp instances/att532.tsp instances/kroB150.tsp instances/bayg29.tsp instances/kroB200.tsp instances/bays29.tsp instances/kroC100.tsp instances/berlin52.tsp instances/kroD100.tsp instances/bier127.tsp instances/kroE100.tsp instances/brazil58.tsp instances/lin105.tsp instances/brg180.tsp instances/lin318.tsp instances/burma14.tsp instances/linhp318.tsp instances/ch130.tsp instances/pcb442.tsp instances/ch150.tsp instances/pr107.tsp instances/d198.tsp instances/pr124.tsp instances/d493.tsp instances/pr136.tsp instances/dantzig42.tsp instances/pr144.tsp instances/eil101.tsp instances/pr152.tsp instances/eil51.tsp instances/pr226.tsp instances/eil76.tsp instances/pr264.tsp instances/fl417.tsp instances/pr299.tsp instances/fri26.tsp instances/pr76.tsp instances/gil262.tsp instances/rat195.tsp instances/gr120.tsp instances/rat99.tsp instances/gr137.tsp instances/rd100.tsp instances/gr17.tsp instances/rd400.tsp instances/gr202.tsp instances/si175.tsp instances/gr21.tsp instances/si535.tsp instances/gr229.tsp instances/st70.tsp instances/gr24.tsp instances/swiss42.tsp instances/gr431.tsp instances/ts225.tsp instances/gr48.tsp instances/tsp225.tsp instances/gr96.tsp instances/u159.tsp instances/hk48.tsp instances/ulysses16.tsp instances/kroA100.tsp instances/ulysses22.tsp 


rebuild: clean tsp
