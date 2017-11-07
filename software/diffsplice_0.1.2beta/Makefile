all: DiffSplice SeparateSAM SeparateFragment GenomeTree DiffExprAnalysis DiffASMAnalysis

OPTFLAGS = -m64


#LDFLAGS = -static

CFLAGS += $(OPTFLAGS)

#CFLAGS += $(LDFLAGS)

DiffSplice:
	g++ $(CFLAGS) -o diffsplice src/diffsplice/diffsplice.cpp

SeparateSAM:
	g++ $(CFLAGS) -o bin/sepSAM src/sepSAM.cpp 

SeparateFragment:
	g++ $(CFLAGS) -o bin/fragment src/parseFragment.cpp 

GenomeTree:
	g++ $(CFLAGS) -o bin/tree src/gtree/GTree.cpp 

DiffExprAnalysis:
	g++ $(CFLAGS) -o bin/diff_expr_analysis src/testexpr/ExpressionAnalysis.cpp src/testexpr/selectSignificance.cpp 

DiffASMAnalysis:
	g++ $(CFLAGS) -o bin/diff_asm_analysis src/testtrans/TranscriptionAnalysis.cpp src/testtrans/selectSignificance_asm.cpp 


