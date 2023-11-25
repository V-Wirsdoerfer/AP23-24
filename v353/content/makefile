auswertung.tex: build/auswertung.pdf


3_auswertung.short: auswertung.short
3_auswertung.long: auswertung.long
auswertung.short: FORCE | build
	lualatex --output-directory=build --interaction=nonstopmode --halt-on-error 3_auswertung.tex
auswertung.long: build/auswertung.long
build/auswertung.long: FORCE |build
	TEXINPUTS=../..: \
	BIBINPUTS=../..: \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	3_auswertung.tex


4_diskussion.short: diskussion.short
4_diskussion.long: diskussion.long
diskussion.short: FORCE | build
	lualatex --output-directory=build --interaction=nonstopmode --halt-on-error 4_diskussion.tex
diskussion.long: build/diskussion.long
build/diskussion.long: FORCE |build
	TEXINPUTS=../..: \
	BIBINPUTS=../..: \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	4_diskussion.tex

2_durchfuehrung.short: durchfuehrung.short
2_durchfuehrung.long: durchfuehrung.long
durchfuehrung.short: FORCE | build
	lualatex --output-directory=build --interaction=nonstopmode --halt-on-error 2_durchfuehrung.tex
durchfuehrung.long: build/durchfuehrung.long
build/durchfuehrung.long: FORCE |build
	TEXINPUTS=../..: \
	BIBINPUTS=../..: \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	2_durchfuehrung.tex

1_theorie.short: theorie.short
1_theorie.long: theorie.long
theorie.short: FORCE | build
	lualatex --output-directory=build --interaction=nonstopmode --halt-on-error 1_theorie.tex
theorie.long: build/theorie.long
build/theorie.long: FORCE |build
	TEXINPUTS=../..: \
	BIBINPUTS=../..: \
	max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	1_theorie.tex


build:
	mkdir -p build
clean:
	rm -rf build

FORCE:

.PHONY: all clean