auswertung.tex: build/auswertung.pdf

build:
	mkdir -p build

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

clean:
	rm -rf build

FORCE:

.PHONY: all clean