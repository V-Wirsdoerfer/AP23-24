auswertung.tex: build/auswertung.pdf

build:
	mkdir -p build

auswertung.short: FORCE | build
	lualatex --output-directory=build --interaction=nonstopmode --halt-on-error auswertung.tex
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
	auswertung.tex


diskussion.short: FORCE | build
	lualatex --output-directory=build --interaction=nonstopmode --halt-on-error diskussion.tex
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
	diskussion.tex


durchfuehrung.short: FORCE | build
	lualatex --output-directory=build --interaction=nonstopmode --halt-on-error durchfuehrung.tex
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
	durchfuehrung.tex


theorie.short: FORCE | build
	lualatex --output-directory=build --interaction=nonstopmode --halt-on-error theorie.tex
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
	theorie.tex

clean:
	rm -rf build

FORCE:

.PHONY: all clean