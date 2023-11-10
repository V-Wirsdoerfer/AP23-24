auswertung.tex: build/auswertung.pdf

build/auswertung.pdf:
	lualatex --output-directory=build --interaction=nonstopmode --halt-on-error auswertung.tex


diskussion.tex: build/diskussion.pdf

build/diskussion.pdf:
	lualatex --output-directory=build --interaction=nonstopmode --halt-on-error diskussion.tex



durchfuehrung.tex: build/durchfuehrung.pdf

build/durchfuehrung.pdf:
	lualatex --output-directory=build --interaction=nonstopmode --halt-on-error durchfuehrung.tex



theorie.tex: build/theorie.pdf

build/theorie.pdf:
	lualatex --output-directory=build --interaction=nonstopmode --halt-on-error theorie.tex


build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
