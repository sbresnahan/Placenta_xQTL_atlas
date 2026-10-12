.PHONY: build prepare clean publish

build:
	Rscript render_book.R

prepare:
	Rscript prepare_book.R

clean:
	rm -rf .staging _book _bookdown_files _main.Rmd *.utf8.md *.knit.md

publish:
	./publish_gh_pages.sh
