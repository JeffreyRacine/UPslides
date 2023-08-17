talk:

	quarto render index.qmd
	Rscript -e 'knitr::purl("index.qmd")'
	Rscript -e 'knitr::purl("density_probability.qmd")'
	Rscript -e 'knitr::purl("regression.qmd")'
	Rscript -e 'knitr::purl("ar_volatility_panel.qmd")'
	Rscript -e 'knitr::purl("semipar.qmd")'
	Rscript -e 'knitr::purl("model_performance.qmd")'
	Rscript -e 'knitr::purl("quarto.qmd")'
	cat index.R density_probability.R regression.R ar_volatility_panel.R semipar.R model_performance.R quarto.R > tmp_R;rm *.R;mv tmp_R index.R
	git add index_files

chrome:

	@open -a Google\ Chrome.app index-speaker.html

present:

	@open -a RStudio.app
	@open -a Google\ Chrome.app index-speaker.html
	@clear
	@echo "Presentation Tips"
	@echo " "
	@echo "- Cmd-F1 to mirror/unmirror external display (first plug in external display)"
	@echo "- Ctrl-Cmd-F for fullscreen of apps (“Fullscreen”)"
	@echo "- Ctrl-Cmd-G for moving an app to external display (“Go”)"
	@echo " "
	@echo "In revealjs presentation"
	@echo " "
	@echo "Z - Zoom"
	@echo "E - Export to PDF"
	@echo "F - Fullscreen"
	@echo "B - Blackboard"
	@echo "C - Crayon"
	@echo " "
	@echo "Start speaker view, push slides to external display (safest venue option) "

clean: 

	rm -rf index_files index_cache *~ *.bak index.html index-speaker.html
