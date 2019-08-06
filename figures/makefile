ALL_FIGS = $(sort $(patsubst figure%.py, fig%, $(wildcard figure*.py)))
current_dir = $(notdir $(shell pwd))
parent_dir = $(notdir ${current_dir}/..)
.PHONY : clean
# setting up suffix rules

##<<Different commands for makefile>>
## all : runs all figures
all : $(ALL_FIGS)

## fig% : runs figure%.py
fig% : figure%.py
	python2 $<


## clean : remove auto-generated files
clean :
	rm -f *.pgf *.pdf *.pyc

help : makefile
	@sed -n 's/^##//p' $<
