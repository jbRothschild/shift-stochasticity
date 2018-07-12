CODE_DIR = mte figures
current_dir = $(notdir $(shell pwd))
parent_dir = $(notdir ${current_dir}/..)
.PHONY : clean zip
# setting up suffix rules

##<<Different commands for makefile>>

## all : runs main code and creates figure
all :
	$(MAKE) -C $(CODE_DIR)

## shiftStoch : runs the numerical calculations only for shift in mte
shift :
	$(MAKE) -C mte $@

## gillespie : runs the numerical calculations only for gillespie algo in mte
gillespie :
	$(MAKE) -C mte $@

## Figs : runs all figures in figures
Figs :
	$(MAKE) -C figures

## fig% : runs figure%.py in figures
fig% :
	$(MAKE) -C figures $@

## zip : zips all files in the folder
zip :
	zip -r ${parent_dir}/${current_dir}_rothschild  ${parent_dir}/${current_dir}

## clean : remove auto-generated files
clean :
	$(MAKE) -C $(CODE_DIR) clean

help : makefile
	@sed -n 's/^##//p' $<
