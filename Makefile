
#--------------------------------
# Variables
#--------------------------------

SIZES := 1 2 3 4 5 6 7 8

#--------------------------------
# Constants
#--------------------------------

PYFILES :=
PYFILES += csf.py
PYFILES += fragment.py
PYFILES += hess.py
PYFILES += makedeps.py
PYFILES += path.py
PYFILES += perm.py
PYFILES += util.py

OUTFILES :=

#--------------------------------
# Top-level targets
#--------------------------------

all: output.py

clean:
	git clean -dfx

test:
	python -m doctest $(PYFILES)

#--------------------------------
# Included makefiles
#--------------------------------

include $(SIZES:%=var/size-%.d)

#--------------------------------
# Internal targets
#--------------------------------

output:
	mkdir output

var:
	mkdir var

output.py: output-preamble.py $(OUTFILES)
	cat output-preamble.py $(OUTFILES) >output.py

var/size-%.d: | var
	python makedeps.py $* >$@

var/csf-size-%: $(PYFILES) | output var
	python csf.py $*
	touch $@

output/hess-%.py: $(PYFILES) | output
	python hess.py $*

.PHONY: all clean test

