
OUTPUTDIR := output

PYFILES :=
PYFILES += csf.py
PYFILES += fragment.py
PYFILES += hess.py
PYFILES += path.py
PYFILES += perm.py
PYFILES += util.py

all: | $(OUTPUTDIR)

clean:
	git clean -dfx

test:
	python -m doctest $(PYFILES)

$(OUTPUTDIR):
	-mkdir $(OUTPUTDIR)

.PHONY: all clean

