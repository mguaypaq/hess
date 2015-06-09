
OUTPUTDIR := output
PYFILES := \
    csf.py \
    fragment.py \
    hess.py \
    path.py \
    perm.py \
    util.py \

all: | $(OUTPUTDIR)

clean:
	git clean -dfx

test:
	python -m doctest $(PYFILES)

$(OUTPUTDIR):
	-mkdir $(OUTPUTDIR)

.PHONY: all clean
.DELETE_ON_ERROR:
