# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXPROJ    = planetplanet
SOURCEDIR     = .
BUILDDIR      = .build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
# Copy title image
	@cp title.png .build/html/_images/title.png
	@cp PPOs.pdf .build/html/PPOs.pdf
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
ifdef PLANETPLANET_DOCS
# Export to server
	#@mkdir -p .build/html/_images/
	#@rsync -aP --times --delete --exclude='*.pdf' --exclude='*.hires.png' .build/html/ $(PLANETPLANET_DOCS)
	#@rsync --times .build/html/PPOs.pdf $(PLANETPLANET_DOCS)
	@echo "Export finished."
endif