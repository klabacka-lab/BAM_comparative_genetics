VENV   := .venv
PYTHON := $(VENV)/bin/python
PIP    := $(VENV)/bin/pip

INPUT     ?= ../HyPhy/HyPhy_FEL_output.txt
SELECTION ?=
ALPHA_MIN ?=
ALPHA_MAX ?=
BETA_MIN  ?=
BETA_MAX  ?=
LRT_MIN   ?=
LRT_MAX   ?=
PARTITION ?=

FLAGS :=
ifdef PVALUE
FLAGS += --pvalue $(PVALUE)
endif
ifdef OMEGA
FLAGS += --omega $(OMEGA)
endif
ifdef SELECTION
FLAGS += --selection $(SELECTION)
endif
ifdef ALPHA_MIN
FLAGS += --alpha-min $(ALPHA_MIN)
endif
ifdef ALPHA_MAX
FLAGS += --alpha-max $(ALPHA_MAX)
endif
ifdef BETA_MIN
FLAGS += --beta-min $(BETA_MIN)
endif
ifdef BETA_MAX
FLAGS += --beta-max $(BETA_MAX)
endif
ifdef LRT_MIN
FLAGS += --lrt-min $(LRT_MIN)
endif
ifdef LRT_MAX
FLAGS += --lrt-max $(LRT_MAX)
endif
ifdef PARTITION
FLAGS += --partition $(PARTITION)
endif

CODON_FILE ?= output/codon_positions.txt
CIF_FILE   ?= output/cif_residues.txt

.PHONY: install wrangle map vis clean

$(VENV)/.installed: requirements.txt
	python3 -m venv $(VENV)
	$(PIP) install -q -r requirements.txt
	touch $(VENV)/.installed

install: $(VENV)/.installed

wrangle: install
ifeq ($(PVALUE)$(OMEGA),)
	$(error At least one filter is required. Usage: make wrangle PVALUE=0.05  or  make wrangle OMEGA=0.5)
endif
	mkdir -p output
	$(PYTHON) parse_fel.py $(INPUT) -o $(CODON_FILE) $(FLAGS)

map: install
	$(PYTHON) map_codons.py $(CODON_FILE) -o $(CIF_FILE)

vis: install
	$(VENV)/bin/streamlit run fel_vis_app.py -- $(INPUT)

clean:
	rm -rf $(VENV) __pycache__ output FEL_sites.png
