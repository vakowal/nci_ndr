# Make recipies for NCI NDR project

PYTHON = python36
ENV = env
BASHLIKE_SHELL_COMMAND := cmd.exe /C
ENV_SCRIPTS = $(ENV)\Scripts
ENV_ACTIVATE = $(ENV_SCRIPTS)\activate
MAKE := make
RM := rm -r
RMDIR := $(RM)

PIP = $(PYTHON) -m pip

env:
	$(PYTHON) -m venv $(ENV)
	$(BASHLIKE_SHELL_COMMAND) "$(ENV_ACTIVATE) && $(PIP) install "Shapely-1.6.4.post2-cp36-cp36m-win_amd64.whl""
	$(BASHLIKE_SHELL_COMMAND) "$(ENV_ACTIVATE) && $(PIP) install "GDAL-2.4.1-cp36-cp36m-win_amd64.whl""
	$(BASHLIKE_SHELL_COMMAND) "$(ENV_ACTIVATE) && $(PIP) install "Rtree-0.8.3-cp36-cp36m-win_amd64.whl""
	$(BASHLIKE_SHELL_COMMAND) "$(ENV_ACTIVATE) && $(PIP) install numpy"
	$(BASHLIKE_SHELL_COMMAND) "$(ENV_ACTIVATE) && $(PIP) install pandas"    
	$(BASHLIKE_SHELL_COMMAND) "$(ENV_ACTIVATE) && $(PIP) install pygeoprocessing"

clean:
	-$(RMDIR) $(ENV)