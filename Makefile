# Make recipies for NCI NDR project

PYTHON = python
ENV = env
BASHLIKE_SHELL_COMMAND := cmd.exe /C
ENV_SCRIPTS = $(ENV)\Scripts
ENV_ACTIVATE = $(ENV_SCRIPTS)\activate
MAKE := make

PIP = $(PYTHON) -m pip

env:
	$(PYTHON) -m virtualenv --system-site-packages $(ENV)
	$(BASHLIKE_SHELL_COMMAND) "$(ENV_ACTIVATE) && $(PIP) install -r requirements.txt"
	$(BASHLIKE_SHELL_COMMAND) "$(ENV_ACTIVATE) && $(MAKE) install"