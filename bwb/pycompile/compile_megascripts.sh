#!/bin/bash

# Build megaplot.spec file
pyinstaller megaplot.py --strip --onefile

# EDIT megaplot.spec (if necessary)

# Now compile using the spec file
pyinstaller megaplot.spec

# Build megasky executable:
pyinstaller megasky.py --additional-hooks-dir ./hooks --strip --onefile
