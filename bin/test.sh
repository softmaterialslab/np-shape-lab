#!/bin/bash

python3 -m venv ENV

source ENV/bin/activate

pip install hublib
pip install jupyterlab
pip install testbook
pip install widgetsnbextension
pip install svgwrite
pip install pypng
pip install pytest
pip install matplotlib
pip install ipywidgets
jupyter nbextension enable --py widgetsnbextension

export ENVIRON_CONFIG_DIRS=","
pytest test || exit 1

