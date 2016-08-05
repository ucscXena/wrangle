## xena query file
wget https://raw.githubusercontent.com/ucscXena/ucsc-xena-server/master/python/xena_query.py -O xena_query.py

## install notebook
install notebook by Anaconda
https://www.continuum.io/downloads

install new kernal instruction
http://ipython.readthedocs.io/en/stable/install/kernel_install.html

somehow you need to change the kernal  on the UI

## run notebook on laptop

* jupyter notebook

find the notebook file, double click to open the nootbook

## using notebook

good tutorial:
https://geosci.uchicago.edu/~rtp1/PrinciplesPlanetaryClimate/Python/NotebookQuickstart/InstantNotebooks.html


## Running a notebook on a public server configuration

see http://stackoverflow.com/questions/31962862/ipython-ipython-notebook-config-py-missing

* jupyter notebook --generate-config

this should create a ~/.jupyter/jupyter_notebook_config.py file with relevant notebook settings ready

* c = get_config()
* c.NotebookApp.ip = '*'
* c.NotebookApp.open_browser = False
* c.NotebookApp.port = 9999

see https://ipython.org/ipython-doc/3/notebook/public_server.html

* AWS port configuration

## integrate with web page

https://github.com/oreillymedia/thebe
