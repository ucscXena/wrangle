## xena query file
wget https://raw.githubusercontent.com/ucscXena/ucsc-xena-server/master/python/xena_query.py -O xena_query.py

## install notebook
install notebook by Anaconda
https://www.continuum.io/downloads

optional install a different kernal instruction.
http://ipython.readthedocs.io/en/stable/install/kernel_install.html

## run notebook on laptop

* jupyter notebook

or

* screen
* jupyter notebook

## using notebook

good tutorial:
https://geosci.uchicago.edu/~rtp1/PrinciplesPlanetaryClimate/Python/NotebookQuickstart/InstantNotebooks.html


## Running a notebook on a public server configuration


http://jupyter-notebook.readthedocs.io/en/latest/public_server.html
-- see http://stackoverflow.com/questions/31962862/ipython-ipython-notebook-config-py-missing

for set up password
see https://ipython.org/ipython-doc/3/notebook/public_server.html

* jupyter notebook --generate-config

this should create a ~/.jupyter/jupyter_notebook_config.py file with relevant notebook settings ready

c.NotebookApp.ip = '*'
c.NotebookApp.open_browser = False
c.NotebookApp.port = 9999
c.NotebookApp.notebook_dir = u'cgDataNew/Nof1Analysis/'
c.NotebookApp.password = u'sha1:bcd...<your hashed password here>'

#c.NotebookApp.certfile = u'/absolute/path/to/your/certificate/mycert.pem'
#c.NotebookApp.keyfile = u'/absolute/path/to/your/certificate/mykey.key'



* AWS ELB port configuration

Not sure, sockets use TCP/SSL protocol not HTTP so you will need to change it to TCP/SSL in order to make web sockets work.
https://www.bountysource.com/issues/5646899-aws-ec2-behind-elb-always-prints-error-unexpected-response-code-400

## integrate with web page

https://github.com/oreillymedia/thebe
