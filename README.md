# MDAnalysis CLI interface

A simple command line interface (CLI) to the analysis classes of [MDAnalysis](https://www.mdanalysis.org) 
using argparse and [numpydoc](https://numpydoc.readthedocs.io/en/latest/). 
This project is in a really early stage and 
work in progress. Currently it is not connected to MDAnalysis.

To test the functionality clone this repository and install the forked version of MDAnalysis

```bash
git clone git@github.com:PicoCentauri/mdanalysis.git
cd mdanalysis/package
git checkout cli
python setup.py install
```

Afterwards you can run

```bash
python cli_main.py -h
```
