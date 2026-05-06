# corticalsim3D development documentation


## Note

Testing and contributing is very welcome, especially if you can contribute new algorithms and features.


# Development

## Building the documentation locally

Documentation is generated with [Sphinx](https://www.sphinx-doc.org/) and [Doxygen](https://www.doxygen.nl/) (for automatically generating the technical documentation), with [Breathe](https://breathe.readthedocs.io/) providing a bridge between the Sphinx and Doxygen documentation systems. To build the documentation locally, follow these steps:

### Linux
Install doxygen:

```
sudo apt install doxygen
```

Then, install Sphinx and the other Python dependencies from the root of the repository:

```
python -m venv .venv
source .venv/bin/activate
pip install -r doc/requirements.txt
```

Then we generate the Doxygen technical documentation XML files and use Sphinx to generate the HTML files with the full documentation:

```
doxygen doc/Doxyfile
sphinx-build doc doc/build
```

And we can view the generated documentation by opening `doc/build/index.html` in a browser.