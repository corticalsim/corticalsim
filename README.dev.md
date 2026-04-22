# corticalsim3D development documentation


## Note

Testing and contributing is very welcome, especially if you can contribute with new algorithms and features.


# Development

## Building the documentation locally

Documentation is generated with [Sphinx](https://www.sphinx-doc.org/) and [Doxygen](https://www.doxygen.nl/) (for automatically generating the technical documentation) with [Breathe](https://breathe.readthedocs.io/) providing a bridge between the Sphinx and Doxygen documentation systems. To build the documentation locally follow these steps:

### Linux
Install doxygen:

```
sudo apt install doxygen
```

Then, from the root of the repository install Sphinx and the other Python dependencies:

```
python -m venv .venv
source .venv/bin/activate
pip install -r doc/requirements.txt
```

Then we generate the doxygen technical documentation xml files and use Sphinx to generate the html files with the full documentation:

```
doxygen doc/Doxyfile
sphinx-build doc doc/build
```

And we can view the generated documentation by opening `doc/build/index.md` in a browser.


# TODO

- [ ] Reproduce the installation instructions for Windows and update them as necessary.
- [ ] Add specific installation instructions for dependencies for different Linux distributions (Ubuntu, Debian, Arch...).
- [ ] Add installation instructions for MacOSX.
- [ ] Meson: Option to install the executable into a custom directory
