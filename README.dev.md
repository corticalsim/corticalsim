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

Then Sphinx and the other Python dependencies.
From the root of the repository:

```
python -m venv ../venv-cs3d
source ../venv-cs3d/bin/activate
pip install -r doc/requirements.txt
```

Then we first generate the doxygen technical documentation xml files:

```
mkdir -p build/docs/doxygen
cd doc
doxygen
cd ..
```

Finally we use Sphinx to generate the html files with the integrated documentation

```
sphinx-build doc build/docs
```

And we can view the generated documentation by opening the file in `build/docs/


# TODO

- [ ] Reproduce the installation instructions for Windows and update them as necessary.
- [ ] Add specific installation instructions for dependencies for different Linux distributions (Ubuntu, Debian, Arch...).
- [ ] Add installation instructions for MacOSX.
- [ ] Meson: Option to install the executable into a custom directory
