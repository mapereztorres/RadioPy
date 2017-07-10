# RadioPy
Python package for synchrotron radio emission in Radio Astrophysics


Welcome to RadioPy : the python package for radio astronomy!

This package is an open source software intended to help researchers
with the data-processing for synchrotron emission at radio frequency,
created by ultrarelativistic electrons. Astrophysical events regarding
this particular form of electromagnetic emission could be
radiogalaxies, quasares, or super novae, among others.

License
-------
RadioPy is open source. See file LICENSE.txt

Installation
------------
RadioPy needs the installation of other packages for full and complete
functioning. More information about those packages and how to install
them can be found in this file:
    /radiopy/docs/installation/radiopy_install_guide.txt



Documentation
-------------
VISUALIZATING RADIOPY DOCUMENTATION

    /radiopy/docs/_build/html/index.html

CREATING YOUR OWN HTML DOCFILES
A guide for RadioPy documentation:

    /radiopy/docs/radiopy_documentation_guide.txt



Features
--------

RadioPy uses the Synchrotron Self-Absorption (SSA) model to calculate
three important magnitudes about astrophysical synchrotron radio
emission data: the radious, the intensity of the magnetic field that
drives the electrons to ultrarelativistic velocities and the density
of the super novaes (SN) circumstellar wind.  RadioPy's version 1 is
specialized in the radio emission in SN, though the SSA model can be
used to describe other phenomena such as far radiogalaxies or quasars,
among others.

The package is able to automathiclly read flux vs frequency data from
a .csv file*, process the information and create a report with the
results of the analysis, both in html and pdf formats. But it is not
only a straight-forward tool since its Object Oriented Programming
(OOP) architecture grants its methods a great flexibility and a high
degree of customization, allowing researches to easily create their
own models and simulations. Thus, the architechture offers modularity,
a desired characteristic since this package is launched to evolve and
to integrate more and better models.

*Such file must have a form specified in the file:
    /radiopy/docs/radiopy_csv_format_guide.txt
