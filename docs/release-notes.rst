Release notes
=============

3.0.7
------

    - File extension will NOT be added automatically for output filenames. User specified filename is uses as it is.
    - Documentation improved. Parameter type hint added to every parameter
    - Chemical shift values in API dump are now converted into float before processing.

3.0.6
------

    - Bug fix: Filtering chemical shifts based on standard deviation was not working properly when you fetch chemical shifts of more than one atom type. It is fixed now.

3.0.5
------

    - New: generic 2D spectrum can optionally include chemical shifts from preceding and next residues
    - New: sequential connectivity can be shown as trace in the generic 2D spectrum
    - Documentation improved

3.0.4
------

    - function added to export peak list in csv or sparky format
    - PEP8 standard implemented
    - pytest improved
    - documentation updated

3.0.3
------

    - Class methods are converted in to Modules
    - plot title added
    - PyNMRSTAR entry object can be given as input for spectral simulation
    - parameter type definition added

3.0.2
------

    - Missing y axis label on histogram fixed
    - Output file name extensions are added after checking

3.0.1
------

    - Documentation improved
    - added few more tests
    - pytest timeout issue fixed


3.0.0
-----
First release of version 3

    - verision 1.x and 2.x will not be updated
    - pybmrb ends support for Python 2.x
    - Use  PyBMRB v3.x which runs on Python 3.6,3.7,3.8 and 3.9



