Project Vision
==============

Problem definition
------------------

Powerful libraries for the analysis of Molecular Dynamics’ trajectories
and derived data are being actively developed; a `Python`_ interface
serves the most known libraries: `MDTraj`_, `MDAnalysis`_, `OpenMM`_, `ProDy`_, `nglview`_, and others.

Many researchers with strong Biochemical and Biophysics background use
Molecular Dynamics (MD) to generate *ab initio* models or to module
their *in vitro* or *in vivo* acquired data. This collective of
researchers use MD as a tool and not necessary have (or need to have)
programming/scripting skills.

*Problem 1:* There is a growing gap in accessibility between the
libraries being developed and the non-developer users that use them. 

*Problem 2:* Data analysis rapidly reach routine. Providing pre-defined
routines for data mining, parameter calculation and plotting,
considerably speeds up analysis of multiple datasets.

Proposed Solutions
------------------

The gap described in *problem 1* can be fulfilled either by teaching
programming skills to this research collective, asking them to be
self-taught, or by providing an interface for them to use the MD
analysis libraries without the need to learn how to program.

*Solution 1:* Tauren-MD aims to provide such interface for non-developer users, i.e., users without
programming skills.

Also for users with programming skills, repetitive data analysis is time consuming, and even
challenging is one considers the time required to assemble and organize routines that automate analysis steps.

*Solution 2:* Tauren-MD aims to provide predefined set routines to
stream-line analysis procedures and, therefore, boost deliverables’
rate.

Implementation
==============

Versioning
----------

Tauren-MD follows `semantic version rules`_.

User Requirements
-----------------

-  Users should be able to access all Tauren-MD routines, hereafter
   named ``actions``, through a configurable text file, hereafter named
   ``config``.
-  to configure the ``config`` file NO programming skills should be
   required. The minimum requirement being the syntax of the ``config``
   file itself.
-  the ``config`` file should be a list-like menu of the ``actions``
   that the user wants to perform.
-  ``actions`` entries should contain two fields: 1) its name and 2) its
   arguments (options).
-  ``actions`` can be flagged (with ``#`` character) to (de)activate its
   execution.
-  ``actions`` can be deactivated by removing them from the list.
-  ``actions`` can be repeated.
-  ``actions`` can be reordered.
-  ``config`` files can optionally provide a ``path`` to the trajectory
   file.
-  ``config`` files can optionally provide a ``path`` to the topology
   file.
-  Tauren-MD execution is provided by a ``bin`` executable file.
-  Tauren-MD executable takes four optional arguments:

   -  a ``path`` to a config file
   -  a ``path`` to a trajectory file
   -  a ``path`` to a topology file
   -  the MD library the user wishes to operate with

-  for optional arguments not provided, a default value is used (see
   bellow).

Architecture
------------

It should be possible to use Tauren-MD as an interface for non-developer users as well as a well documented Python library. With that in mind:

Installation and dependency management
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tauren-MD installation should be possible via the `Tree-of-Life project`_, in order to facilitate this step to  non-developers.

Additional methods cen eventually be implemented, for example, PiPy or Conda distribution. Nonetheless, installation via Tree-of-Life should always be possible.

Tauren-MD executable file
~~~~~~~~~~~~~~~~~~~~~~~~~

Tauren-MD executable file should:

-  Use ``paths`` from ``--trajectory`` and ``--topology`` arguments when
   provided;

   -  otherwise use ``paths`` provided in the ``config`` file,
   -  if none provided, ``raise error``.

-  use ``config`` file provided in arguments;

   -  if none provided use ``config`` distributed by default with
      package,
   -  default ``config`` should just load a ``trajectory``.

-  read over the ``config``\ ’s list of ``actions`` and execute them in
   order with the respective ``action``\ ’s arguments.

   -  communication between ``config`` file ``action``\ ’s name and the
      ``action``\ ’s function should be provided by a dictionary in
      ``tauren._interface`` module.

Tauren-MD actions
~~~~~~~~~~~~~~~~~

Tauren-MD ``actions``:

-  are functions defined in Tauren-MD libs.
-  can be configured with ``kwargs``;

   -  these ``kwargs`` can be used in the ``config`` file ``action``
      entries.

Lib organization
~~~~~~~~~~~~~~~~

General project organization:

-  ``TaurenMD/``: *main* library folder

   -  ``tauren.py``: contains Tauren Trajectory classes.
    
      - ``TaurenTraj`` is an abstract class that provides an interface between the user and the MD analysis libs. **Polymorphism** *reigns* here.
      - Classes implementing methods from the specific MD analysis libs should inherit from ``TaurenTraj`` and implement its *abstractmethods* accordingly. Additional dunders can and should be implemented if needed.
      - There should be a dunder *abstractmethod* in ``TaurenTraj`` for every method there implemented that is an interface to a specific MD lib operation.
      - Documented interface should be provided **only** by ``TaurenTraj``, from here, the program workflow is directed to the corresponding hidden methods of the specific subclasses thta configure the algorithms for the different MD analysis libs implemented. In this way, subclasses inheriting from ``TaurenTraj`` should provide **only** dunder methods. Exceptions are allowed for **@properties**.
      - In this way, documentation provided in ``TaurenTraj`` methods should cover all subclass implementations.
      - If a given subclass does not provides a routine for a given method in ``TaurenTraj``, a ``Not implemented`` message should be logged to warn the user and that action ignored.
      - **@staticmethods** and **@classmethods** should be testable by their *return* value.

   -  ``core.py``: variables and functions used system wide, commons, decorators,
      helpers…
   -  ``_interface.py``: interface between the Tauren-MD configuration file and
      Tauren-MD library.
   -  ``logger.py``: the Tauren-MD logging interface.
   -  ``load.py``: module with functions used to load data from outside
      Tauren-MD.
   -  ``produce.py``: this module provides functions that concatenate related operations in logical series, for example, calculating, exporting and plotting a given type of data.
   -  ``plot.py``: plotting templates.

      -  plotting routines should be functions in modules and NOT
         methods in classes.
      -  parameters to plots should be provided by kwargs.
         
.. _Python: https://www.python.org/
.. _MDTraj: https://github.com/mdtraj/mdtraj
.. _MDAnalysis: https://www.mdanalysis.org/
.. _Prody: http://prody.csb.pitt.edu/index.html
.. _nglview: https://github.com/arose/nglview
.. _OpenMM: https://github.com/pandegroup/openmm
.. _Tree-of-Life project: https://github.com/joaomcteixeira/Tree-of-Life
.. _semantic version rules: https://semver.org/
