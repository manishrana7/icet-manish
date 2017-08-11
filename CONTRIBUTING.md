Contribution guidelines
=======================

Please use spaces
-----------------

While you are entitled to [your own
opinion](http://lea.verou.me/2012/01/why-tabs-are-clearly-superior/) this
project uses spaces instead of tabs. Even if you are geeky enough to care and
like [Silicon Valley](https://www.youtube.com/watch?v=SsoOG6ZeyUI) you should
know that [developers who use spaces make more
money](https://stackoverflow.blog/2017/06/15/developers-use-spaces-make-money-use-tabs/).
Also the use of spaces is strongly recommended by our beloved
[pep8](https://www.python.org/dev/peps/pep-0008/) standard.


C++
---

This project adopts a more concise style when writing C++ code that is in the
spirit of [K&R](https://en.wikipedia.org/wiki/Indent_style) and
[pep7](https://www.python.org/dev/peps/pep-0007/).

Any functions/functionality *must* be properly documented. The API documentation
is generated using [doxygen](http://www.stack.nl/~dimitri/doxygen/). You should
therefore include comment blocks that document your code and are formatted to
comply with the doxygen markup style.

For most functions, class members, etc. that can be comprehensively described
using a single line one can use the triple-slash form, e.g.,:

   /// Space group number according to ITCA.
   int spacegroup;
   
   /// Returns the space group number.
   int getSpacegroup() { return spacegroup; }

For more substantial functions, classes, or other elements (such as ``enum``-s)
adopt the extended documentation form:

   /**
   @brief Write a structure to file.
   @details This function writes an atomic structure to file using different formats.
   @param struct   The atomic configuration
   @param filename The name of the output file.
   @param format   The output file format; possible values: 'vasp', 'xyz'
   */
   void writeStructureToFile(AtomicStructure *struct, string::string filename, string::string format) {
     ...
}

More examples can of course be found in the code.

Please use [CamelCase](https://en.wikipedia.org/wiki/Camel_case) and expressive
names for functions, classes, and members (avoiding unnecessary and non-standard
abbreviations), e.g. ``writeStructureToFile``, ``AtomicStructure``. Private and
protected class members should be preceded by an underscore as in
``_myPrivateVariable``.

Python
------

Code should be [pep8](https://www.python.org/dev/peps/pep-0008/)
compliant and pass
[pyflakes](https://pypi.python.org/pypi/pyflakes). (Eventually,
pep8/pyflakes will be added to the CI, at which point code *must* be
compliant.)

Any functions/functionality *must* be properly documented. This
includes [docstrings](https://en.wikipedia.org/wiki/Docstring) for
functions, classes, and modules that clearly describe the task
performed, the interface (where necessary), the output, and if
possible an example. atomicrex uses [NumPy Style Python
Docstrings](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html).

When in doubt ask the main developers. Also [the coding conventions from
ASE](https://wiki.fysik.dtu.dk/ase/development/python_codingstandard.html)
provide useful guidelines.


Commits/issues/merge requests
-----------------------------

When writing commit messages, generating issues, or submitting merge
requests, please use the following codes to identify the category of
your commit/issue/merge request.

* BLD: change related to building  
* BUG: bug fix
* DATA: general data
* DOC: documentation  
* ENH: enhancement  
* MAINT: maintenance commit (refactoring, typos, etc.)  
* STY: style fix (whitespace, PEP8)  
* TST: addition or modification of tests  
* REL: related to releases  

Less common:
* API: an (incompatible) API change  
* DEP: deprecate something, or remove a deprecated object  
* DEV: development tool or utility  
* FIG: images and figures
* REV: revert an earlier commit  

The first line should not exceed 78 characters. If you require more
space, insert an empty line after the "title" and add a longer message
below. In this message, you should again limit yourself to 78
characters *per* line.

Hint: If you are using emacs you can use ``Meta``+``q``
[shortcut](https://shortcutworld.com/en/Emacs/23.2.1/linux/all) to
"reflow the text". In sublime you can achieve a similar effect by using
``Alt``+``q`` (on Linux)/ ``Alt``+``Cmd``+``q`` (on MacOSX).
