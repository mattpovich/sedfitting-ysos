# sedfitting-ysos
Software package for identifying and modeling young stellar objects with mid-infrared excess emission by fitting infrared SED with the models of [Robitaille et al. (2017)](https://zenodo.org/record/166732#.XNsojy_MwWo).

This package currently consists of a custom IDL library (`pro` directory),
a couple of custom `python` wrapper modules for Tom Robitaille's [sedfitter](https://github.com/astrofrog/sedfitter), and  a set of data analysis `recipes`. 

The IDL library has been verified using IDL 8.7.1 and requires the [IDL Astronomy User's Library](https://github.com/wlandsman/IDLAstro) and [idl-coyote library](https://github.com/idl-coyote/coyote).

The Python codes have been tested on versions 3.6 and 3.7 with no adverse effects.

To get started, see the step-by-step instructions and explanations included with this package in `recipes/sedfitting_procedure_ysos.md`. 
