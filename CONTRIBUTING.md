# Development of FindBounce package

If you would like to modify _FindBounce_ package yourself or contribute back to the original repository,
then the following instructions can help you.
First you need to install [Git](https://git-scm.com/) version control system and
[clone](https://help.github.com/articles/cloning-a-repository/) the project
from its GitHub homepage to your local computer.

## Required tools

The main required tool is [Mathematica](https://www.wolfram.com/mathematica/) version 11.0 or later,
although version 12.0 is recommended because it produces better looking documentation.

You will need command line tool [WolframScript](https://www.wolfram.com/wolframscript/)
for running the test suite and building the package.
On most systems it already comes bundled with Mathematica installation.

Also [Wolfram Workbench](https://www.wolfram.com/workbench/) is needed to
process documentation notebooks for proper integration into Mathematica
documentation framework.
Wolfram Workbench is a plugin for [Eclipse](https://www.eclipse.org/), a popular and robust IDE.
Please follow these [instructions](http://support.wolfram.com/kb/27221) for installation
and pay attention to choose "Eclipse IDE for Java Developers". The procedure was succesfully
carried out in the Eclipse version 4.13.0.

## Testing code

It is considered good practice that every (public) function in this package includes its own set of unit tests.
A bunch of them is collected in `Tests/Tests.wl` file, using the Mathematica testing
[framework](https://reference.wolfram.com/language/guide/SystematicTestingAndVerification.html).
It is recommended that you run them periodically during development and especially before every commit.
This can be done by calling script file `Tests/RunTests.wls` in command line
(from repository root directory) or by evaluating whole notebook `Tests/RunTests.nb`.

Make sure that all tests pass on the earliest supported version of Mathematica (currently version 10.0).

### Integration of tests in Git hook

Unit test can be run automatically before every commit via Git client-side
[hooks](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks).
File `pre-commit` should contain call to `Tests/RunTests.wls` script,
which exits with value 0 if all tests pass and aborts the commit otherwise.
Minimal example of `pre-commit` file content is:

    #!/bin/sh
    ./Tests/RunTests.wls

## Building the package

Package building is a __two step__ procedure.
First, documentation notebooks have to be processed with Wolfram Workbench.
Then, the processed documentation is packaged with other source files into
a `.paclet` file which can be installed by the users.
This second step is done by the build script included in the repository,
which can be run from command line (terminal) with `Build.wls`.
On Unix-like operating systems it may be necessary to set executable permissions first
with `chmod a+x Build.wls`.
For more information on how to run scripts see official `wolframscript`
[documentation](https://reference.wolfram.com/language/ref/program/wolframscript.html).

### Documentation processing

The package documentation comes in a form of notebooks with textual explanation
and executable code examples.
Source documentation notebooks are created from templates in Wolfram Workbench and
have to be processed for proper integration into Mathematica documentation framework.
During processing, for example [`Graphics`](https://reference.wolfram.com/language/ref/Graphics.html)
objects are rasterized to minimize documentation file size.

Before manipulating documentation for the first time, you have to import
the project into Workbench. This procedure needs to be done only once.

* From the menu bar, select __File > Open Project from File System__.
* In the Import wizard, select _FindBounce_ repository root directory for
the input field "Import source:".
* You can leave all checkboxes as default and click __Finish__.
* _FindBounce_ directory tree now appears in "Package Explorer" tab (on the left side).

For more information about adding and modifying documentation notebooks, please see official
Wolfram Workbench [documentation](https://reference.wolfram.com/workbench/index.jsp).

To process documentation follow the steps bellow.
This has to be done for every package build.

* From the Workbench menu bar, select __Window > Show View > Application Tools__.
(If the last choice is not present, find it under section
__Window > Show View > Others > Wolfram > Application Tools__.)
* In the "Application Tools" tab select _FindBounce_ in project drop down list.
* Click on button "Build", just bellow the drop down list.
* Procedure may take a few minutes and it copies processed notebooks to
"build" subdirectory in repository root directory.
* You can also see the result of procedure after it is completed
by clicking button "Preview" in "Application Tools" tab.

### Creating a .paclet file

Running the build script `Build.wls` in command line (terminal) from repository root directory
will automatically create a `FindBounce-X.Y.Z.paclet` file in the "build"
folder.
Keep in mind that you have to process documentation with Workbench __before__
running the build script.

Script option `--install` will also install it locally to `$UserBasePacletsDirectory`.

### Release checklist

The following tasks should be performed for successful release of a new public version of the package.

* Check if all known bugs have been fixed.
* Check if all changes in functionality since the previous release have been documented
(use `git log --oneline` to refresh memory).
* Run a full suite of automated tests.
* Build the package along with documentation and test if installed version works ok (e.g. documentation examples).
You can send this version to other test users for additional feedback.
* Update [PacletInfo.m]( PacletInfo.m ) with a new version number.
It should be determined according to [semantic versioning](https://semver.org/) principles.
* Update [README.md]( README.md ) and [CHANGELOG.md]( CHANGELOG.md ) with relevant changes.
* Make a new [`git tag`](https://git-scm.com/book/en/v2/Git-Basics-Tagging)
with corresponding version number and push it to remote repository.
* Build the package again. This is the public version.
* Create a new release on GitHub [releases](https://github.com/vguada/FindBounce/releases)
page and attach there resulting `.paclet` file.
* Inform users about a new version (email, website, etc).
* Celebrate.
