# Development of FindBounce package

If you would like to modify _FindBounce_ package yourself or contribute back to the original repository,
then the following instructions can help you.
First you need to install [Git](https://git-scm.com/) and
[clone](https://help.github.com/articles/cloning-a-repository/) the project
from its GitHub homepage to your local computer.

You will also need [Mathematica](https://www.wolfram.com/mathematica/) version 11.1 or later.
[WolframScript](https://www.wolfram.com/wolframscript/) is recommended for easier running of testing suite from command line.
On most systems it already comes bundled with Mathematica installation.

## Testing code

It is considered good practice that every (public) function in this package includes its own set of unit tests.
A bunch of them is collected in `Tests/Tests.wl` file, using the Mathematica testing
[framework](https://reference.wolfram.com/language/guide/SystematicTestingAndVerification.html).
It is recommended that you run them periodically during development and especially before every commit.
This can be done by calling script file `Tests/RunTests.wls` in command line
(from repository root directory) or by evaluating whole notebook `Tests/RunTests.nb`.

### Integration of tests in Git hook

Unit test can be run automatically before every commit via Git client-side
[hooks](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks).
File `pre-commit` should contain call to `Tests/RunTests.wls` script,
which exits with value 0 if all tests pass and aborts the commit otherwise.
Minimal example of `pre-commit` file content is:

    #!/bin/sh
    ./Tests/RunTests.wls

## Building the package

Package building covers procedures to convert contents of repository to `.paclet` file which can be installed by the users.
This is done by the build script included in the repository, which can be run from command line (terminal) with `wolfram -script Build.wls` .
Use option `--help` to see all available script options.
See also [console interface](https://reference.wolfram.com/language/ref/program/wolfram.html) documentation page for more information.

### Documentation

The package documentation comes in a form of a notebook with textual  explanation and executable examples.
This approach with a single notebook ensures maximal backward compatibility between versions on Mathematica.
Development version of documentation notebook contains only `"Text"` and `"Input"` cells.
`"Output"` cells are produced when the whole notebook is evaluated during documentation build procedure.
This is started with build script option `--docs` .

### Creating a .paclet file

Running the build script in command line from repository root directory
will automatically create a `FindBounce-X.Y.Z.paclet` file in the "build"
folder and also install it to `$UserBasePacletsDirectory`.
Script option `--release` additionally processes documentation notebooks and should be used for proper public releases.

### Release checklist

The following tasks should be performed for successful release of a new public version of the package.

- Check if all known bugs have been fixed.
- Check if all changes in functionality since the previous release have been documented (use `git log --oneline` to refresh memory).
- Run a full suite of automated tests.
- Build the package along with documentation and test if installed version works ok (e.g. documentation examples).
You can send this version to other test users for additional feedback.
- Update [PacletInfo.m]( PacletInfo.m ) with a new version number.
It should be determined according to [semantic versioning](https://semver.org/) principles.
- Update [README.md]( README.md ) and [CHANGELOG.md]( CHANGELOG.md ) with relevant changes.
- Make a new [`git tag`](https://git-scm.com/book/en/v2/Git-Basics-Tagging)
with corresponding version number and push it to remote repository.
- Build a public version with `--release` option of build script.
- Create a new release on GitHub [releases](https://github.com/vguada/FindBounce/releases)
page and attach there resulting `.paclet` file.
- Inform users about a new version (email, website, etc).
- Celebrate.
