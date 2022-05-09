# Contributing to ESPResSo
[contributing-to-espresso]: #contributing-to-espresso

Thank you for your interest in contributing to ESPResSo! There
are many ways to contribute and we appreciate all of them. This
document is quite lengthy, so here is a list of links to the
major sections:

* [New Features](#new-features)
* [Bug Reports](#bug-reports)
* [Pull Requests](#pull-requests)
* [Writing Documentation](#writing-documentation)

If you have questions, please make a post on [the users
mailing list][mailing-list].

[mailing-list]: https://lists.nongnu.org/mailman/listinfo/espressomd-users

## New Features
[new-features]: #new-features

If you plan to create a new feature for ESPResSo, it is
*absolutely essential* that you first post your plan on the
[users mailing list][mailing-list].  This way we can give you
useful tips on how to go about a successful integration and can
point out possibly reusable code.  Maybe your feature even exists
already and you just didn't notice it, so we want to save you the
effort.

## Bug Reports
[bug-reports]: #bug-reports

We cannot fix bugs we don't know about, so please report liberally.
If you are unsure whether something is a bug or not, feel free to
ask on the [users mailing list][mailing-list] or on the
[discussion forum](https://github.com/espressomd/espresso/discussions).

To fix your problem quickly it is essential that we have an easy
way to reproduce your bug.  Therefore you should include in your
report the *desired behavior*, a *specific problem or error* and
*the shortest script necessary* to reproduce it.  We have handy
step-by-step instructions in the wiki page [Filing bug reports
](https://github.com/espressomd/espresso/wiki/Filing-bug-reports).

## Pull Requests
[pull-requests]: #pull-requests

Pull requests are the primary mechanism we use to change ESPResSo.
GitHub itself has some [great documentation][gh-pull-requests]
on using the Pull Request feature.  We use the "fork and pull"
model [described here][gh-development-models], where contributors
push changes to their personal fork and create pull requests to
bring those changes into the source repository.

[gh-pull-requests]: https://help.github.com/articles/about-pull-requests/
[gh-development-models]: https://help.github.com/articles/about-collaborative-development-models/

Please make pull requests against the `python` branch.

We have a continuous integration system for [automated testing](
https://github.com/espressomd/espresso/wiki/Code-review#automated-code-review),
which runs on all pull requests to see whether new additions play
nicely with the existing code.  Because it takes a long time for
all the tests to be completed, you might want to run `make check`
with the settings from `maintainer/configs/maxset.hpp` locally
first.

All pull requests are reviewed by one or more of the ESPResSo
core team members, as explained in the wiki page [Code review
](https://github.com/espressomd/espresso/wiki/Code-review#Human_code_review).

## Writing Documentation
[writing-documentation]: #writing-documentation

Documentation improvements are very welcome.  The source of it is
located in `doc/sphinx` in the tree, and standard API documentation
is generated from the source code itself.

To find documentation-related issues, filter using the
[Documentation label][issues-doc].

[issues-doc]: https://github.com/espressomd/espresso/issues?q=is%3Aissue+is%3Aopen+label%3ADocumentation
