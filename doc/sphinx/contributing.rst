Getting involved
================

Up to date information about the development of can be found at the web
page http://espressomd.org As the important information can change in
time, we will not describe its contents in detail but rather request the
reader to go directly to the URL. Among other things, one can find
information about the following topics there:

-  FAQ

-  Latest stable release of and older releases

-  Obtaining development version of

-  Archives of both developers’ and users’ mailing lists

-  Registering to mailing lists

-  Submitting a bug report

Community support and mailing lists
-----------------------------------

If you have any questions concerning which you cannot resolve by
yourself, you may post a message to the mailing list. Instructions on
how to register to the mailing lists and post messages can be found on
the homepage http://espressomd.org. Before posting a question and
waiting for someone to answer, it may be useful to search the mailing
list archives or FAQ and see if you can get the answer immediately. For
several reasons it is recommended to send all questions to the mailing
lists rather than to contact individual developers:

-  All registered users get your message and you have a higher
   probability that it is answered soon.

-  Your question and the answers are archived and the archives can be
   searched by others.

-  The answer may be useful also to other registered users.

-  There may not be a unique answer to your problem and it may be useful
   to get suggestions from different people.

Please remember that this is a community mailing list. It is other users
and developers who are answering your questions. They do it in their
free time and are not paid for doing it.

Contributing your own code
--------------------------

If you are planning to make an extension to or already have a piece of
your own code which could be useful to others, you are very welcome to
contribute it to the community. Before you start making any changes to
the code, you should obtain the current development version of it. For
more information about how to obtain the development version, refer to
the homepage http://espressomd.org.

It is also generally a good idea to contact the mailing lists before you
start major coding projects. It might be that someone else is already
working on the problem or has a solution at hand.

Developers’ guide
-----------------

Besides the User guide, also contains a Developers’ guide which is a
programmer documentation automatically built from comments in the source
code and using Doxygen. It provides a cross-referenced documentation of
all functions and data structures available in source code. It can be
built by typing

make doxygen

in the build directory. Afterwards it can be found in the subdirectory
of the build directory: ``doc/doxygen/html/index.html``.

A recent version of this guide can also be found on the homepage
http://espressomd.org/html/dox/.

User’s guide
------------

If, while reading this user guide, you notice any mistakes or badly (if
at all) described features or commands, you are very welcome to
contribute to the guide and have others benefit from your knowledge.

For this, you should also checkout the development version as described
on the homepage. You can then build the user guide by typing

make sphinx

