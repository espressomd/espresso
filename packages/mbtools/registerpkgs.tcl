# 
# tcl script to create pkgIndex for all packages inside the packages directory
#
# If you add a new package you should ensure that the package gets its
# own directory and that you add the name of that directory to the
# list "pkgs" below. 
#
# Inside the directory you should always have one main tcl file (with
# same name as directory) that contains the package provide command
# and gives the version.  The minor version number should represent
# bug fixes.  The middle version should represent new features and the
# major version should be changed if there are compatibility issues
# between versions. 
#


set pkgs { system_generation analysis utils }

foreach package $pkgs {
    pkg_mkIndex -verbose $package 
}


