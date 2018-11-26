# Release checklist

## May change repository

   * Update CI containers to recent version of libs and distributions
   * Check Python modules version in [`/requirements.txt`](python/requirements.txt)
   * Check tutorials and samples
   * Check and fix copyright headers
   * Check contributors in [`/AUTHORS`](python/AUTHORS)
   * Check CMake version in [`/CMakeLists.txt`](python/CMakeLists.txt)
   * Spell and grammar check comments and docs
   * Make and test Tarball
   * Test installation (`make check_cmake_install`)

## No changes after this

   * Tag release
   * Fork to version branch
   * Build and deploy static Sphinx documentation to [espressomd.org/html/doc](http://espressomd.org/html/doc/)
