# Release checklist

## May change repository

   * Update CI containers to recent version of libs and distributions
   * Check Python modules version in [`/requirements.txt`](python/requirements.txt)
   * Check CMake version in [`/CMakeLists.txt`](python/CMakeLists.txt)
   * Check contributors in [`/AUTHORS`](python/AUTHORS)
      + relevant script: [`find_potentially_missing_authors`](python/maintainer/find_potentially_missing_authors)
   * Check tutorials and samples
   * Spell and grammar check comments and docs
      + relevant thread: #2216
   * Check and fix copyright headers
      + relevant threads: #2198, #1421
   * Test installation (`make check_cmake_install`)
   * Test builds on i386 and some big-endian architectures
   * Make and test Tarball

## No changes after this

   * Tag release
   * Fork to version branch
   * Build and deploy static Sphinx documentation to [espressomd.org/html/doc](http://espressomd.org/html/doc/)
