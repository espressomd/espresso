#!/usr/bin/env sh
#
# Copyright (C) 2012-2022 The ESPResSo project
# Copyright (C) 2012 Olaf Lenz
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


# List all files that are eligible for a copyright header.

git ls-files --exclude-standard |
grep -vE '\.(blk|gz|npz|data|dat|tab|chk|jpg|png|pdf|fig|gif|xcf|css|bib|vtf|vtk|svg|ico|eps|rst|ipynb)$' |
grep -vE '^testsuite/configs/|^old/|^cmake/|^libs/' |
grep -vE '(ChangeLog|AUTHORS|COPYING|NEWS|README|INSTALL|README\.md|CONTRIBUTING\.md)' |
grep -vE '(\.gitmodules|\.github|\.gitignore|\.codecov\.yml|\.gitlab-ci\.yml|kodiak\.toml)' |
grep -vE '(\.clang-format|\.cmake-format|\.clang-tidy|\.coveragerc|\.pylintrc|ubsan\.supp)' |
grep -vE '(\.lgtm\.yml|\.pre-commit-config\.yaml|requirements\.txt)' |
grep -vE '(Doxyfile|latexmk\.1|latexmkrc|assemble_quickref\.awk)' |
grep -vE '(src/config/features\.def)' |
grep -vE '(featurelist)' |
grep -vE '(\.cproject|\.project|\.settings)'
