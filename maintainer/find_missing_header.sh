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
# check for missing GPL and copyright headers
#

files=`git ls-files --exclude-standard |
    egrep -v '\.(gz|data|dat|tab|chk|jpg|png|pdf|fig|gif|xcf|bib)$' |
    egrep -v '^testsuite/configs/|^old/' |
    egrep -v '(ChangeLog|AUTHORS|INSTALL|Doxyfile|latexmk.1|latexmkrc|bootstrap.sh)' |
    egrep -v '(\.gitignore|assemble_quickref.awk)'
    `

num=`echo $files | wc -w`

echo "Examining $num files."

# current_year=`date +%Y`

echo "Checking for missing copyright disclaimer"
echo "-----------------------------------------"
egrep -L "Copyright" $files
echo

# echo "Checking for missing GPL header"
# echo "--------------------------------"
# egrep -L "ESPResSo is free software" $files
# echo

#echo "Checking for missing Copyright line with current year $current_year:"
#egrep -rL "Copyright.*$current_year" *.[ch]
    
