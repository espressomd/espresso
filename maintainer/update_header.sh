# Copyright (C) 2012,2013 Olaf Lenz
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
# Add current year to copyright header
#
# 1. Find files that do have copyright disclaimers
# 1.1. Find files that do not have the current year in the disclaimer
# 1.1.1. Find files that have "The ESPResSo project" in the line
# 1.1.1.1. Append 2012 to the years
# 1.1.2. Insert new line "Copyright (C) <current_year> The ESPResSo
#        project" before the copyright disclaimer
#
files=`sh maintainer/files_with_header.sh`
num_files=`echo $files | wc -w`
current_year=`date +%Y`

echo "Examining $num_files files."

echo "Files with copyright disclaimer(s)..."
files=`egrep -l "Copyright" $files`
num_files=`echo $files | wc -w`
echo "  $num_files files."

echo "Files that are missing the current year ($current_year) in the copyright disclaimer(s)..."
files=`egrep -L "Copyright.*$current_year" $files`
for file in $files; do
    echo "  $file"
done

noyear_files=`egrep -l "Copyright.*The ESPResSo project" $files`
echo "  Adding current year to project copyright disclaimer..."
echo "    \"$current_year\""
tmpfile=`mktemp`
for file in $noyear_files; do
    echo "    $file"
    perl -pe "if (!\$done) { s/(Copyright\s+(?:\(C\)\s+)?\d+(?:\s*,\s*\d+)*)(\s+The ESPResSo project)/\1,$current_year\2/ and \$done=1; }" $file > $tmpfile
    cat $tmpfile > $file
done

noproject_files=`egrep -L "Copyright.*The ESPResSo project" $files`
echo "Files that are missing the project copyright disclaimer..."
num_files=`echo $noproject_files | wc -w`
echo "  $num_files files."
echo "  Adding project copyright disclaimer..."
disclaimer="Copyright (C) $current_year The ESPResSo project"
echo "    \"$disclaimer\""
for file in $noproject_files; do
    echo "    $file"
    perl -pe "if (!\$done) { s/^(.*)Copyright/\1$disclaimer\n\1Copyright/ and \$done=1; }" $file > $tmpfile
    cat $tmpfile > $file
done
rm $tmpfile

