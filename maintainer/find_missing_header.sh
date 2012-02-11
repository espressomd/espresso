# Copyright (C) 2012 The ESPResSo project
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

files=`sh maintainer/files_with_header.sh`
num_files=`echo $files | wc -w`

echo "Examining $num_files files."

current_year=`date +%Y`
echo "Checking for copyright disclaimer missing the current year $current_year"
echo "-------------------------------------------------------------------"
egrep -L "Copyright.*$current_year" $files
echo

echo "Checking for missing copyright disclaimer"
echo "-----------------------------------------"
egrep -L "Copyright" $files
echo

echo "Checking for missing GPL/simple header"
echo "--------------------------------------"
no_gpl_files=`egrep -L "(ESPResSo|This program) is free software" $files`
egrep -L "Copying and distribution of this file" $no_gpl_files
echo

    
