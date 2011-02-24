#!/bin/sh

cd $ESPRESSO_SOURCE
cd XCodeEspresso
mkdir XCodeEspresso.xcodeproj


sed -e "s#<ESPRESSO_SOURCE>#$ESPRESSO_SOURCE#g" XCodeGen/project.pbxproj > XCodeEspresso.xcodeproj/project.pbxproj
sed -e "s#<ESPRESSO_SOURCE>#$ESPRESSO_SOURCE#g" XCodeGen/user.mode1 > XCodeEspresso.xcodeproj/`whoami`.mode1
sed -e "s#<ESPRESSO_SOURCE>#$ESPRESSO_SOURCE#g" XCodeGen/user.pbxuser > XCodeEspresso.xcodeproj/`whoami`.pbxuser
