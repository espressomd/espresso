#!/bin/sh

sed -i '/ostringstream msg;/d' $1
sed -i 's/msg <</runtimeErrorMsg() <</g' $1
sed -i '/runtimeError(msg);/d' $1

