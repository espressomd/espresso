#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    export EF_ALLOW_MALLOC_0=1; exec tclsh8.4 $0 $*;

set places [list pckr27 pckr28 pckr33 pckrls pckr37 pckr38 pckr39 pckr40 pckr41 pckr42 pckr43 pckr44 pckr45 pckr46 pckr47 pckr48 pckr49 pckr50]
set user [exec whoami]

exec echo "renice +19 `ps U $user | grep lamd | awk '{print \$1}'`" > renice.sh
exec echo "renice +19 `ps U $user | grep Espresso | awk '{print \$1}'`" >> renice.sh
exec chmod 755 renice.sh
exec mv -f renice.sh /people/thnfs/homes/$user/

puts "\nHello $user, thanks for re-nicing your jobs!"
foreach place $places {
    puts "Entering '$place' to renice..."
    catch { exec ssh $place renice.sh }
}
exec rm -f /people/thnfs/homes/$user/renice.sh
puts "All done - Have a (re-)nice(d) day!\n"

