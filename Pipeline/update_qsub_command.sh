rsync -Cuavz qsub_command.py ~/scripts/qsub_command.py
rsync -Cuavz qsub_command.py daveminh@login.osgconnect.net:~/scripts/qsub_command.py
rsync -Cuavz qsub_command.py dminh@ccb.tbc.iit.edu:~/scripts/qsub_command.py
ssh dminh@ccb.tbc.iit.edu "rsync -Cuavz ~/scripts/qsub_command.py dminh@otsgrid.iit.edu:~/scripts/qsub_command.py"
rsync -Cuavz qsub_command.py dminh@bridges.psc.xsede.org:~/scripts/qsub_command.py
rsync -Cuavz qsub_command.py daveminh@comet.sdsc.xsede.org:~/scripts/qsub_command.py
# This one needs a password and can be cancelled sometimes
rsync -Cuavz qsub_command.py minh1@syrah.llnl.gov:/g/g19/minh1/scripts/qsub_command.py
