## Automated bug sweep

* This folder contains a Claude Code workflow that sweeps the entire src/ for bugs. 
* For each find, 2 verifiers try to refute it before it is reported
* The code is scanned in chunks by parallel agents.
* Separate agents validate the call path from Python to the core
* The outcome is a prioritized list of validated findes.
* Calling the workflow is expensive (~250 individual sub-agents). As of June 2026, it takes approximately 1 to 1.5 session limits on the extended seat in the team plan.


