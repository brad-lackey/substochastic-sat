runsolver Copyright (C) 2010-2013 Olivier ROUSSEL

This is runsolver version 3.3.5 (svn: 2013)

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

command line: ./runsolver -d 3 -C cpu=-4 -W -4 -v out.v -w out.w ./ssmc 

running on 32 cores: 0-7,16-23,8-15,24-31

Enforcing CPUTime limit (soft limit, will send SIGTERM then SIGKILL): 0 seconds
Enforcing CPUTime limit (hard limit, will send SIGXCPU): 30 seconds
Enforcing wall clock limit (soft limit, will send SIGTERM then SIGKILL): 18446744073709551612 seconds
Current StackSize limit: 8192 KiB


[startup+0 s]
/proc/loadavg: 0.01 0.02 1.21 3/1208 23756
/proc/meminfo: memFree=123770928/131984384 swapFree=10816/1652244
[pid=23756] ppid=23755 vsize=0 CPUtime=0 cores=0-31
/proc/23756/stat : 23756 (ssmc) Z 23755 23756 21533 34833 23754 4227084 91 0 0 0 0 0 0 0 20 0 1 0 14703973 0 0 18446744073709551615 0 0 0 0 0 0 0 0 16384 1 0 0 17 23 0 0 0 0 0 0 0 0 0 0 0 0 512
/proc/23756/statm: 0 0 0 0 0 0 0

Solver just ended. Dumping a history of the last processes samples

Child status: 2
Real time (s): 0.023665
CPU time (s): 0
CPU user time (s): 0
CPU system time (s): 0
CPU usage (%): 0
Max. virtual memory (cumulated for all children) (KiB): 0

getrusage(RUSAGE_CHILDREN,...) data:
user time used= 0
system time used= 0
maximum resident set size= 1648
integral shared memory size= 0
integral unshared data size= 0
integral unshared stack size= 0
page reclaims= 91
page faults= 0
swaps= 0
block input operations= 0
block output operations= 0
messages sent= 0
messages received= 0
signals received= 0
voluntary context switches= 1
involuntary context switches= 1

runsolver used 0.004 second user time and 0.02 second system time

The end
