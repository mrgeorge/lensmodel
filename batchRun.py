#! env python

# python script to send a bunch of jobs to the queue that each run driver.py for a different set of survey parameters
# Uses the qsub STDIN option to avoid writing lots of .qsub job files

import os

survey=["sdss","lsst"]
target=["galaxy","cluster"]
Rmin=[20.,40.]
magFrac=[0.5,1.0]
concPriorType=["true","low","high","free"]

queue="batch"

if((queue=="batch") | (queue=="fast")):
    nThreads=8
elif(queue=="big"):
    nThreads=48

for ss in range(len(survey)):
    for tt in range(len(target)):
        for rr in range(len(Rmin)):
            for mm in range(len(magFrac)):
                for cc in range(len(concPriorType)):
                    command="python $HOME/sdsslens/lensmodel/driver.py {} {} {} {} {} {}".format(survey[ss],target[tt],Rmin[rr],magFrac[mm],concPriorType[cc],nThreads)
                    jobName="j{}{}{}{}{}".format(ss,tt,rr,mm,cc)
                    os.system("echo {} | qsub -q {} -l nodes=1:ppn={} -j oe -m bea -M mgeorge@astro.berkeley.edu -V -N {}".format(command,queue,nThreads,jobName))
