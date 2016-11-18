import os, sys, glob, commands, subprocess, datetime

topdir = "/home/server/server2/projects/Gromacs"
resultsdir = '/home/server/server2/.results'
now = datetime.datetime.now()
file = str(now.month) + "-" + str(now.day) + "-" + str(now.year) + ".dat"
f = open(os.path.join(resultsdir, file), 'w+')

all_ns = 0
projects = []

for i in os.listdir(topdir):
    if i.startswith("p") and os.path.isdir(os.path.join(topdir, i)):
        projects.append(i)

for i in projects:
    os.chdir(os.path.join(topdir, i))
    mdp=[]
    for file in os.listdir("."):
        if file.endswith(".mdp") and file != "mdout.mdp" and file != "mdpOut.mdp":
            mdp.append(file)
#    if len(mdp) > 1:
#        print "project", i, "has too many mdp files to be parsed."
#        print "Using", mdp[0]
    if len(mdp) == 0:
#        print "There is no mdp file in", os.path.join(topdir, i)
        continue
    mdp = mdp[0]
    try: # extracting number of steps from mdp file
        cmd = 'cat ' + mdp + ' | grep nsteps | grep -o "[0-9].*" | sed "s/ *;//" | cut -d" " -f1'
        ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        steps = int(ps.communicate()[0])
    except ValueError, e:
        print("Cannot parse mdp file in", os.path.join(topdir, i), "for number of steps.")
    
    ns_per_frame = steps * 0.002 / 1000
    Continue = True 
    lengths = []

    proj_number = str(filter(str.isdigit, i))
    projdir = "/home/server/server2/data/SVR166219/PROJ" + proj_number
    if os.path.exists(projdir):
#    print 'frame\tlength(ns)\tN(frame)'
        frame = 0
        maxframe = 1000
        n = int(commands.getoutput('find %s | grep xtc | wc -l' %projdir ))
        if n > 0:
            total_ns = ns_per_frame*n
            f.write(str(proj_number) + "\t%3.1f\n" %total_ns)
            all_ns += total_ns

f.write("total\t" + str(all_ns))


