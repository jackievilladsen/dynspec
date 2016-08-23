def load_dirs():
    # load list of directories from /data/jrv/casa_utils/dynspec/P_dir.txt
    fname = '/data/jrv/casa_utils/dynspec/P_dir.txt' 
    f = open(fname)
    lines = f.readlines()
    f.close()
    dirs = []
    for l in lines:
        if l[0]!='#': # allows commenting with '#'
            dirs.append(l.rstrip())
    return dirs

dirs = load_dirs()
for d in dirs:
    os.chdir(d)
    print 'Running flagPband.py in directory', d
    execfile('/data/jrv/casa_utils/dynspec/flagPband.py')
