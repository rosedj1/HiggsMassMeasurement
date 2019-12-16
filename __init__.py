import sys, os

def addCWDtoPyPath():
    '''
    PURPOSE:
        Solve the import issue with Python.
        If CWD is not in PYTHON PATH, add it and all parent directories of CWD.
    '''
    if os.getcwd() not in sys.path:
        print "Adding CWD and all parent directories to PYTHONPATH."

        pwd = os.getcwd()

        for count in range( len(pwd.split('/'))-1 ):
            if count == 0:
                sys.path.append(pwd)
            else:
                pwd,_ = os.path.split(pwd)
                sys.path.append(pwd)
