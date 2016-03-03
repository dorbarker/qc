from getpass import getuser
from socket import gethostname
import subprocess
import functools
import os

class Git(object):
    
    def __init__(self, args):

        self.arguments = self._format_message(str(args))

    def run(self):

        if self._is_repo():
            self._update()
        else:
            self._create()

    def _subproc(self, cmd):

        #exitcode = subprocess.call(cmd, stdout = subprocess.DEVNULL,
        #                           stderr = subprocess.DEVNULL)

        exitcode = subprocess.call(cmd)
        return exitcode

    def _is_repo(self):
        
        return self._subproc(['git', 'status']) == 0

    def _create(self):

        user = getuser()
        email = '{}@{}'.format(user, gethostname())

        self._subproc(['git', 'init'])
        self._subproc(['git', 'config', 'user.name', user])
        self._subproc(['git', 'config', 'user.email', email])
        self._subproc(['git', 'add', '-A'])
        self._subproc(['git', 'commit', '-m', '"Original unfiltered"'])

    def _update(self):

        self._subproc(['git', 'add', '-u'])
        self._subproc(['git', 'commit', '-m', '"{}"'.format(self.arguments)])
    
    def _format_message(self, msg):

        o = msg.replace('Namespace(', '').replace("'",'').replace('"','')[:-1]
        return o

def version(git_dir):
    def decorate(f):
        @functools.wraps(f)
        def decorator(args):
            original_dir = os.path.abspath(os.curdir)
            os.chdir(git_dir)
            f(args)
            os.chdir(original_dir)
        return decorator
    return decorate

def commit(git_dir, args):
    
    @version(git_dir)
    def commit_to_repository(args):

        git = Git(args)
        git.run()

    commit_to_repository(args)

