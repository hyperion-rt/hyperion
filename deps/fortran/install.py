#!/usr/bin/env python

from __future__ import print_function

import os
import platform
import tempfile
import sys
import tarfile
import hashlib
import shlex
from distutils import version
import subprocess

try:
    from urllib import urlopen
except ImportError:
    from urllib.request import urlopen

try:
    from urllib import urlencode
except ImportError:
    from urllib.parse import urlencode

TEST_FC_F90 = '''
program test
  print *,'hello, world!'
end program test
'''

TEST_MPIF90_F90 = '''
program test

  use mpi

  implicit none

  integer :: ierr, rank, nproc

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, nproc, ierr)
  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_finalize(ierr)

  print *, 'MPI Installation Successful!'

end program test
'''

TEST_HDF5_F90 = '''
program test

  use hdf5
  use h5tb

  print *, 'HDF5 Installation Successful!'

end program test
'''

TEST_ARCH_C = "int main() {}"

TEST_ARCH_F90 = '''
program test
end program test
'''

INSTALL_HDF5 = True
INSTALL_MPICH2 = True

VERBOSE = False

fc = None
cc = None
cxx = None

prefix = None


def usage():
    print("")
    print(" Usage: python install.py [options] installation_path")
    print("")
    print(" Available options include:")
    print("")
    print("  --only-hdf5             only compile the HDF5 library")
    print("  --only-mpich2           only compile the MPICH2 library")
    print("  --fc-compiler=<value>   force script to use a specific fortran compiler")
    print("  --cc-compiler=<value>   force script to use a specific C compiler")
    print("  --cxx-compiler=<value>  force script to use a specific C++ compiler")
    print("  --verbose               output standard output to terminal")
    print("  --help                  display this message")
    print("")
    sys.exit(1)

if '--help' in sys.argv[1:]:
    usage()

for arg in sys.argv[1:]:

    if arg == '--only-hdf5':
        INSTALL_MPICH2 = False

    if arg == '--only-mpich2':
        INSTALL_HDF5 = False

    if arg == '--verbose':
        VERBOSE = True

    if arg.startswith('--fc-compiler'):
        if '=' in arg:
            fc = arg.split('=')[1].replace('"', '').replace("'", "")
        else:
            usage()

    if arg.startswith('--cc-compiler'):
        if '=' in arg:
            cc = arg.split('=')[1].replace('"', '').replace("'", "")
        else:
            usage()

    if arg.startswith('--cxx-compiler'):
        if '=' in arg:
            cxx = arg.split('=')[1].replace('"', '').replace("'", "")
        else:
            usage()

    if not arg.startswith('-'):
        if prefix is None:
            prefix = arg
        else:
            print("ERROR: only one non-optional argument can be provided (the installation path)")
            usage()

if prefix is None:
    print("ERROR: please specify installation directory with 'python install.py directory'")
    sys.exit(1)


def run(command, logfile):

    import subprocess

    if VERBOSE:
        status = subprocess.call(command + ' 2>&1 | tee ' + logfile, shell=True, executable="/bin/bash")
    else:
        status = subprocess.call(command + ' >& ' + logfile, shell=True, executable="/bin/bash")

    if status != 0:
        par = {}
        par['api_dev_key'] = 'd3b3e1a0b666fbcbe383162be949f81e'
        par['api_option'] = 'paste'
        par['api_paste_code'] = open(logfile).read()
        par['api_paste_name'] = command
        par['api_paste_format'] = 'bash'
        u = urlopen('http://pastebin.com/api/api_post.php',
                           data=urlencode(par))
        url = u.read()
        print("=" * 72)
        print("The installation failed. The log of the failed command has been sent")
        print("to: " + url)
        print("=" * 72)
        sys.exit(1)
    else:
        return

print("=" * 72)
print("Determining system setup")
print("=" * 72)

start_dir = os.path.abspath('.')

prefix = os.path.abspath(prefix)

ZLIB_URL = "http://downloads.sourceforge.net/project/libpng/zlib/1.2.7/zlib-1.2.7.tar.gz"
ZLIB_SHA1 = '4aa358a95d1e5774603e6fa149c926a80df43559'

HDF5_URL = 'http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.9/src/hdf5-1.8.9.tar.gz'
HDF5_SHA1 = '4ba3ede947b1571e9586fdeb8351d6585a56133c'

MPICH2_URL = 'http://www.mpich.org/static/tarballs/1.5/mpich2-1.5.tar.gz'
MPICH2_SHA1 = 'be7448227dde5badf3d6ebc0c152b200998421e0'

print(" -> changing to work directory... ")
work_dir = tempfile.mkdtemp()
os.chdir(work_dir)
print("    %s" % work_dir)

print(" -> determining platform... ", end=' ')

# Determine system
system = platform.uname()[0]

# Determine system version
if system == 'Darwin':

    # Determine MacOS X version
    system_version = "MacOS " + platform.mac_ver()[0]

else:

    system_version = ''

print(system, '/', system_version)

if fc is None:

    print(" -> determining best fortran compiler... ", end=' ')

    # Determine available fortran compilers
    FORTRAN_COMPILERS = ['ifort',
                         'gfortran',
                         'gfortran-mp-4.3',
                         'gfortran-mp-4.4',
                         'gfortran-mp-4.5',
                         'gfortran-mp-4.6',
                         'g95',
                         'pgfortran',
                         'pgf95']

    # Create test script
    open('test_fc.f90', 'w').write(TEST_FC_F90)

    fc = None
    for compiler in FORTRAN_COMPILERS:

        if system_version.startswith('MacOS 10.5') and compiler == 'ifort':
            compiler += " -m32"

        return_code = subprocess.call(compiler + ' test_fc.f90 -o test_fc', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if return_code == 0:
            fc = compiler
            break

    print(fc)

else:

    print(" -> using user-specified fortran compiler...", fc)

if fc is None:
    print("ERROR: none of the known Fortran compilers were found")
    sys.exit(1)

if cc is None:

    print(" -> determining best C compiler...", end=' ')

    # TODO: implement search
    cc = 'gcc'

    print(cc)

else:

    print(" -> using user-specified C compiler...", cc)

if cc is None:
    print("ERROR: none of the known C compilers were found")
    sys.exit(1)

if cxx is None:

    print(" -> determining best C compiler...", end=' ')

    # TODO: implement search
    cxx = 'g++'

    print(cxx)

else:

    print(" -> using user-specified C++ compiler...", cxx)

if cxx is None:
    print("ERROR: none of the known C++ compilers were found")
    sys.exit(1)

# The following section deals with issues that occur when using the intel
# fortran compiler with gcc.

# Check whether C compiler is gcc
p = subprocess.Popen(shlex.split(cc + ' --version'), stdout=subprocess.PIPE)
output = (p.communicate()[0]).decode('ascii').strip().splitlines()[0]
is_gcc = 'GCC' in output

# Check whether Fortran compiler is ifort
p = subprocess.Popen(shlex.split(fc + ' --version'), stdout=subprocess.PIPE)
output = (p.communicate()[0]).decode('ascii').strip().splitlines()[0]
is_ifort = '(IFORT)' in output

# Check whether Fortran compiler is g95
p = subprocess.Popen(shlex.split(fc + ' --version'), stdout=subprocess.PIPE)
output = (p.communicate()[0]).decode('ascii').strip().splitlines()[0]
is_g95 = 'g95' in output

# Check whether Fortran compiler is pgfortran
p = subprocess.Popen(shlex.split(fc + ' --version'), stdout=subprocess.PIPE)
output = (p.communicate()[0]).decode('ascii').strip().splitlines()[0]
is_pgfortran = 'pgfortran' in output or \
               'pgf95' in output

# On MacOS X, when using gcc 4.5 or 4.6, the fortran compiler needs to link
# with libgcc_eh that is in the gcc library directory. This is not needed if
# using gfortran 4.5 or 4.6, but it's easier to just add it for all
# compilers.
if system == 'Darwin' and is_gcc:

    # Determine gcc version
    p = subprocess.Popen(shlex.split(cc + ' -dumpversion'), stdout=subprocess.PIPE)
    gcc_version = version.LooseVersion((p.communicate()[0]).decode('ascii').splitlines()[0])

    if gcc_version >= version.LooseVersion('4.5.0'):

        p = subprocess.Popen(shlex.split(cc + ' -print-search-dirs'),
                             stdout=subprocess.PIPE)
        output = (p.communicate()[0]).decode('ascii').splitlines()[0]
        if output.startswith('install:'):
            libs = output.split(':', 1)[1].strip()
            libs = ' -L' + libs.replace(':', '-L')
            libs += ' -lgcc_eh'
            fc += libs
            print(" -> SPECIAL CASE: adjusting fortran compiler:", fc)
        else:
            print("ERROR: unexpected output for %s -print-search-dirs: %s" % (cc, output))
            sys.exit(1)

    # Check whether the C and Fotran compiler give different architecture builds by default

    open('test_arch.c', 'w').write(TEST_ARCH_C)
    subprocess.Popen(shlex.split(cc + ' test_arch.c -o test_arch_c')).wait()
    p = subprocess.Popen(shlex.split('file test_arch_c'), stdout=subprocess.PIPE)
    output = (p.communicate()[0]).decode('ascii').splitlines()[0].strip()
    if output == 'test_arch_c: Mach-O 64-bit executable x86_64':
        arch_c = 64
    elif output == 'test_arch_c: Mach-O executable i386':
        arch_c = 32
    else:
        arch_c = None

    open('test_arch.f90', 'w').write(TEST_ARCH_F90)
    subprocess.Popen(shlex.split(fc + ' test_arch.f90 -o test_arch_f90')).wait()
    p = subprocess.Popen(shlex.split('file test_arch_f90'), stdout=subprocess.PIPE)
    output = (p.communicate()[0]).decode('ascii').splitlines()[0].strip()
    if output == 'test_arch_f90: Mach-O 64-bit executable x86_64':
        arch_f90 = 64
    elif output == 'test_arch_f90: Mach-O executable i386':
        arch_f90 = 32
    else:
        arch_f90 = None

    if arch_c is None or arch_f90 is None:
        pass  # just be safe and don't assume anything
    elif arch_c != arch_f90:
        if arch_c == 32:
            cc += '- m64'
            cxx += ' -m64'
            print(" -> SPECIAL CASE: adjusting C compiler:", cc)
        else:
            fc += ' -m64'
            print(" -> SPECIAL CASE: adjusting fortran compiler:", fc)

if INSTALL_HDF5:

    print("=" * 72)
    print("Installing ZLIB (for HDF5)")
    print("=" * 72)

    zlib_file = os.path.basename(ZLIB_URL)
    if os.path.exists(zlib_file):
        sha1 = hashlib.sha1(open(zlib_file, 'rb').read()).hexdigest()
        if sha1 == ZLIB_SHA1:
            print(" -> file exists, skipping download")
        else:
            print(" -> file exists but incorrect SHA1, re-downloading")
            open(zlib_file, 'wb').write(urlopen(ZLIB_URL).read())
    else:
        print(" -> downloading")
        open(zlib_file, 'wb').write(urlopen(ZLIB_URL).read())

    print(" -> expanding tarfile")
    t = tarfile.open(zlib_file, 'r:gz')
    t.extractall()
    os.chdir(zlib_file.replace('.tar.gz', ''))

    print(" -> configuring")
    run('./configure --prefix={prefix}'.format(prefix=prefix), 'log_configure')

    print(" -> making")
    run('make', 'log_make')

    print(" -> installing")
    run('make install', 'log_make_install')

    os.chdir(work_dir)

    print("=" * 72)
    print("Installing HDF5")
    print("=" * 72)

    hdf5_file = os.path.basename(HDF5_URL)
    if os.path.exists(hdf5_file):
        sha1 = hashlib.sha1(open(hdf5_file, 'rb').read()).hexdigest()
        if sha1 == HDF5_SHA1:
            print(" -> file exists, skipping download")
        else:
            print(" -> file exists but incorrect SHA1, re-downloading")
            open(hdf5_file, 'wb').write(urlopen(HDF5_URL).read())
    else:
        print(" -> downloading")
        open(hdf5_file, 'wb').write(urlopen(HDF5_URL).read())

    print(" -> expanding tarfile")
    t = tarfile.open(hdf5_file, 'r:gz')
    t.extractall()
    os.chdir(hdf5_file.replace('.tar.gz', ''))

    # SPECIAL CASE - g95 requires patching
    if is_g95:
        print(" -> SPECIAL CASE: patching for g95")
        conf = open('config/gnu-fflags', 'rb').read()
        conf = conf.replace('-Wconversion -Wunderflow ', '')
        open('config/gnu-fflags', 'w').write(conf)

    print(" -> configuring")
    run('./configure FC="{fc}" --enable-fortran --enable-hl --with-zlib={prefix}/include,{prefix}/lib --prefix={prefix}'.format(fc=fc, prefix=prefix), 'log_configure')

    print(" -> making")
    run('make', 'log_make')

    print(" -> installing")
    run('make install', 'log_make_install')

    os.chdir(work_dir)

    print(" -> testing installation")
    open('test_hdf5.f90', 'w').write(TEST_HDF5_F90)
    run('{prefix}/bin/h5fc test_hdf5.f90 -o test_hdf5'.format(prefix=prefix), 'log_test_hdf5_compile')
    run('./test_hdf5', 'log_test_hdf5_run')

if INSTALL_MPICH2:

    print("=" * 72)
    print("Installing MPICH2")
    print("=" * 72)

    mpich2_file = os.path.basename(MPICH2_URL)
    if os.path.exists(mpich2_file):
        sha1 = hashlib.sha1(open(mpich2_file, 'rb').read()).hexdigest()
        if sha1 == MPICH2_SHA1:
            print(" -> file exists, skipping download")
        else:
            print(" -> file exists but incorrect SHA1, re-downloading")
            open(mpich2_file, 'wb').write(urlopen(MPICH2_URL).read())
    else:
        print(" -> downloading")
        open(mpich2_file, 'wb').write(urlopen(MPICH2_URL).read())

    print(" -> expanding tarfile")
    t = tarfile.open(mpich2_file, 'r:gz')
    t.extractall()

    print(" -> configuring")
    os.chdir(mpich2_file.replace('.tar.gz', ''))
    run('./configure F77="{fc}" FC="{fc}" CC="{cc}" CXX="{cxx}" --enable-fc --prefix={prefix}'.format(fc=fc, cc=cc, cxx=cxx, prefix=prefix), 'log_configure')

    print(" -> making")
    run('make', 'log_make')

    print(" -> installing")
    run('make install', 'log_make_install')

    os.chdir(work_dir)

    print(" -> testing installation")
    open('test_mpif90.f90', 'w').write(TEST_MPIF90_F90)
    run('{prefix}/bin/mpif90 test_mpif90.f90 -o test_mpif90'.format(prefix=prefix), 'log_test_mpif90_compile')
    run('./test_mpif90', 'log_test_mpif90_run')

print("=" * 72)

# Go back to starting directory
os.chdir(start_dir)

# Print out message regarding environment variables
print("Installation succeeded! You will now need to add the installation\n" \
    + "directory to your $PATH environment variable. For example, in bash,\n" \
    + "you should add the following to your ~/.bash_profile:\n" \
    + "\n" \
    + "export PATH=%s/bin:$PATH\n" % prefix \
    + "\n" \
    + "Check that your installation is being correctly picked up:" \
    + "\n" \
    + "$ which h5fc\n" \
    + "%s/bin/h5fc\n" % prefix\
    + "$ which mpif90\n" \
    + "%s/bin/mpif90" % prefix)

print("=" * 72)
