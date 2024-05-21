import subprocess
import sys
import os
import tempfile
import shutil
import site
try:
    import git
except ModuleNotFoundError:
    print("CHIRIPA cannot be installed due to git library is not accesible")
    print("Please install GitPython library using pip install GitPython")
    sys.exit()

from setuptools import setup, find_packages, Extension
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler
from Cython.Build import cythonize

# Install packages from pip ==============================================================
def install_with_pip(pack, vers = None):

    # sys.executable gives the path to the python interpreter
    if vers is None:
        print("CHIRIPA: Installing {}".format(pack))
        subprocess.call([sys.executable, "-m", "pip", "install", "{0}".format(pack)])
    else:
        print("CHIRIPA: Installing {}=={}".format(pack, vers))
        subprocess.call([sys.executable, "-m", "pip", "install", "{0}=={1}".format(pack, vers)])

# Disabling-output-when-compiling-with-distutil =================================================
def hasfunction(cc, funcname, include=None, extra_postargs=None):
    # From http://stackoverflow.com/questions/
    #            7018879/disabling-output-when-compiling-with-distutils
    tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
    devnull = oldstderr = None
    try:
        try:
            fname = os.path.join(tmpdir, 'funcname.c')
            with open(fname, 'w') as fout:
                if include is not None:
                    fout.write('#include {0!s}\n'.format(include))
                fout.write('int main(void) {\n')
                fout.write('    {0!s};\n'.format(funcname))
                fout.write('}\n')
            # Redirect stderr to /dev/null to hide any error messages
            # from the compiler.
            # This will have to be changed if we ever have to check
            # for a function on Windows.
            devnull = open('/dev/null', 'w')
            oldstderr = os.dup(sys.stderr.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())
            objects = cc.compile([fname], output_dir=tmpdir,
                                 extra_postargs=extra_postargs)
            cc.link_executable(objects, os.path.join(tmpdir, "a.out"))
        except Exception:
            return False
        return True
    finally:
        if oldstderr is not None:
            os.dup2(oldstderr, sys.stderr.fileno())
        if devnull is not None:
            devnull.close()
        shutil.rmtree(tmpdir)

# Does this compiler support OpenMP parallelization?""" ==============================================================
def detect_openmp():

    print("CHIRIPA: Attempting to autodetect OpenMP support... ", end="")
    compiler = new_compiler()
    customize_compiler(compiler)
    compiler.add_library('gomp')
    include = '<omp.h>'
    extra_postargs = ['-fopenmp']
    hasopenmp = hasfunction(compiler, 'omp_get_num_threads()', include=include,
                            extra_postargs=extra_postargs)
    if hasopenmp:
        print("CHIRIPA: Compiler supports OpenMP")
    else:
        print("CHIRIPA: Did not detect OpenMP support.")

    return hasopenmp

# Install indigox-bond software ======================================================================================
def install_indigox_bond():

    """
    Installing the indigo-bond library if is not present in the python enviorement.
    """

    giturl = 'https://github.com/allison-group/indigo-bondorder.git'
    install_dir = 'thirdparty/indigo-bondorder'

    try:
        import indigox as ix
        print("CHIRIPA: indigo-bondorder is already installed in your system. {}".format(giturl))
    except ModuleNotFoundError:
        print("CHIRIPA: indigo-bondorder is not installed in your system")
        print("CHIRIPA: Installing from git... {}".format(giturl))

        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            print("================= ERROR INSTALL ================")
            print("CHIRIPA: Cannot find CMake executable")
            print("CHIRIPA: The installation is aborted")
            print("================= ERROR INSTALL ================")
            sys.exit()


        # Look at thirdparty directory
        if os.path.isdir("thirdparty"):
            pass
        else:
            os.makedirs("thirdparty")

        fullpath_cmake = os.path.abspath(install_dir)

        # Check if exists a distribution of indigox in the thirdparty directory
        if os.path.isdir("thirdparty/indigo-bondorder"):
            pass
        else:
            #git clone https://github.com/allison-group/indigo-bondorder.git
            try:
                git.Repo.clone_from(giturl, install_dir)
            except git.GitCommandError:
                if not os.path.isdir(install_dir):
                    print("================= ERROR INSTALL ================")
                    print("CHIRIPA: The github repository for indigo-bondorder is not valid or not exists.!!!")
                    print("CHIRIPA: giturl     : {}".format(giturl))
                    print("CHIRIPA: install_dir: {}".format(install_dir))
                    print("CHIRIPA: Indigo-bondorder cannot be installed")
                    print("CHIRIPA: The installation is aborted")
                    print("================= ERROR INSTALL ================")
                    sys.exit()
                else:
                    pass

            subprocess.call(["rm", "-rf", "thirdparty/indigo-bondorder/build"])
            subprocess.call(["mkdir", "thirdparty/indigo-bondorder/build"])
            os.chdir("thirdparty/indigo-bondorder/build")
            cmake_arguments = [ "-DCMAKE_INSTALL_PREFIX={}".format(fullpath_cmake)]
            subprocess.check_call(["cmake", "{}".format(fullpath_cmake)]+cmake_arguments)
            subprocess.call("make")
            subprocess.call(["make", "install"])
            os.chdir("../../")
            subprocess.call(["tar","cvfz", "indigo-bondorder.tar.gz", "indigo-bondorder"])
            subprocess.call(["rm", "-rf", "indigo-bondorder ./git"])
            os.chdir("..")

        print("The *.so library has been installed in: {envdir}/lib/python3.8/site-packages/indigox/pyindigox.cpython-36m-x86_64-linux-gnu.so")
        print("                                        {envdir}/lib/python3.8/site-packages/indigox/__init__.py")

# Install eigen library software ======================================================================================
def install_eigen():

    """
    Installing the eigen library which is needed in openbabel.
    """

    import git

    giturl = 'https://gitlab.com/libeigen/eigen.git'
    install_dir = 'thirdparty/eigen'

    if not os.path.isdir("thirdparty/eigen/eigen_library/include"):
        print("CHIRIPA: eigen is not installed in your system")
        print("http://eigen.tuxfamily.org/index.php?title=Main_Page")
        print("CHIRIPA: Installing from git... {}".format(giturl))


    try:
        subprocess.check_output(['cmake', '--version'])
    except OSError:
        print("================= ERROR INSTALL ================")
        print("CHIRIPA: Cannot find CMake executable")
        print("CHIRIPA: The installation is aborted")
        print("================= ERROR INSTALL ================")
        sys.exit()

    # Look at thirdparty directory
    if os.path.isdir("thirdparty"):
        pass
    else:
        os.makedirs("thirdparty")

    fullpath_cmake = os.path.abspath(install_dir)

    # Check if exists a distribution of indigox in the thirdparty directory
    if os.path.isdir("thirdparty/eigen/eigen_library/include"):
        pass
    else:
        #git clone https://gitlab.com/libeigen/eigen.git
        try:
            git.Repo.clone_from(giturl, install_dir)
        except git.GitCommandError:
            if not os.path.isdir(install_dir):
                print("================= ERROR INSTALL ================")
                print("CHIRIPA: The github repository for openbabel is not valid or not exists.!!!")
                print("CHIRIPA: giturl     : {}".format(giturl))
                print("CHIRIPA: install_dir: {}".format(install_dir))
                print("CHIRIPA: openbabel cannot be installed")
                print("CHIRIPA: The installation is aborted")
                print("================= ERROR INSTALL ================")
                sys.exit()
            else:
                pass

        subprocess.call(["rm", "-rf", "thirdparty/eigen/build"])
        subprocess.call(["mkdir", "thirdparty/eigen/build"])
        subprocess.call(["mkdir", "thirdparty/eigen/eigen_library"])
        os.chdir("thirdparty/eigen/build")
        cmake_arguments1 = [ "-DCMAKE_INSTALL_PREFIX={}".format(fullpath_cmake+"/eigen_library")]
        subprocess.check_call(["cmake", "{}", "{}".format(fullpath_cmake)]+cmake_arguments1)
        subprocess.call("make")
        subprocess.call(["make", "install"])
        os.chdir("../../")
        subprocess.call(["tar","cvfz", "eigen.tar.gz", "eigen"])
        os.chdir("..")

        print("The eigen library has been installed in: thirdparty/eigen/eigen_library")

# Install indigox-bond software ======================================================================================
def install_openbabel():

    """
    Installing the openbabel library if is not present in the python enviorement.
    """
    #
    # giturl = 'https://github.com/openbabel/openbabel.git'
    # install_dir = 'thirdparty/openbabel'
    #
    # try:
    #     import openbabel
    #     print("CHIRIPA: openbabel is already installed in your system. {}".format(giturl))
    # except ModuleNotFoundError:
    #     print("CHIRIPA: openbabel is not installed in your system")
    #     print("CHIRIPA: Installing from git... {}".format(giturl))
    #
    #     try:
    #         subprocess.check_output(['cmake', '--version'])
    #     except OSError:
    #         print("================= ERROR INSTALL ================")
    #         print("CHIRIPA: Cannot find CMake executable")
    #         print("CHIRIPA: The installation is aborted")
    #         print("================= ERROR INSTALL ================")
    #         sys.exit()
    #
    #     # Look at thirdparty directory
    #     if os.path.isdir("thirdparty"):
    #         pass
    #     else:
    #         os.makedirs("thirdparty")
    #
    #     fullpath_cmake = os.path.abspath(install_dir)
    #
    #     # Check if exists a distribution of indigox in the thirdparty directory
    #     #cJ if os.path.isdir("thirdparty/openbabel"):
    #     if os.path.isdir("thirdparty/openbabel"):
    #          pass
    #     else:
    #         #git clone https://github.com/openbabel/openbabel.git
    #         try:
    #             git.Repo.clone_from(giturl, install_dir)
    #         except git.GitCommandError:
    #             if not os.path.isdir(install_dir):
    #                 print("================= ERROR INSTALL ================")
    #                 print("CHIRIPA: The github repository for openbabel is not valid or not exists.!!!")
    #                 print("CHIRIPA: giturl     : {}".format(giturl))
    #                 print("CHIRIPA: install_dir: {}".format(install_dir))
    #                 print("CHIRIPA: openbabel cannot be installed")
    #                 print("CHIRIPA: The installation is aborted")
    #                 print("================= ERROR INSTALL ================")
    #                 sys.exit()
    #             else:
    #                 pass
    #
    #         subprocess.call(["rm", "-rf", "thirdparty/openbabel/build"])
    #         subprocess.call(["mkdir", "thirdparty/openbabel/build"])
    #         os.chdir("thirdparty/openbabel/build")
    #         cmake_arguments1 = [ "-DCMAKE_INSTALL_PREFIX={}".format(fullpath_cmake)]
    #         cmake_arguments2 = [ "-DPYTHON_BINDINGS=ON"]
    #         subprocess.check_call(["cmake", "{}", "{}".format(fullpath_cmake)]+cmake_arguments1+cmake_arguments2)
    #         subprocess.call("make")
    #         subprocess.call(["make", "install"])
    #         os.chdir("../../")
    #         subprocess.call(["tar","cvfz", "openbabel.tar.gz", "openbabel"])
    #         subprocess.call(["rm", "-rf", "openbabel ./git"])
    #         os.chdir("..")
    #
    #     print("The *.so library has been installed in: {envdir}/lib/python3.6/site-packages/indigox/pyindigox.cpython-36m-x86_64-linux-gnu.so")
    #     print("                                        {envdir}/lib/python3.6/site-packages/indigox/__init__.py")

    import git

    giturl = 'https://github.com/openbabel/openbabel.git'
    # Look at thirdparty directory
    if not os.path.isdir("thirdparty"):
        os.makedirs("thirdparty")

    install_dir = 'thirdparty'
    eigen_dir = 'thirdparty/eigen/eigen_library/include/eigen3'
    fullpath_cmake = os.path.abspath(install_dir)
    fullpath_eigen = os.path.abspath(eigen_dir)

    if not os.path.isdir(fullpath_cmake+"/openbabel"):
        print("Downloading: ... openbabel-3-1-1(Wait for one minute)")
        #urllib.request.urlretrieve(giturl, fullpath_cmake+"/openbabel-3-1-1.zip")
        git.Repo.clone_from(giturl, fullpath_cmake+"/openbabel")

    try:
        subprocess.check_output(['cmake', '--version'])
    except OSError:
        print("================= ERROR INSTALL ================")
        print("CHIRIPA: Cannot find CMake executable")
        print("CHIRIPA: The installation is aborted")
        print("================= ERROR INSTALL ================")
        sys.exit()

    if not os.path.isdir(fullpath_cmake+"/openbabel/build"):

        subprocess.call(["rm", "-rf", "thirdparty/openbabel/build"])
        subprocess.call(["mkdir", "thirdparty/openbabel/build"])
        os.chdir("thirdparty/openbabel/build")
        cmake_arguments1 = [ "-DCMAKE_INSTALL_PREFIX={}".format(fullpath_cmake+"/openbabel")]
        cmake_arguments2 = [ "-DPYTHON_BINDINGS=ON"]
        cmake_arguments3 = [ "-DRUN_SWIG=ON"]
        cmake_arguments4 = [ "-DEIGEN3_INCLUDE_DIR={}".format(fullpath_eigen)]
        subprocess.check_call(["cmake", "{}", "{}".format(fullpath_cmake+"/openbabel")]+
                              cmake_arguments1+cmake_arguments2+cmake_arguments3+cmake_arguments4)
        subprocess.call(["make", "-j4"])
        subprocess.call(["make", "install"])
        os.chdir("../../")
        os.chdir("..")

    # Copy the library to the root site-engines of the python distribution
    dir_env_python = site.getsitepackages()[0]
    dir_openbabel_installed = fullpath_cmake+"/openbabel/lib/python3.8/site-packages/openbabel"
    cmd = 'cp -rf {} {}'.format(dir_openbabel_installed, dir_env_python)
    print(cmd)
    os.system(cmd)


# Main setup
if __name__ == '__main__':

    # Install requirements ===================================
    with open('requirements.txt') as f:
        required = f.read().splitlines()
    for ipack in required:
        try:
            pkg, version = ipack.split(">=")[0:2]
            install_with_pip(pkg, version)
        except ValueError:
            pkg = ipack
            install_with_pip(pkg)

    # Install third-party software ===========================
    install_indigox_bond()
    install_eigen()
    install_openbabel()

    # Compile external extensions in C =======================
    debug_cflags = False
    use_openmp = True
    has_openmp = detect_openmp()
    if use_openmp and not has_openmp:
        print('No openmp compatible compiler found default to serial build.')

    parallel_libraries = ['gomp'] if has_openmp and use_openmp else []
    mathlib = ['m']
    define_macros = []
    extra_compile_args = ['-std=c99', '-ffast-math', '-O3', '-funroll-loops', '-Wno-cpp']
    if debug_cflags:
        extra_compile_args.extend(['-Wall', '-pedantic'])
        define_macros.extend([('DEBUG', '1')])

    parallel_args = ['-fopenmp'] if has_openmp and use_openmp else []
    parallel_libraries = ['gomp'] if has_openmp and use_openmp else []
    parallel_macros = [('PARALLEL', None)] if has_openmp and use_openmp else []

    # Extensions ==============================================
    print("CHIRIPA: Need to compile external libraries:\n"+"\text_libc/c_distances.pyx\n"
                                                          +"\text_libc/c_distances_openmp.pyx\n"
                                                          +"\text_libc/c_discriminant.pyx\n")
    extensions = [
        Extension("ext_libc.c_distances",["ext_libc/c_distances.pyx"],
                  libraries=mathlib,
                  define_macros=define_macros,
                  extra_compile_args=extra_compile_args,),

        Extension("ext_libc.c_distances_openmp",
                  ["ext_libc/c_distances_openmp.pyx"],

                  libraries=mathlib+parallel_libraries,
                  define_macros=define_macros+parallel_macros,
                  extra_compile_args=parallel_args + extra_compile_args,
                  extra_link_args=parallel_args),

        Extension("ext_libc.c_discriminant", ["ext_libc/c_discriminant.pyx"],
                   libraries=mathlib,
                   define_macros=define_macros,
                   extra_compile_args=extra_compile_args,
                   library_dirs=["ext_libc"],)
        ]

    # Setup Chiripa ===========================================
    setup(
        name = 'Chiripa',
        version = '0.1',
        description='Python library to calculate the chi interaction parameter between monomers.',
        author="Javier Ramos",
        author_email="jrdcasa@gmail.com",
        # This automatically detects the packages in the specified
        # (or current directory if no directory is given).
        packages=find_packages(exclude=['tests', 'examples', 'docs']),
        ext_modules=cythonize(extensions),
    )

