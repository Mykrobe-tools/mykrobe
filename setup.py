#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
import logging
import os
import requests
import shutil
import subprocess
import tarfile
import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup
from setuptools.command.install import install as DistutilsInstall


class MyInstall(DistutilsInstall):

    def run(self):
        self._install_mccortex()

    def _get_mykrobe_data(self):
        data_tarball_url = "https://ndownloader.figshare.com/files/20996829"
        dir_of_this_file = os.path.dirname(os.path.realpath(__file__))
        mykrobe_dir = os.path.join(dir_of_this_file, "src", "mykrobe")
        assert os.path.exists(mykrobe_dir)
        data_dir = os.path.join(mykrobe_dir, "data")
        if os.path.exists(data_dir):
            shutil.rmtree(data_dir)
        extracted_name = "mykrobe-data"
        tarball_filename = "mykrobe_data.tar.gz"
        request = requests.get(data_tarball_url, allow_redirects=True)
        with open(tarball_filename, 'wb') as t:
            t.write(request.content)
        if os.path.exists(extracted_name):
            shutil.rmtree(extracted_name)
        with tarfile.open(tarball_filename, mode="r") as t:
            t.extractall()
        assert os.path.exists(extracted_name)
        os.rename(extracted_name, data_dir)
        os.unlink(tarball_filename)

    def _install_mccortex(self):
        dir_of_this_file = os.path.dirname(os.path.realpath(__file__))

        # This is for the appveyor testing with tox. Building mccortex required
        # some hacking, so we can't use the simple call to `make` here.
        # Tox runs somewhere else, in an isolated dir, which means it doesn't
        # see the build mccortex31 binary. Use an environment variable to
        # find the built checkout of mccortex, so it doesn't try to build again,
        # which will just fail
        if "TOX_INI_DIR" in os.environ:
            mccortex_git_dir = os.path.join(os.environ["TOX_INI_DIR"], "mccortex")
        else:
            mccortex_git_dir = os.path.join(dir_of_this_file, "mccortex")

        if not os.path.exists(mccortex_git_dir):
            subprocess.call(
                ["git", "clone", "--recursive", "-b", "geno_kmer_count", "https://github.com/Mykrobe-tools/mccortex", mccortex_git_dir], cwd=dir_of_this_file)

        mccortex_build_binary = os.path.join(mccortex_git_dir, "bin", "mccortex31")
        if not os.path.exists(mccortex_build_binary):
            subprocess.call(
                ["make", "clean"], cwd=mccortex_git_dir)
            subprocess.call(
                ["make"], cwd=mccortex_git_dir)

        mccortex_install_dir = os.path.join(dir_of_this_file, "src", "mykrobe", "cortex")
        mccortex_install_binary = os.path.join(mccortex_install_dir, "mccortex31")
        assert os.path.exists(mccortex_install_dir)
        shutil.copy(mccortex_build_binary, mccortex_install_binary)

        DistutilsInstall.run(self)


def read(*names, **kwargs):
    return io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ).read()


setup(
    name='mykrobe',
    version='0.8.2',
    license='MIT',
    description='Antibiotic resistance prediction in minutes',
    # long_description='%s\n%s' % (
    #     re.compile('^.. start-badges.*^.. end-badges',
    #                re.M | re.S).sub('', read('README.md')),
    #     re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst'))
    # ),
    author='Phelim Bradley',
    author_email='wave@phel.im',
    url='https://github.com/Mykrobe-tools/mykrobe',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob(
        'src/*.py')]+[splitext(basename(path))[0] for path in glob('src/*/*.py')]+[splitext(basename(path))[0] for path in glob('src/*/*/*.py')],
    include_package_data=True,
    package_data={"mykrobe": ["cortex/mccortex31"]},
    zip_safe=False,
    classifiers=[
        'Intended Audience :: Developers',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Utilities',
    ],
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    install_requires=[
        'Biopython',
        "pyvcf",
        'mongoengine',
        'requests',
        # eg: 'aspectlib==1.1.1', 'six>=1.7',
    ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    entry_points={
        'console_scripts': [
            'mykrobe = mykrobe.cli:main',
        ]
    },
    cmdclass={'install': MyInstall}
)
