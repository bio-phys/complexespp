from __future__ import absolute_import

from os.path import join as pjoin, exists, isfile
import pytest


class TestDataDir(object):
    """
    Simple class to access a directory with test data
    """

    def __init__(self, folder, data_folder, file_only=False):
        self.folder = pjoin(folder, data_folder)
        self.file_only = file_only

    def __getitem__(self, file):
        data_filename = pjoin(self.folder, file)
        if self.file_only:
            if not isfile(data_filename):
                raise RuntimeError(
                    "no file '{}' found in folder '{}'".format(file, self.folder)
                )

        if exists(data_filename):
            return data_filename
        else:
            raise RuntimeError(
                "no file/folder '{}' found in folder '{}'".format(file, self.folder)
            )


@pytest.fixture
def datafiles(request):
    """access test directory in a pytest. This works independent of where tests are
    started"""
    return TestDataDir(request.fspath.dirname, "data", file_only=True)


@pytest.fixture
def data(request):
    """access test directory in a pytest. This works independent of where tests are
    started"""
    return TestDataDir(request.fspath.dirname, "data")
