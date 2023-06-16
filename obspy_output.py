from obspy import read, read_events, read_inventory
import fnmatch
import os 
import sys


class ObspyfiedOutput:
    def __init__(self, obspyfied_path:str = None, mseed_file_path:str = None):
        mseed_file_path, cat_file_path, inv_file_path = self._find_obspyfied_files(obspyfied_path, mseed_file_path)
        self.stream = read(mseed_file_path)
        self.cat = read_events(cat_file_path)
        self.inv = read_inventory(inv_file_path)

        self.mseed_file_name = mseed_file_path.split('/')[-1]
        self.inv_file_name = inv_file_path.split('/')[-1]

    def _find_obspyfied_files(self, obspyfied_path:str = None, mseed_file_path:str = None) -> list:
        if obspyfied_path is not None and mseed_file_path is None:
            mseed_file_path = self._find_mseed_files(obspyfied_path)
            cat_file_path = self._find_cat_files(obspyfied_path)
            inv_file_path = self._find_inv_files(obspyfied_path)
        elif obspyfied_path is None and mseed_file_path is not None:
            obspyfied_path = os.path.dirname(mseed_file_path)
            cat_file_path = self._find_cat_files(obspyfied_path)
            inv_file_path = self._find_inv_files(obspyfied_path)
        return [mseed_file_path, cat_file_path, inv_file_path]

    def _find_inv_files(self, obspyfied_path):
        inv_files = self.search_files(obspyfied_path, 'inv.xml')
        if len(inv_files) > 1:
            raise FileExistsError('Multiple mseed files were found')
            sys.exit(1)
        elif len(inv_files) == 0:
            raise FileNotFoundError('No mseed files were found')
            sys.exit(1)
        else:
            inv_file_path = inv_files[0]
        return inv_file_path

    def _find_cat_files(self, obspyfied_path):
        cat_files = self.search_files(obspyfied_path, 'cat.xml')
        if len(cat_files) > 1:
            raise FileExistsError('Multiple mseed files were found')
            sys.exit(1)
        elif len(cat_files) == 0:
            raise FileNotFoundError('No mseed files were found')
            sys.exit(1)
        else:
            cat_file_path = cat_files[0]
        return cat_file_path


    def _find_mseed_files(self, obspyfied_path):
        mseed_files = self.search_files(obspyfied_path, '.mseed')
        if len(mseed_files) > 1:
            raise FileExistsError('Multiple mseed files were found')
            sys.exit(1)
        elif len(mseed_files) == 0:
            raise FileNotFoundError('No mseed files were found')
            sys.exit(1)
        else:
            mseed_file_path = mseed_files[0]
        return mseed_file_path


    def search_files(self, directory, keyword, include_subdirectories=False):
        """
        Search for files containing a specific keyword in a directory.

        Args:
            directory (str): The directory to search in.
            keyword (str): The specific keyword to search for in file names.
            include_subdirectories (bool, optional): Determines whether to include subdirectories in the search.
                Defaults to True.

        Returns:
            list: A list of file paths that contain the specified keyword.
        """
        matches = []
        for root, dirnames, filenames in os.walk(directory):
            if not include_subdirectories and root != directory:
                break
            for filename in filenames:
                if fnmatch.fnmatch(filename, '*' + keyword + '*'):
                    matches.append(os.path.join(root, filename))
        
        return matches