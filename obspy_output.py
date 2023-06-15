from obspy import read, read_events, read_inventory


class ObspyfiedOutput:
    def __init__(self, obspyfied_path: str):
        mseed_file_path, cat_file_path, inv_file_path = self._find_obspyfied_files(obspyfied_path)
        self.mseed_file = read(mseed_file_path)
        self.cat_file = read_events(cat_file_path)
        self.inv_file = read_inventory(inv_file_path)

    
    def _find_obspyfied_files(obspyfied_path: str) -> list:
        pass

    def search_files(directory, keyword, include_subdirectories=True):
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