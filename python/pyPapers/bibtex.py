#!/usr/bin/python

import bibtexparser

import bibtexparser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.bibdatabase import BibDatabase

class BibDatabaseManager:
    def __init__(self, bibtex_file):
        self.bibtex_file = bibtex_file
        self.entries = self._load_bibtex()

    def _load_bibtex(self):
        with open(self.bibtex_file, 'r') as bib_file:
            bib_database = bibtexparser.load(bib_file)
            return bib_database.entries

    def search(self, field, query):
        results = []
        for entry in self.entries:
            if field in entry and query.lower() in entry[field].lower():
                results.append(entry)
        return results

    def get_full_record(self, entry):
        # Create a BibDatabase object and add the entry
        db = BibDatabase()
        db.entries = [entry]

        # Create a BibTexWriter to convert it back to a BibTeX string
        writer = BibTexWriter()
        return writer.write(db)



# if main
if __name__ == "__main__":
    bib_db = BibDatabaseManager('/home/prokop/Mendeley Desktop/library.bib')
    #search_results = bib_db.search('title', 'machine learning')
    search_results = bib_db.search('title', 'electron force field')

    for result in search_results:
        full_record = bib_db.get_full_record(result)
        print(full_record)

